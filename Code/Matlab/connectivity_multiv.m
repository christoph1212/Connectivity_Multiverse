function connectivity_multiv(dir_Preproc, dir_Log, dir_Connect, CONNECTIVITY, combine_conn_files, Overwrite)
%% Connectivity Analysis
% Run:
%   (1) Morlet Wavelet Decomposition
%   (2) Connectivity Multiverse
%
% Output:
%   (1) .mat file(s) with the unthresholded adjacency matrix (one for phase
%       and one for oAEC based methods)
%
% Inputs:
%   dir_Preproc:        String pointing to the folder where preproc Data is
%   dir_Log:            String pointing to where Log-Files should be saved
%   dir_Connect         String pointing to where connectivity data should
%                       be saved
%   CONNECTIVITY:       Struct with analysis settings
%   combine_conn_files: Boolean. Should connectivity files be combined?
%                       Default: false
%   Overwrite:          Boolean. Should calculation be recalculated? 
%                       Default: false
%
% Created by: Christoph Frühlinger
% Last edited: April 2026

%% get from function input
if nargin < 5
    combine_conn_files = false;
end

if nargin < 6
    Overwrite = false;
end

fprintf(['\n%s\n' ...
         'Starting Connectivity Analysis' ...
         '\n%s\n'], ...
         repmat('=', 1, 100), ...
         repmat('=', 1, 100));

%% Prepare List of Files to be Processed
Preproc_Files   = dir(fullfile(dir_Preproc, '**/*.set'));  

%% Directory where file should be saved
dir_Log         = fullfile(dir_Log, 'Connectivity');

if ~isfolder(dir_Connect)
    mkdir(dir_Connect)
end

if ~isfolder(dir_Log)
    mkdir(dir_Log)
end

fprintf('\nAnalyzing %d Files. \n', length(Preproc_Files));
fprintf(['InputFolder:  %s \n' ...
         'OutputFolder: %s \n' ...
         'LogFolder:    %s \n\n'], ...
         dir_Preproc , dir_Connect, dir_Log);

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
if isempty(CONNECTIVITY.nWorkers)
    parpool("Processes")
else
    % parpool("Processes", CONNECTIVITY.nWorkers);
    poolAvailable = false;

    while ~poolAvailable
        try
            parpool('Processes', CONNECTIVITY.nWorkers);
            poolAvailable = true;
        catch ME
            if contains(ME.message, 'license') || contains(ME.message, 'License')
                fprintf('Parpool currently not available. Waiting 5 Minutes...\n');
                pause(5*60);
            else
                rethrow(ME);
            end
        end
    end
end

nFiles          = length(Preproc_Files);

%% Morlet-Wavelet Decomposition
freqbands_struct = struct(...
    'delta',    [2 4], ...
    'theta',    [4 8], ...
    'alpha1',   [8 10], ...
    'alpha2',   [10 13], ...
    'beta',     [13 30] ...
    );

%% Loop over Files
q = parallel.pool.DataQueue;
afterEach(q, @(msg) fprintf('%s', msg));

parfor i_File = 1:nFiles

    try

        filesplits = split(Preproc_Files(i_File).name, '_');
        OutputFileName = [filesplits{6} '_', strjoin(filesplits(1:5), '_') '_connectivity.mat'];
        OutputFile = fullfile(dir_Connect, OutputFileName);
    
        if isfile(OutputFile) && ~Overwrite
            continue
        end
    
        % Load Preprocessed Data
        InputFile = fullfile(Preproc_Files(i_File).folder, Preproc_Files(i_File).name);
        EEG = pop_loadset(InputFile);
        bands = fieldnames(freqbands_struct);
    
        % Initialise Output Struct
        connectivity_data = struct();
    
        % Get Frequency Band(s) of Interest
        switch CONNECTIVITY.Bands %#ok
    
            case 'all'
            
                for i_band = 1:numel(bands)
                    current_band = bands{i_band};
                    freqband = freqbands_struct.(current_band);
    
                    % Analyze Data using Connectivity Measure(s) of Interest
                    connectivity_data = compute_measures(EEG, CONNECTIVITY.Measures, current_band, freqband, connectivity_data);
    
                end % i_band 
    
            otherwise
                band = CONNECTIVITY.Bands;
                if ~isfield(freqbands_struct, band)
                    error(['Unknown frequency band: %s. Possible options:' ...
                        'delta, theta, alpha1, alpha2, beta, or all'])
                end
                current_band = band;
                freqband = freqbands_struct.(current_band);
    
                % Analyze Data using Connectivity Measure(s) of Interest
                connectivity_data = compute_measures(EEG, CONNECTIVITY.Measures, current_band, freqband, connectivity_data);
                
        end % CONNECTIVITY.Bands
        
        % Save File as struct
        save(OutputFile, '-fromstruct', connectivity_data)
        send(q, sprintf('[%d/%d] %s saved\n', i_File, nFiles, OutputFileName));

    catch e
        ErrorMessage = string(e.message);
        for ierrors = 1:length(e.stack)
            ErrorMessage = strcat(ErrorMessage, " // ", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
        end
    
        fprintf('***Error in File: %s;\n%s.\n', InputFile, ErrorMessage);
    
        % make error log
        [~, ErrorFile, ext] = fileparts(InputFile);
        ErrorFile = [ErrorFile ext];
        fprintf('Problem executing File: %s\n',ErrorFile);
        fprintf('The Error Message is: \n%s \n',ErrorMessage);
        [~, ErrorFile, ~] = fileparts(ErrorFile);
        ErrorFile = fullfile(dir_Log, ['Error_Connect_', ErrorFile, '.txt']);
        fid1 = fopen(ErrorFile, 'wt');
        fprintf(fid1, 'Error-Subject: %s \nThe returned Error Message is: \n\n%s \n', InputFile,  ErrorMessage);
        fclose(fid1);
    
    end % try-catch

end % i_File

% Combine Connectivity Files into one File
if combine_conn_files && strcmp(CONNECTIVITY.Measures, 'all')
    fprintf("\n\nConnectivity Analysis completed.\n\nCombining files. Please wait...\n")

    connectivity_files = dir(fullfile(dir_Connect, '*.mat'));
    [~, sort_idx] = sort({connectivity_files.name});
    connectivity_files = connectivity_files(sort_idx);

    phase_files = connectivity_files(contains({connectivity_files.name}, 'phase'));

    for i_File = 1:length(phase_files)
        % Get Sub-ID from 1st File
        Filename_phase = phase_files(i_File).name;
        parts = strsplit(Filename_phase, '_');
        subID = strjoin(parts(1:5), '_');
         
        OutputFilename = [subID '_full_connectivity.mat'];
        % Check if file exists
        if isfile(fullfile(dir_Connect, OutputFilename))
            continue
        end
    
        oaec_pattern = [strjoin(parts(1:5), '_') '_oAEC_*'];
        oaec_match   = dir(fullfile(dir_Connect, oaec_pattern));

        if isempty(oaec_match)
            warning("No matching oAEC file found for: %s", Filename_phase)
            continue
        elseif length(oaec_match) > 1
            warning("Multiple oAEC matches for %s – skipping.", subID)
            continue
        end

        fprintf("Combining file %d/%d: %s\n", i_File, length(phase_files), subID)

        file_phase = fullfile(phase_files(i_File).folder, phase_files(i_File).name);
        file_oaec = fullfile(oaec_match(1).folder, oaec_match(1).name);

        connectivity_matrix       = load(file_phase);
        oaec_data                 = load(file_oaec);
        connectivity_matrix.oaec  = oaec_data.oaec;

        OutputFile = fullfile(dir_Connect, OutputFilename);
        save(OutputFile, '-fromstruct', connectivity_matrix)
        delete(file_phase, file_oaec)
        fprintf("%s saved - Deleted single files\n", OutputFilename)

    end
    fprintf("\nCombining files completed.\n")
end
fprintf(['\n%s\n' ...
        'Finished Connectivity Analysis' ...
        '\n%s\n'], ...
        repmat('=', 1, 100), ...
        repmat('=', 1, 100));
end

function connectivity_data = compute_measures(EEG, measures, current_band, freqband, connectivity_data)
    switch measures
    
        case 'all'
            % oAEC and Phase-based measure are stored in different files (different epoch lengths)
            if contains(EEG.filename, 'oAEC')
                connectivity_data.oaec.(current_band).unthresh = compute_oAEC(EEG, current_band, freqband);
            else
                AS = compute_analytic_signal(EEG, current_band, freqband);
                connectivity_data.imcoh.(current_band).unthresh = compute_imcoh(EEG, current_band, freqband, AS);
                connectivity_data.wpli.(current_band).unthresh = compute_wpli(EEG, current_band, freqband, AS);
                connectivity_data.pli.(current_band).unthresh = compute_pli(EEG, current_band, freqband, AS);
                connectivity_data.pcoh.(current_band).unthresh = compute_pcoh(EEG, current_band, freqband);
            end
    
        case 'imcoh'
            connectivity_data.imcoh.(current_band).unthresh = compute_imcoh(EEG, current_band, freqband);
    
        case 'wpli'
            connectivity_data.wpli.(current_band).unthresh = compute_wpli(EEG, current_band, freqband);
    
        case 'pli'
            connectivity_data.pli.(current_band).unthresh = compute_pli(EEG, current_band, freqband);
    
        case 'pcoh'
            connectivity_data.pcoh.(current_band).unthresh = compute_pcoh(EEG, current_band, freqband);
    
        case 'oaec'
            connectivity_data.oaec.(current_band).unthresh = compute_oAEC(EEG, current_band, freqband);
    
    end % CONNECTIVITY.Measures
end