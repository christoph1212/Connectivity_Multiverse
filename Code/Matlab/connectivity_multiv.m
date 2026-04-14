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
%   dir_Preproc:        String pointing to the folder where Raw Data is
%   dir_Log:            String pointing to where Log-Files should be saved
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
dir_Log     = fullfile(dir_Log, 'Connectivity');

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
parpool("Processes", CONNECTIVITY.nWorkers);

nFiles              = length(Preproc_Files);

%% Morlet-Wavelet Decomposition
freqbands_struct = struct(...
    'delta',    [1 4], ...
    'theta',    [4 8], ...
    'alpha1',   [8 10], ...
    'alpha2',   [10 13], ...
    'beta',     [13 30] ...
    );

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
        ErrorFile = fullfile(dir_Log, ['Error_PreProc_', ErrorFile, '.txt']);
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

    for i_File = 1:2:length(connectivity_files)
        % Get Sub-ID from 1st File
        Filename1 = connectivity_files(i_File).name;
        subID1 = strsplit(Filename1, '_');
        subID1 = strjoin(subID1(1:5), '_');

        % Check if file exists
        if isfile(fullfile(dir_Connect, [subID1 '_full_connectivity.mat']))
            continue
        end
    
        % Get Sub-ID from 2nd File
        Filename2 = connectivity_files(i_File+1).name;
        subID2 = strsplit(Filename2, '_');
        subID2 = strjoin(subID2(1:5), '_');

        if strcmp(subID1, subID2)
            fprintf("Combining File %d/%d\n", i_File, length(connectivity_files)/2)
            file1 = fullfile(connectivity_files(i_File).folder, connectivity_files(i_File).name);
            file2 = fullfile(connectivity_files(i_File+1).folder, connectivity_files(i_File+1).name);

            if contains(file1, 'phase')
                % Load and Combine Data
                connectivity_matrix = load(file1);
                oaec_data = load(file2);
                connectivity_matrix.oaec = oaec_data.oAEC;

                % Save new File as struct
                newOutputFileName = strrep(Filename1, 'phase', 'full');
                OutputFile = fullfile(dir_Connect, newOutputFileName);
                delete(file1, file2)
                save(OutputFile, '-fromstruct', connectivity_matrix)
                fprintf("%s saved - Deleted single files\n", newOutputFileName)

            elseif contains(file2, 'phase')
                connectivity_matrix = load(file2);
                oaec_data = load(file1);
                connectivity_matrix.oaec = oaec_data.oAEC;

                % Save new File as struct
                newOutputFileName = strrep(Filename1, 'oAEC', 'full');
                OutputFile = fullfile(dir_Connect, newOutputFileName);
                delete(file1, file2)
                save(OutputFile, '-fromstruct', connectivity_matrix)
                fprintf("%s saved - Deleted single files\n", newOutputFileName)
            end
            
        else
            warning("No matching files found for: %s", Filename1)
            continue

        end

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
                connectivity_data.oaec.(current_band) = compute_oAEC(EEG, current_band, freqband);
            else
                connectivity_data.imcoh.(current_band) = compute_imcoh(EEG, current_band, freqband);
                connectivity_data.wpli.(current_band) = compute_wpli(EEG, current_band, freqband);
                connectivity_data.pli.(current_band) = compute_pli(EEG, current_band, freqband);
                connectivity_data.pcoh.(current_band) = compute_pcoh(EEG, current_band, freqband);
            end
    
        case 'imcoh'
            connectivity_data.imcoh.(current_band) = compute_imcoh(EEG, current_band, freqband);
    
        case 'wpli'
            connectivity_data.wpli.(current_band) = compute_wpli(EEG, current_band, freqband);
    
        case 'pli'
            connectivity_data.pli.(current_band) = compute_pli(EEG, current_band, freqband);
    
        case 'pcoh'
            connectivity_data.pcoh.(current_band) = compute_pcoh(EEG, current_band, freqband);
    
        case 'oaec'
            connectivity_data.oaec.(current_band) = compute_oAEC(EEG, current_band, freqband);
    
    end % CONNECTIVITY.Measures
end