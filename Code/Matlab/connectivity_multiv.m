function connectivity_multiv(dir_Preproc, dir_Log, dir_Connect, CONNECTIVITY, Overwrite)
%% Connectivity Analysis
% Run:
%   (1) Morlet Wavelet Decomposition
%   (2) Connectivity Multiverse
%
% Output:
%   (1) ...
%
% Inputs:
%   dir_Preproc:  String pointing to the folder where Raw Data is
%   dir_Log:      String pointing to where Log-Files should be saved
%   CONNECTIVITY: Struct with analysis settings
%   Overwrite:    Boolean. Should calculation be recalculated? 
%                 Default: false
%
% Created by: Christoph Frühlinger
% Last edited: March 2026

%% get from function input
if nargin < 5
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
% LogFilename = fullfile(dir_Log, 'All_Logs.csv');
% dir_Log     = fullfile(dir_Log, 'Connectivity');

if ~isfolder(dir_Connect)
    mkdir(dir_Connect)
end

% if ~isfolder(dir_Log)
%     mkdir(dir_Log)
% end

fprintf('\n%d Files found. \n', length(Preproc_Files));
fprintf('InputFolder is %s \nOutputFolder is %s \nLogFolder is %s \n\n', ...
    dir_Preproc , dir_Connect, dir_Log);

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
parpool("Processes", CONNECTIVITY.nWorkers);

Files_Connect       = dir(dir_Connect);
fileNames_Connect   = {Files_Connect.name};
nFiles              = length(Preproc_Files);

%% Morlet-Wavelet Decomposition
freqbands_struct = struct(...
    'delta',    [1 4], ...
    'theta',    [4 8], ...
    'alpha1',   [8 10], ...
    'alpha2',   [10 13], ...
    'beta',     [13 30] ...
    );

parfor i_File = 1:nFiles

    filesplits = split(Preproc_Files(i_File).name, '_');
    OutputFileName = [filesplits{6} '_', strjoin(filesplits(1:5), '_') '_connectivity.mat'];

    if ismember(OutputFileName, fileNames_Connect) && ~Overwrite
        continue
    end

    % Load Preprocessed Data
    InputFile = fullfile(Preproc_Files(i_File).folder, Preproc_Files(i_File).name);
    EEG = pop_loadset(InputFile);
    bands = fieldnames(freqbands_struct);

    % Initialise Output Struct
    connectivity_data = struct();

    switch CONNECTIVITY.Bands %#ok

        case 'all'
        
            for i_band = 1:numel(bands)
                current_band = bands{i_band};
                freqband = freqbands_struct.(current_band);

                switch CONNECTIVITY.Measures

                    case 'all'
                        if contains(EEG.filename, 'oAEC')
                            connectivity_data.oAEC.(current_band) = compute_oAEC(EEG, current_band, freqband);
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
                        connectivity_data.oAEC.(current_band) = compute_oAEC(EEG, current_band, freqband);

                end % CONNECTIVITY.Measures

            end % i_band 

        otherwise
            band = CONNECTIVITY.Bands;
            if ~isfield(freqbands_struct, band)
                error(['Unknown frequency band: %s. Possible options:' ...
                    'delta, theta, alpha1, alpha2, beta, or all'])
            end
            current_band = band;
            freqband = freqbands_struct.(current_band);

            switch CONNECTIVITY.Measures

                case 'all'
                    if contains(EEG.filename, 'oAEC')
                        connectivity_data.oAEC.(current_band) = compute_oAEC(EEG, current_band, freqband);
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
                    connectivity_data.oAEC.(current_band) = compute_oAEC(EEG, current_band, freqband);

            end % CONNECTIVITY.Measures
            
    end % CONNECTIVITY.Bands
    
    OutputFile = fullfile(dir_Connect, OutputFileName);
    save(OutputFile, '-fromstruct', connectivity_data)
    fprintf("%s saved\n", OutputFileName)

end % i_File

end