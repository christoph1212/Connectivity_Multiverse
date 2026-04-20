function thresholding_multiv(dir_Connect, dir_Log, THRESHOLD, Overwrite)
%% Connectivity Analysis - Thresholding Multiverse
% Run:
%   (1) Thresholding Multiverse
%
% Output:
%   (1) updated .mat file with the thresholded adjacency matrix 
%
% Inputs:
%   dir_Connect:        String pointing to the folder where unthresholded
%                       Data is stored and where data will be updated
%   dir_Log:            String pointing to where Log-Files should be saved
%   THRESHOLD:          Struct with analysis settings
%   Overwrite:          Boolean. Should calculation be recalculated? 
%                       Default: false
%
% Created by: Christoph Frühlinger
% Last edited: April 2026

%% get from function input
if nargin < 4    
    Overwrite = false;
end

fprintf(['\n%s\n' ...
         'Starting Thresholding' ...
         '\n%s\n'], ...
         repmat('=', 1, 100), ...
         repmat('=', 1, 100));

%% Prepare List of Files to be Processed
Connect_Files   = dir(fullfile(dir_Connect, '*.mat'));  

%% Directory where file should be saved
dir_Log         = fullfile(dir_Log, 'Connectivity');

if ~isfolder(dir_Log)
    mkdir(dir_Log)
end

fprintf('\nThresholding %d Files. \n', length(Connect_Files));
fprintf(['InputFolder:  %s \n' ...
         'OutputFolder: %s \n' ...
         'LogFolder:    %s \n\n'], ...
         dir_Connect , dir_Connect, dir_Log);

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
if isempty(THRESHOLD.nWorkers)
    parpool("Processes")
else
    parpool("Processes", THRESHOLD.nWorkers);
end

nFiles          = length(Connect_Files);

%% Loop over Files
q = parallel.pool.DataQueue;
afterEach(q, @(msg) fprintf('%s', msg));

parfor i_File = 1:nFiles

    try

        % Load Connectivity Data
        InputFile = fullfile(Connect_Files(i_File).folder, Connect_Files(i_File).name);
        connectivity_data = load(InputFile);
        
        % Check if Data has been thresholded before
        measures = fieldnames(connectivity_data);
        first_measure = measures{1};
        bands = fieldnames(connectivity_data.(first_measure));
        first_band = bands{1};

        if numel(fieldnames(connectivity_data.(first_measure).(first_band))) > 1 && ~Overwrite
            continue
        end      

        for i_measure = 1:numel(measures)
            current_measure = measures{i_measure};
            bands = fieldnames(connectivity_data.(current_measure));

            for i_band = 1:numel(bands)
                current_band = bands{i_band};
                
                % Get unthresholded adjacency matrix
                am = connectivity_data.(current_measure).(current_band).untresh;

                switch THRESHOLD.Method %#ok

                    case 'all'

                        % Density-based Thresholding
                        connectivity_data.(current_measure).(current_band).dens = density_thresholding(am);
        
                        % orthogonalized Minimum Spanning Tree Thresholding
                        omst_path = fullfile(dir_Root, 'Code', 'Matlab', 'OMST');
                        addpath(omst_path);
                        plot_gce = 0;                       
                        [~, CIJtree, ~,  ~, ~, ~] = threshold_omst_gce_wu(am, plot_gce);
                        connectivity_data.(current_measure).(current_band).omst = double(CIJtree > 0);
        
                        % Efficiency Cost Optimization Thresholding
                        eco_path = fullfile(dir_Root, 'Code', 'Matlab', 'ECO');
                        addpath(eco_path);
                        directed = 0;
                        connectivity_data.(current_measure).(current_band).eco = ECOfilter(am, directed);
        
                        % Minimum Connected Component Thresholding
                        connectivity_data.(current_measure).(current_band).mcc = mcc_thresholding(am);

                    case 'dens'

                        % Density-based Thresholding
                        connectivity_data.(current_measure).(current_band).dens = density_thresholding(am);

                    case 'omst'

                        % orthogonalized Minimum Spanning Tree Thresholding
                        omst_path = fullfile(dir_Root, 'Code', 'Matlab', 'OMST');
                        addpath(omst_path);
                        plot_gce = 0;                       
                        [~, CIJtree, ~,  ~, ~, ~] = threshold_omst_gce_wu(am, plot_gce);
                        connectivity_data.(current_measure).(current_band).omst = double(CIJtree > 0);

                    case 'eco'

                        % Efficiency Cost Optimization Thresholding
                        eco_path = fullfile(dir_Root, 'Code', 'Matlab', 'ECO');
                        addpath(eco_path);
                        directed = 0;
                        connectivity_data.(current_measure).(current_band).eco = ECOfilter(am, directed);

                    case 'mcc'

                        % Minimum Connected Component Thresholding
                        connectivity_data.(current_measure).(current_band).mcc = mcc_thresholding(am);

                end % THRESHOLD.Method 
            end % i_band
        end % i_measure

        % Save File as struct
        OutputFile = InputFile;
        save(OutputFile, '-fromstruct', connectivity_data)
        send(q, sprintf('[%d/%d] %s saved\n', i_File, nFiles, InputFile));

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
        ErrorFile = fullfile(dir_Log, ['Error_Thresh_', ErrorFile, '.txt']);
        fid1 = fopen(ErrorFile, 'wt');
        fprintf(fid1, 'Error-Subject: %s \nThe returned Error Message is: \n\n%s \n', InputFile,  ErrorMessage);
        fclose(fid1);
    
    end % try-catch

end % i_File
    
