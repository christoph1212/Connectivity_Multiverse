function graph_metrics(dir_Connect, dir_Log, GRAPH, Overwrite)
%% Connectivity Analysis - Graph Theory
% Run:
%   (1) Calculate graph-theoretic metrics
%
% Output:
%   (1) updated .mat file with the thresholded adjacency matrix 
%
% Inputs:
%   dir_Connect:        String pointing to the folder where connectivity
%                       Data is stored (also Output Folder)
%   dir_Log:            String pointing to where Log-Files should be saved
%   GRAPH:              Struct with analysis settings
%   Overwrite:          Boolean. Should calculation be recalculated? 
%                       Default: false
%
%
% Created by: Christoph Frühlinger
% Last edited: May 2026

%% get from function input
if nargin < 4    
    Overwrite = false;
end

fprintf(['\n%s\n' ...
         'Calculating Graph Metrics' ...
         '\n%s\n'], ...
         repmat('=', 1, 100), ...
         repmat('=', 1, 100));

%% Prepare List of Files to be Processed
Connect_Files   = dir(fullfile(dir_Connect, '*.mat'));  
nFiles          = length(Connect_Files);

%% Directory where file should be saved
dir_Log         = fullfile(dir_Log, 'Connectivity');

if ~isfolder(dir_Log)
    mkdir(dir_Log)
end

addpath("SmallWorldNess-master");

fprintf('\nCalculating Metrics for %d Files. \n', nFiles);
fprintf(['InputFolder:  %s \n' ...
         'OutputFolder: %s \n' ...
         'LogFolder:    %s \n\n'], ...
         dir_Connect , dir_Connect, dir_Log);

OutputFile = fullfile(dir_Connect, 'Graph_metrics.csv');

% Check if file exists
if isfile(OutputFile) && ~Overwrite
    fprintf('Graph metrics already calculated. Skipping.\n');
    return
end

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
if isempty(GRAPH.nWorkers)
    parpool("Processes")
else
    parpool("Processes", GRAPH.nWorkers);
end

%% Loop over Files
q = parallel.pool.DataQueue;
afterEach(q, @(msg) fprintf('%s', msg));

all_results = cell(nFiles, 1);

parfor i_File = 1:nFiles

    try

        % Load Connectivity Data        
        InputFile = fullfile(Connect_Files(i_File).folder, Connect_Files(i_File).name);
        connectivity_data = load(InputFile);        

        % Get info from connectivity data
        switch GRAPH.metrics %#ok
            case 'all'
                metrics = {'cc', 'pathl', 'eglob', 'eloc', 'smallworld'}';
            otherwise
                metrics = {GRAPH.metrics};
        end
        measures = fieldnames(connectivity_data);
        bands    = fieldnames(connectivity_data.(measures{1}));
        thresh   = fieldnames(connectivity_data.(measures{1}).(bands{1}));
        thresh   = thresh(~strcmp('unthresh', thresh));

        % Create empty table
        headers = {'ID', 'Run', 'Condition'};
        for i_meas = 1:numel(measures)
            for i_band = 1:numel(bands) 
                headers{end+1} = sprintf('Percol_Thresh_%s_%s', measures{i_meas}, bands{i_band}); %#ok
            end
        end

        for i_metric = 1:numel(metrics)
            for i_meas = 1:numel(measures)
                for i_band = 1:numel(bands)                    
                    for i_thresh = 1:numel(thresh)
                        headers{end+1} = sprintf('%s_%s_%s_%s', metrics{i_metric}, measures{i_meas}, bands{i_band}, thresh{i_thresh}); %#ok
                    end
                end
            end
        end        
                
        splits = strsplit(Connect_Files(i_File).name, '_');
        results = table();
        results.ID        = {splits{1}}; %#ok
        results.Run       = {splits{2}}; %#ok
        results.Condition = {[splits{4} '_' splits{5}]};
        for i_h = 4:numel(headers)
            results.(headers{i_h}) = NaN;
        end

        % Loop through the whole structure file
        for i_meas = 1:numel(measures)

            for i_band = 1:numel(bands)

                for i_thresh = 1:numel(thresh)

                    am = connectivity_data.(measures{i_meas}).(bands{i_band}).(thresh{i_thresh});
                    is_dens = strcmp(thresh{i_thresh}, 'dens');

                    switch GRAPH.metrics

                        case 'all'
                    
                            if is_dens
                                [AUC_CC, AUC_PL, AUC_Eglob, AUC_Eloc, AUC_SW, ~, percol_thresh] = density_metrics_auc(am, GRAPH);

                                header = sprintf('cc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_CC;

                                header = sprintf('pathl_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_PL;

                                header = sprintf('eglob_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_Eglob;
                                
                                header = sprintf('eloc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_Eloc;

                                header = sprintf('smallworld_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_SW;

                                header = sprintf('Percol_Thresh_%s_%s', measures{i_meas}, bands{i_band});
                                results.(header)(1) = percol_thresh;
                                
                            else
                                
                                % Small World Index
                                header = sprintf('smallworld_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                n = size(am,1);
                                k = sum(am);
                                m = sum(k)/2;
                                Num_ER_repeats = 100;
                                FLAG_Cws = 1;
                                
                                [Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
                                Lrand_mean = mean(Lrand(Lrand < inf));
                                [S_ws_MC,C,L] = small_world_ness(am,Lrand_mean,mean(CrandWS),FLAG_Cws);

                                results.(header)(1) = S_ws_MC;
                                
                                % Clustering Coefficient
                                header = sprintf('cc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = C;

                                % Characteristic Path Length
                                header = sprintf('pathl_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});                                
                                results.(header)(1) = L;

                                % Global Efficiency
                                header = sprintf('eglob_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = efficiency_bin(am);

                                % Local Efficiency
                                header = sprintf('eloc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = mean(efficiency_bin(am,1));
                        
                            end

                        case 'cc'
                            if is_dens
                                [~, ~, ~, ~, ~, AUC_metric, percol_thresh] = density_metrics_auc(am, GRAPH);

                                header = sprintf('cc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_metric;

                                header = sprintf('Percol_Thresh_%s_%s', measures{i_meas}, bands{i_band});
                                results.(header)(1) = percol_thresh;

                            else
                                header = sprintf('cc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = mean(clustering_coef_bu(am));
                            end

                        case 'pathl'
                            if is_dens
                                [~, ~, ~, ~, ~, AUC_metric, percol_thresh] = density_metrics_auc(am, GRAPH);

                                header = sprintf('pathl_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_metric;
                                
                                header = sprintf('Percol_Thresh_%s_%s', measures{i_meas}, bands{i_band});
                                results.(header)(1) = percol_thresh;

                            else
                                header = sprintf('pathl_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                [~,D] = reachdist(am);
                                results.(header)(1) = mean(D(:));
                            end

                        case 'eglob'
                            if is_dens
                                [~, ~, ~, ~, ~, AUC_metric, percol_thresh] = density_metrics_auc(am, GRAPH);

                                header = sprintf('eglob_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_metric;
                                
                                header = sprintf('Percol_Thresh_%s_%s', measures{i_meas}, bands{i_band});
                                results.(header)(1) = percol_thresh;

                            else
                                header = sprintf('eglob_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = efficiency_bin(am);
                            end

                        case 'eloc'
                            if is_dens
                                [~, ~, ~, ~, ~, AUC_metric, percol_thresh] = density_metrics_auc(am, GRAPH);
                                header = sprintf('eloc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_metric;
                                
                                header = sprintf('Percol_Thresh_%s_%s', measures{i_meas}, bands{i_band});
                                results.(header)(1) = percol_thresh;
                            else
                                header = sprintf('eloc_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = mean(efficiency_bin(am,1));
                            end

                        case 'smallworld'
                            if is_dens
                                [~, ~, ~, ~, ~, AUC_metric, percol_thresh] = density_metrics_auc(am, GRAPH);
                                header = sprintf('smallworld_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                results.(header)(1) = AUC_metric;
                                
                                header = sprintf('Percol_Thresh_%s_%s', measures{i_meas}, bands{i_band});
                                results.(header)(1) = percol_thresh;
                            else
                                 header = sprintf('smallworld_%s_%s_%s', measures{i_meas}, bands{i_band}, thresh{i_thresh});
                                n = size(am,1);
                                k = sum(am);
                                m = sum(k)/2;
                                Num_ER_repeats = 100;
                                FLAG_Cws = 1;
                                
                                [Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
                                Lrand_mean = mean(Lrand(Lrand < inf));
                                [S_ws_MC,~,~] = small_world_ness(am,Lrand_mean,mean(CrandWS),FLAG_Cws);

                                results.(header)(1) = S_ws_MC;
                            end

                    end % switch

                end % i_thresh
            end % i_band
        end % i_meas           
                
        
        all_results{i_File} = results;        
        send(q, sprintf('[%d/%d] %s done\n', i_File, nFiles, Connect_Files(i_File).name));

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
        ErrorFile = fullfile(dir_Log, ['Error_Graph_', ErrorFile, '.txt']);
        fid1 = fopen(ErrorFile, 'wt');
        fprintf(fid1, 'Error-Subject: %s \nThe returned Error Message is: \n\n%s \n', InputFile,  ErrorMessage);
        fclose(fid1);

    end % try
end % for i_File

all_results = vertcat(all_results{:});

writetable(all_results, OutputFile);

fprintf("\nGraph Metrics saved to %s\n", OutputFile)

fprintf(['\n%s\n' ...
         'Finished Calculating Graph Metrics' ...
         '\n%s\n'], ...
         repmat('=', 1, 100), ...
         repmat('=', 1, 100)); 
end