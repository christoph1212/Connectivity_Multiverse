function preprocess_data(dir_Raw, dir_Log, dir_Preproc, PREPROC, Overwrite)
%% Preprocessing EEG files
% Run preprocessing 
%
% Output:
%   (1) a .set file of the preprocessed data, one per subject and condition
%       and run (six files per subject)
%   (2) a .csv Log-file of the preprocessing
%   (3) a potential .txt Error-file
%
% Inputs:
%   dir_Raw:    String pointing to the folder where Raw Data is
%   dir_Root:   String pointing to the project's parent folder
%   dir_Log:    String pointing to where Log-Files should be saved
%   PREPROC:    Struct with preprocessing settings
%   Overwrite:  Boolean. Should calculation be recalculated? 
%               Default: false
%
% Created by: Christoph Frühlinger
% Last edited: April 2026

%% get from function input
if nargin < 5
    Overwrite = false;
end

fprintf(['\n%s\n' ...
         'Starting Preprocessing' ...
         '\n%s\n'], ...
         repmat('=', 1, 100), ...
         repmat('=', 1, 100));

%% Prepare List of Files to be Processed
% get *.set files
Raw_Files   = dir(fullfile(dir_Raw, '**/*.set'));

%% Directory where file should be saved
LogFilename = fullfile(dir_Log, 'All_Logs.csv');
dir_Log     = fullfile(dir_Log, 'Preproc');

if ~isfolder(dir_Preproc)
    mkdir(dir_Preproc)
end

if ~isfolder(dir_Log)
    mkdir(dir_Log)
end

fprintf('\nPreprocessing %d Files. \n', length(Raw_Files));
fprintf(['InputFolder:  %s \n' ...
         'OutputFolder: %s \n' ...
         'LogFolder:    %s \n\n'], ...
         dir_Raw , dir_Preproc, dir_Log);

%% Infos on Triggers and Conditions
SplitStruct = struct('Trigger', {11, 12 21 22 31 32}, ...
    'Condition', { ...
    'first_run_eyes_open' 'first_run_eyes_closed' ...
    'second_run_eyes_open' 'second_run_eyes_closed' ...
    'third_run_eyes_open' 'third_run_eyes_closed' ...
    });

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
if isempty(PREPROC.nWorkers)
    parpool("Processes")
else
    % parpool("Processes", PREPROC.nWorkers);
    poolAvailable = false;

    while ~poolAvailable
        try
            parpool('Processes', PREPROC.nWorkers);
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

FileName            = '';
InputFile           = '';
Cond_FileName       = '';

Files_PreProc       = dir(dir_Preproc);
fileNames_PreProc   = {Files_PreProc.name};
nSubs               = length(Raw_Files);

%% Looped preprocessing
q = parallel.pool.DataQueue;
afterEach(q, @(msg) fprintf('%s', msg));

parfor i_Sub = 1:nSubs

    % Check if Subject has been preprocessed 
    if sum(contains(fileNames_PreProc,Raw_Files(i_Sub).name)) == 4 && ...
            ~Overwrite 
        continue
    end
    % ignore subs with unknown reference
    if contains(Raw_Files(i_Sub).name, "sub-SS08EL29") || ...
        contains(Raw_Files(i_Sub).name, "sub-EL07EL18")
        continue
    end

    run_preproc_silent(Raw_Files(i_Sub), dir_Preproc, dir_Log, SplitStruct, ...
        FileName, InputFile, Cond_FileName, PREPROC, Overwrite)

    send(q, sprintf('\n%d/%d preprocessed\n', i_Sub, nSubs));

end

fprintf('\nCompiling Log Files...\n')

% Set path and load files
log_files = dir(fullfile(dir_Log, 'Log_*.csv'));

% Initialize empty table
log_cell = cell(length(log_files), 1);
for i = 1:length(log_files)
    fprintf('Reading Log file %d/%d\n', i, length(log_files))
    file = fullfile(dir_Log, log_files(i).name);
    log_cell{i} = readtable(file, VariableNamingRule="preserve");
end
all_logs = vertcat(log_cell{:});

% Adapt Table
all_logs.ID = string(all_logs.ID);
all_logs.Run = string(all_logs.Run);
all_logs.Condition = string(all_logs.Condition);

% Save Table as csv File
writetable(all_logs, LogFilename);
fprintf('Compiled Log File saved to %s\n', LogFilename);

fprintf(['\n%s\n' ...
     'Finished Preprocessing' ...
     '\n%s\n'], ...
     repmat('=', 1, 100), ...
     repmat('=', 1, 100));
end


function run_preproc_silent(Raw_File, dir_Preproc, dir_Log, SplitStruct, ...
    FileName, InputFile, Cond_FileName, PREPROC, Overwrite)

try
    %% Step 1: Import Data
    InputFile = fullfile(Raw_File.folder, Raw_File.name);

    fprintf("Standardizing File: %s\n", Raw_File.name)

    evalc("EEG = pop_loadset(InputFile);");    

    %% Step 2: remove unwanted channels
    % AFZ and FPZ are used as grounds in some labs
    Common_Channels =  {'FP1', 'FP2', 'AF7', 'AF8', 'AF3', 'AF4', ...
        'F1','F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'FT7', 'FT8', ...
        'FC5', 'FC6', 'FC3', 'FC4', 'FC1', 'FC2', 'C1', 'C2', 'C3', ...
        'C4', 'C5', 'C6', 'T7', 'T8', 'TP7', 'TP8', 'CP5', 'CP6', ...
        'CP3', 'CP4', 'CP1', 'CP2', 'P1', 'P2', 'P3', 'P4', 'P5', ...
        'P6', 'P7', 'P8', 'PO7', 'PO8', 'PO3', 'PO4', 'O1', 'O2', ...
        'OZ', 'POZ', 'PZ', 'CPZ', 'CZ', 'FCZ', 'FZ', 'VOGabove', ...
        'VOGbelow', 'HOGl', 'HOGr'};
    
    Common_Channels = Common_Channels(ismember(Common_Channels, ...
        {EEG.chanlocs.labels})); %#ok
    evalc("EEG = pop_select(EEG, 'channel', Common_Channels);");
    
    %% Step 3: rereference all Files to FCZ (differs between Labs)
    if strcmp('CMS/DRL', EEG.Info_Lab.Reference)
        evalc("EEG = pop_reref(EEG, 'FCZ');");
        
    elseif strcmp('CZ', EEG.Info_Lab.Reference)
        evalc("EEG = pop_reref(EEG, 'FCZ', 'refloc', struct('labels',{'CZ'}, 'type',{'EEG'}, 'ref', [], 'urchan', [], 'theta',{0}, 'radius',{0},'X',{6.12e-17},'Y',{0},'Z',{1},'sph_theta', {0},'sph_phi',{90},'sph_radius',{1}));");
        
    elseif strcmp('FCZ', EEG.Info_Lab.Reference)
        % do nothing

    else
        msg = ['An error occured while rereferencing to FCz. ' ...
            'EEG.Info_Lab.Reference does not contain CMS/DRL, ' ...
            'Cz or FCz.'];
        error(msg)
    end

    % append new reference channel 'FCz' to chanlocs (as no data channel)
    RefInfo = {EEG.nbchan 'labels' 'FCZ' 'theta' 0 'radius' 0.127 ...
        'X' 0.388 'Y' 0 'Z' 0.922 'sph_theta' 0 'sph_phi' 67.2 ...
        'sph_radius' 1 'type' 'EEG' 'datachan' 0 }; %#ok
    evalc("EEG = pop_chanedit(EEG, 'append', RefInfo);");
    

    % Chanlocs without Ref - save chanlocs for interpolation:
    evalc("EEG_interp = pop_select( EEG, 'chantype','EEG');");

    % some datasets miss OZ 
    if length(EEG_interp.chanlocs) ~= 58 %#ok
        if length(EEG_interp.chanlocs) == 57 && ...
            ~ismember('OZ', upper({EEG_interp.chanlocs.labels}))
            % Add OZ to chanlocs interpolation template
            Oz.labels     = 'OZ';
            Oz.theta      = 180;
            Oz.radius     = 0.5000;
            Oz.X          = -1;
            Oz.Y          = -1.2246e-16;
            Oz.Z          = 0;
            Oz.sph_theta  = -180;
            Oz.sph_phi    = 0;
            Oz.sph_radius = 1;
            Oz.type       = 'EEG';
            Oz.ref        = '';
            Oz.urchan     = [];

            EEG_interp.chanlocs(end+1) = Oz;
            EEG_interp.nbchan = EEG_interp.nbchan + 1;
        else
            msg = ['Less than 57 channels. Number of channels in dataset is ' ...
                num2str(length(EEG_interp.chanlocs)) '.'];
            error(msg);
        end
    end
    
    %% Step 4: Downsampling
    % Different Sampling Rate per Site
    if PREPROC.Downsample
        new_srate = EEG.srate/2; %#ok
        evalc("EEG = pop_resample(EEG, new_srate);");
    end

    %% Step 5: Filtering
    % High-Pass Filter
    evalc("EEG = pop_eegfiltnew(EEG, 'locutoff', PREPROC.HP_Filter);");

    % Low-Pass Filter
    evalc("EEG = pop_eegfiltnew(EEG, 'hicutoff', PREPROC.LP_Filter);");
            
    %% Step 6: Separate into Conditions
    EEG_Complete = EEG; %#ok

    % Select triggers and conditions
    if contains(Raw_File.name, 'run-1')
        Rel_SplitStruct.Trigger = [SplitStruct(1:2).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(1:2).Condition};

    elseif contains(Raw_File.name, 'run-2')
        Rel_SplitStruct.Trigger = [SplitStruct(3:4).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(3:4).Condition};

    elseif contains(Raw_File.name, 'run-3')
        Rel_SplitStruct.Trigger = [SplitStruct(5:6).Trigger];
        Rel_SplitStruct.Condition = {SplitStruct(5:6).Condition};
    end

    % Preprocess within each file
    for i_cond = 1:length(Rel_SplitStruct.Trigger)

        try

            Cond_FileName = Rel_SplitStruct.Condition{i_cond};
            FileName = [Cond_FileName,'_',Raw_File.name];
            % oAEC Files
            FileName_oAEC = [Rel_SplitStruct.Condition{i_cond},'_oAEC_',Raw_File.name];
            % Phase-based Files
            FileName_phase = [Rel_SplitStruct.Condition{i_cond},'_phase_',Raw_File.name];

            fprintf("Preprocessing File: %s, Condition: %s\n", Raw_File.name, Cond_FileName);

            % check if files exist already
            if isfile(fullfile(dir_Preproc, FileName_oAEC)) && ...
                    isfile(fullfile(dir_Preproc, FileName_phase)) && ...
                    ~Overwrite
                fprintf('Previously finished Subject %s.\n', Raw_File.name);
                continue
            end
            
            % Check if file is not empty
            if EEG.pnts < EEG.srate
                fprintf('Dataset empty (%s: Condition: %s)\n', Raw_File.name, Cond_FileName);
                
                MissingDataFile = fullfile(dir_Log, ['Error_NotEnoughData_for_Preproc_', FileName, '.txt']);
                fid1 = fopen( MissingDataFile, 'wt' );
                fprintf(fid1, 'Missing Data-File= %s \n Data only includes %i Points. \n', FileName,  EEG.pnts);
                fclose(fid1);
                continue                
                
            end            

            % Select only condition-specific data
            evalc("EEG = pop_epoch( EEG_Complete, {num2str( Rel_SplitStruct.Trigger(i_cond))}, [0  60], 'epochinfo', 'yes');");
            evalc("EEG = eeg_epoch2continuous(EEG);");

            %% Step 7: Bad Channels
           
            % Define Settings based on PREPROC
            if PREPROC.BadChans
                args = {...
                    'FlatlineCriterion',    5, ...
                    'ChannelCriterion',     0.85, ...
                    'LineNoiseCriterion',   4, ...
                    'BurstCriterion',       'off', ...
                    'WindowCriterion',      'off', ...
                    'Highpass',             'off'
                }; %#ok

                evalc("EEG_cleaned = clean_artifacts(EEG, args{:});");

                Clean_Channel_Mask = EEG_cleaned.etc.clean_channel_mask; %#ok
                evalc("EEG = pop_select( EEG, 'channel', {EEG.chanlocs(find(Clean_Channel_Mask)).labels});");

            end

            %% Step 8: Artifact Rejection with ASR
            if PREPROC.Artifacts
                args = {}; %#ok
                args = {...
                    'FlatlineCriterion',    'off', ...
                    'ChannelCriterion',     'off', ...
                    'LineNoiseCriterion',   'off', ...
                    'BurstCriterion',       50, ...
                    'WindowCriterion',      'off', ...
                    'Highpass',             'off', ...
                    'BurstRejection',       'on', ...
                    'Distance',             'Euclidian'
                }; %#ok

                evalc("EEG_asr = pop_clean_rawdata(EEG, args{:});");
                Clean_Segment_Mask = EEG_asr.etc.clean_sample_mask';
                
                retain_data_intervals = reshape(find(diff([false Clean_Segment_Mask' false])),2,[])';
                retain_data_intervals(:,2) = retain_data_intervals(:,2)-1; %#ok
                evalc("EEG = pop_select(EEG, 'point', retain_data_intervals);");
                
            end

            %% Step 9: wICA
            if PREPROC.wICA
                % Run ICA
                addpath("wICA-master")
                [IC_activations, A, W] = fastica(double(EEG.data), 'verbose','off');
                EEG.icaweights = W;
                EEG.icasphere  = eye(size(EEG.icaweights,2));
                EEG.icawinv    = A;
                evalc("EEG = eeg_checkset(EEG, 'ica');");

                % ICLabel
                evalc("EEG = iclabel(EEG);");
                classifications = EEG.etc.ic_classification.ICLabel.classifications;
                artifact_ICs = find(classifications(:,3)>0.7 & classifications(:,1)<0.7);

                % Remove strong Artifacts with wICA
                Kthr = 1.25;
                ArtefThreshold = 4;
                IC_activations_clean = RemoveStrongArtifacts(IC_activations, artifact_ICs, Kthr, ArtefThreshold, EEG.srate, 'off');

                % Project to Original Data
                Data_clean = A * IC_activations_clean;
                EEG.data = Data_clean;
                evalc("EEG = eeg_checkset(EEG);");
            end

            %% Step 10: Interpolation
            evalc("EEG = pop_select( EEG, 'chantype','EEG');");
            n_interp = length(EEG_interp.chanlocs) - length(EEG.chanlocs);

            % Check number of n_interp and exclude
            if n_interp > 5
                fprintf('Too many interpolated Channels (%d) - %s: Condition: %s\n', n_interp, Raw_File.name, Cond_FileName);

                ExludeDataFile = fullfile(dir_Log, ['Error_TooManyExcludedChannels_', FileName, '.txt']);
                fid1 = fopen( ExludeDataFile, 'wt' );
                fprintf(fid1, 'Excluded Data-File= %s \n Too many excluded channels %i. \n', FileName,  n_interp);
                fclose(fid1);
                continue                

            end 

            evalc("EEG = pop_interp(EEG, EEG_interp.chanlocs, 'spherical');");
            evalc("EEG = eeg_checkset(EEG);");

            %% Step 11: Re-referencing
            if PREPROC.CAV_Reference
                evalc("EEG = pop_reref(EEG,[],'refloc',struct('labels',{'FCZ'},'type',{'EEG'},'ref', [], 'urchan', [], 'theta',{0},'radius',{0.127},'X',{0.388},'Y',{0},'Z',{0.922},'sph_theta',{0},'sph_phi',{67.2},'sph_radius',{1}, 'sph_theta_besa', [], 'sph_phi_besa', []));");
            end            

            %% Step 12: Surface Laplacian
            if PREPROC.Surf_Lap
                evalc("EEG = pop_currentdensity(EEG, 'method', 'finite');");
            end

            %% Step 13: Check Channel Number
            if length(EEG.chanlocs) ~= 59
                msg = ['Incorrect Channel Number: ' num2str(length(EEG.chanlocs))];
                error(msg)
            end

            % Order Channels and Data
            chanOrder = { ...
                'FP1','FP2', ...
                'AF7','AF3','AF4','AF8', ...
                'F7','F5','F3','F1','FZ','F2','F4','F6','F8', ...
                'FT7','FC5','FC3','FC1','FCZ','FC2','FC4','FC6','FT8', ...
                'T7','C5','C3','C1','CZ','C2','C4','C6','T8', ...
                'TP7','CP5','CP3','CP1','CPZ','CP2','CP4','CP6','TP8', ...
                'P7','P5','P3','P1','PZ','P2','P4','P6','P8', ...
                'PO7','PO3','POZ','PO4','PO8', ...
                'O1','OZ','O2' ...
            };
            [~, chIdx] = ismember(lower(chanOrder), lower({EEG.chanlocs.labels}));
            EEG.data = EEG.data(chIdx, :);
            EEG.chanlocs = EEG.chanlocs(chIdx);
            evalc("EEG = eeg_checkset(EEG);");

            %% Step 14: Epoching, Post-Epoching Artifacts, and Saving Data
            % Remove boundary events ?!?!?!
            EEG.event = EEG.event(~strcmp({EEG.event.type}, 'boundary'));
            switch PREPROC.Epoching
                case 'all'
                    % 6-seconds with 2-seconds overlap
                    evalc("EEG_oAEC = eeg_regepochs(EEG, 4, [0 6], 0, 'X', {}, 'on');");
                    n_epochs_oaec_before = size(EEG_oAEC.data, 3);
                    
                    % 12-seconds with 4-seconds overlap
                    evalc("EEG_phase = eeg_regepochs(EEG, 8, [0 12], 0, 'X', {}, 'on');");
                    n_epochs_phase_before = size(EEG_phase.data, 3);
                    
                    if PREPROC.Artifacts2
                        % Artifact Removal
                        evalc("EEG_oAEC = post_epoch_artifacts(EEG_oAEC);");
                        n_epochs_oaec_after = size(EEG_oAEC.data, 3);
                        removed_epochs_oaec = n_epochs_oaec_before - n_epochs_oaec_after;
                        evalc("EEG_oAEC = pop_saveset(EEG_oAEC, 'filename', FileName_oAEC, 'filepath', char(dir_Preproc), 'savemode', 'onefile');");
                        
                        evalc("EEG_phase = post_epoch_artifacts(EEG_phase);");
                        n_epochs_phase_after = size(EEG_phase.data, 3);
                        removed_epochs_phase = n_epochs_phase_before - n_epochs_phase_after;
                        evalc("EEG_phase = pop_saveset(EEG_phase, 'filename', FileName_phase, 'filepath', char(dir_Preproc), 'savemode', 'onefile');");
                    end

                case 'oAEC'
                    evalc("EEG_oAEC = eeg_regepochs(EEG, 4, [0 6], 0, 'X', {}, 'on');");
                    n_epochs_oaec_before = size(EEG_oAEC.data, 3);

                    if PREPROC.Artifacts2
                        % Artifact Removal
                        evalc("EEG_oAEC = post_epoch_artifacts(EEG_oAEC);");
                        n_epochs_oaec_after = size(EEG_oAEC.data, 3);
                        removed_epochs_oaec = n_epochs_oaec_before - n_epochs_oaec_after;
                        evalc("EEG_oAEC = pop_saveset(EEG_oAEC, 'filename', FileName_oAEC, 'filepath', char(dir_Preproc), 'savemode', 'onefile');");
                    end

                case 'phase'
                    evalc("EEG_phase = eeg_regepochs(EEG, 8, [0 12], 0, 'X', {}, 'on');");
                    n_epochs_phase_before = size(EEG_phase.data, 3);
                    
                    if PREPROC.Artifacts2
                        % Artifact Removal
                        evalc("EEG_phase = post_epoch_artifacts(EEG_phase);");
                        n_epochs_phase_after = size(EEG_phase.data, 3);
                        removed_epochs_phase = n_epochs_phase_before - n_epochs_phase_after;
                        evalc("EEG_phase = pop_saveset(EEG_phase, 'filename', FileName_phase, 'filepath', char(dir_Preproc), 'savemode', 'onefile');");
                    end

            end

            %% Step 15: Log for Quality Assessment
            if ~exist('EEG_oAEC','var')
                EEG_oAEC.trials = NaN;
            end

            if ~exist('EEG_phase','var')
                EEG_phase.trials = NaN;
            end

            if ~exist('removed_epochs_oaec','var')
                removed_epochs_oaec = NaN;
            end

            if ~exist('removed_epochs_phase','var')
                removed_epochs_phase = NaN;
            end

            if ~exist('artifact_ICs','var')
                artifact_ICs = NaN;
            end

            ID = strsplit(Raw_File.name, '_');
            ID = ID{1};

            cond_splits = strsplit(Cond_FileName, '_');
            run = cond_splits{1};
            condition = [cond_splits{3} '_' cond_splits{4}]; 

            Log_table = table({ID}, {run}, {condition}, EEG_oAEC.trials, removed_epochs_oaec, EEG_phase.trials, removed_epochs_phase, numel(artifact_ICs), n_interp, ...
                'VariableNames', {'ID', 'Run', 'Condition', 'oAEC_Epochs', 'Removed_oAEC_Epochs', 'Phase_Epochs', 'Removed_Phase_Epochs', 'ICs_removed', 'Interpolated_Channels'});
                        
            log_filename = fullfile(dir_Log, ['Log_', ID, '_', Cond_FileName, '.csv']);
            writetable(Log_table, log_filename);

        catch e
            % If error ocurrs, create ErrorMessage
            ErrorMessage = string(e.message);
            for ierrors = 1:length(e.stack)
                ErrorMessage = strcat(ErrorMessage, " // ", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
            end
            
            fprintf('***Error in File: %s;\n%s.\n', FileName, ErrorMessage);
            
            ErrorFile = fullfile(dir_Log, ['Error_PreProc', '_', strrep(FileName, '.set', '.txt')]);
            fid1 = fopen( ErrorFile, 'wt' );
            fprintf(fid1, 'Error-File: %s \nThe returned Error Message is: \n\n%s \n', FileName,  ErrorMessage);
            fclose(fid1);
            
        end % try-catch
        
    end % for condition

catch e % If error ocurrs, create ErrorMessage
    
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
    ErrorFile = fullfile(dir_Log, ['Error_PreProc_', ErrorFile, '_', Cond_FileName , '.txt']);
    fid1 = fopen(ErrorFile, 'wt');
    fprintf(fid1, 'Error-Subject: %s \nThe returned Error Message is: \n\n%s \n', InputFile,  ErrorMessage);
    fclose(fid1);

end % try-catch

end % function definition

function EEG_out = post_epoch_artifacts(EEG_in)
    Clean_Epochs_Mask = ones(EEG_in.trials, 1);
    threshold_DB = 90;
    threshold_SD = 3.29;
    % Demean Data otherwise Thresholding does not work
    EEG_in.data = EEG_in.data - mean(EEG_in.data,2);
    
    % Frequency Spectrum
    [~, bad_Spectrum] = pop_rejspec(EEG_in, 1, 'elecrange', 1:EEG_in.nbchan, 'threshold', [-threshold_DB threshold_DB], 'freqlimits', [1 30]);
    Clean_Epochs_Mask(bad_Spectrum) = 0;
    
    % Kurtosis
    bad_Kurtosis = pop_rejkurt(EEG_in, 1, 1:EEG_in.nbchan,  threshold_SD,threshold_SD,0,0,0);
    bad_Kurtosis = bad_Kurtosis.reject.rejkurt;
    Clean_Epochs_Mask(bad_Kurtosis) = 0;
    
    % Probability
    bad_Probability = pop_jointprob(EEG_in, 1, 1:EEG_in.nbchan,  threshold_SD, threshold_SD,0,0,0);
    bad_Probability = bad_Probability.reject.rejjp;
    Clean_Epochs_Mask(bad_Probability) = 0;

    if  sum(Clean_Epochs_Mask) == 0
        e.message = 'All Trials marked as bad (100%!!) .';
        error(e.message);
    end
    EEG_out = pop_select(EEG_in, 'trial',find(Clean_Epochs_Mask));
end
