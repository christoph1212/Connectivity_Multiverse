%% Main Analysis Script
%
% This script is the main file for the study "Robustness of EEG Functional 
% Brain Networks Associated with Fluid Intelligence: A Multiverse Analysis 
% of Connectivity and Thresholding Methods". 
% 
% It runs:
% 1. Preprocessing
% 2. Connectivity Analysis (Multiverse)
% 3. Thresholding (Multiverse)
% 4. Graph Theory Computation
%
% Make sure to set the folders according to your workspace and define
% analysis configurations.
%
% Created by: Christoph Frühlinger
% Last edited: May 2026

%% Housekeeping
clear
clc
close

%% Setup
current_file        = mfilename('fullpath');
[current_dir, ~, ~] = fileparts(current_file);
dir_Root            = fileparts(fileparts(current_dir));          % path to project
dir_Raw             = fullfile(dir_Root, 'Data', 'RawData');      % path to raw data
dir_Log             = fullfile(dir_Root, 'Data', 'Log');          % path where log data should be stored
dir_Preproc         = fullfile(dir_Root, 'Data', 'Preprocessed'); % path where preprocessed data should be stored (will be created)
dir_Connect         = fullfile(dir_Root, 'Data', 'Connectivity'); % path where connectivity data should be stored (will be created)
combine_conn_files  = true;                                       % combine connectivity files into one file? Single files will be deleted.
Overwrite           = true;                                       % overwrite existing files?
nWorkers            = 12;                                         % for parfor - use [] for max

% Start EEGLAB
dir_eeglab          = fullfile(dir_Root, 'Code', 'Matlab', 'eeglab2026.0.0'); % adjust accordingly
addpath(dir_eeglab);
eeglab nogui
clc

% Check for necessary Plugins
plugins             = ["scd", "Fieldtrip-lite"];

for i = 1:numel(plugins)
    if ~ismember(plugins(i), string({PLUGINLIST.plugin}))
        error("Error: Plugin %s is not installed. Please install.", plugins(i));
    end
end

fprintf(['%s\n' ...
         'Your folder settings:\n\n' ...
         'Root-Folder:                %s\n' ...
         'Raw-Folder:                 %s\n' ...
         'Log-Folder:                 %s\n' ...
         'Overwrite:                  %d\n' ...
         'Combine Connectivity Files: %d\n' ...
         '%s\n'], ...
         repmat('=', 1, 100), ...
         dir_Root, dir_Raw, dir_Log, Overwrite, combine_conn_files, ...
         repmat('=', 1, 100));

% Check if Brain Connectivity Toolbox exists
if exist(fullfile(userpath, '2019_03_03_BCT'), "dir") == 7
    addpath(fullfile(userpath, '2019_03_03_BCT'))
else
    error("Brain Connectivity Toolbox not found.\nPlease download and place in: %s\n\n", userpath)
end

%% Analysis Configurations
PREPROC = struct(...
    'nWorkers',         nWorkers, ... % Number of Workers for parfor
    'Downsample',       false, ...    % Downsample to SR/2
    'HP_Filter',        1, ...        % High-Pass Filter
    'LP_Filter',        30, ...       % Low-Pass Filter
    'BadChans',         true, ...     % Bad Channel Detection
    'Artifacts',        true, ...     % Artifact Rejection
    'wICA',             true, ...     % Wavelet-Enhanced ICA
    'Interpolate',      true, ...     % Interpolate Bad Channels    
    'CAV_Reference',    true, ...     % Re-reference to Common Average
    'Surf_Lap',         true, ...     % Apply Surface Laplacian
    'Epoching',         'all', ...    % Epoch Data for 'oAEC', 'phase' Based Measures or 'all' 
    'Artifacts2',       true ...      % Post Epoching Artifact Rejection
);

CONNECTIVITY = struct(...
    'nWorkers',         nWorkers, ... % Number of Workers for parfor
    'Bands',            'all', ...    % Frequency Bands: 'delta', 'theta', 'alpha1', 'alpha2', 'beta', or 'all'
    'Measures',         'all'  ...    % Connectivity Measures: 'imcoh', 'wpli', 'pli', 'pcoh', 'oaec', or 'all'
);

THRESHOLD = struct(...
    'nWorkers',         nWorkers, ... % Number of Workers for parfor
    'Method',           'all' ...     % Thresholding Method: 'dens', 'omst', 'eco', 'mcc', or 'all'
);

GRAPH = struct(...
    'nWorkers',         nWorkers, ... % Number of Workers for parfor
    'metrics',          'all' ...     % Graph Theory Metrics: 'cc', 'pathl', 'eglob', 'eloc', 'smallworld', or 'all'
);

fprintf([repmat('=', 1, 100), '\nYour analysis settings:\n\n   <strong>Preprocessing</strong>\n'])
disp(PREPROC)
fprintf('   <strong>Connectivity</strong>\n')
disp(CONNECTIVITY)
fprintf('   <strong>Thresholding</strong>\n')
disp(THRESHOLD)
fprintf('   <strong>Graph Theory</strong>\n')
disp(GRAPH)
fprintf([repmat('=', 1, 100), '\n']);

%% Preprocessing
preprocess_data(dir_Raw, dir_Log, dir_Preproc, PREPROC, Overwrite)
% TEST_preprocess_data(dir_Raw, dir_Log, dir_Preproc, PREPROC, Overwrite)

%% Connectivity Multiverse
connectivity_multiv(dir_Preproc, dir_Log, dir_Connect, CONNECTIVITY, combine_conn_files, Overwrite)

%% Thresholding Multiverse
thresholding_multiv(dir_Root, dir_Connect, dir_Log, THRESHOLD, Overwrite)

%% Graph Theory
graph_metrics(dir_Connect, dir_Log, GRAPH, Overwrite)