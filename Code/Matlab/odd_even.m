function odd_even(dir_Preproc, dir_Log, PREPROC)
%% Odd-Even Split of Preprocessed Data
%
% Output:
%   (1) a .set file of the odd-even splitted preprocessed data
%
% Inputs:
%   dir_Preproc:        String pointing to the folder where preproc Data is
%                       and will be saved in extra folder
%   dir_Log:            String pointing to where Log-Files should be saved
%   PREPROC:            Struct with analysis settings
%
% Created by: Christoph Frühlinger
% Last edited: May 2026

%% Prepare List of Files to be Processed
Preproc_Files   = dir(fullfile(dir_Preproc, '**/*.set'));  

%% Directory where file should be saved
dir_Log         = fullfile(dir_Log, 'Connectivity');
dir_odd_even    = fullfile(dir_Preproc, 'OddEven');

if ~isfolder(dir_Log)
    mkdir(dir_Log)
end

if ~isfolder(dir_odd_even)
    mkdir(dir_odd_even)
end

fprintf(['\n%s\n' ...
         'Odd Even Splitting' ...
         '\n%s\n'], ...
         repmat('=', 1, 100), ...
         repmat('=', 1, 100));
fprintf('\nAnalyzing %d Files. \n', length(Preproc_Files));
fprintf(['InputFolder:  %s \n' ...
         'OutputFolder: %s \n' ...
         'LogFolder:    %s \n\n'], ...
         dir_Preproc , dir_odd_even, dir_Log);

%% Increase calculation speed by running multiple subjects in parallel
delete(gcp('nocreate')); % make sure that previous pooling is closed
if isempty(PREPROC.nWorkers)
    parpool("Processes")
else
    parpool("Processes", PREPROC.nWorkers);
end

nFiles          = length(Preproc_Files);


%% Loop over Files
q = parallel.pool.DataQueue;
afterEach(q, @(msg) fprintf('%s', msg));

parfor i_File = 1:nFiles

    try

        % Load Preprocessed Data
        InputFile = fullfile(Preproc_Files(i_File).folder, Preproc_Files(i_File).name);
        EEG = pop_loadset(InputFile);
        % Even Odd Split
        EEG_odd = pop_select(EEG, 'trial', 1:2:EEG.trials);
        EEG_even = pop_select(EEG, 'trial', 2:2:EEG.trials);
        
        % Save data
        [~, name, ~] = fileparts(Preproc_Files(i_File).name);
        OutputFileName_odd = [name, '_odd_connectivity.mat'];
        OutputFileName_even = [name, '_even_connectivity.mat'];
        
        pop_saveset(EEG_odd, 'filename', OutputFileName_odd, 'filepath', char(dir_odd_even), 'savemode', 'onefile');
        pop_saveset(EEG_even, 'filename', OutputFileName_even, 'filepath', char(dir_odd_even), 'savemode', 'onefile');
        send(q, sprintf('[%d/%d] %s finished\n', i_File, nFiles, name));

    catch e
        ErrorMessage = string(e.message);
        for ierrors = 1:length(e.stack)
            ErrorMessage = strcat(ErrorMessage, " // ", e.stack(ierrors).name, ", Line: ",  num2str(e.stack(ierrors).line));
        end
    
        fprintf('***Error in File: %s;\n%s.\n', InputFile, ErrorMessage);
    
    end % try-catch

end % i_File
