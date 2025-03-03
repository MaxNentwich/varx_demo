%% Load file correspondence sheet (TDT blocks and eyetracking files)

function [ieeg_names, eye_names] = file_corresp(data_dir, patient, task)

    % Load list to match ECoG and eyetracking files
    T = readtable(sprintf('%s/%s/%s_ExperimentNotes.xlsx', data_dir, patient, patient));

    eye_names = T.EyetrackingFile(cellfun(@(C) ~strcmp(C, ''), T.EyetrackingFile));

    ieeg_names = T.EEGFolder(cellfun(@(C) ~strcmp(C, ''), T.EyetrackingFile));

    % Select only files of the task being processed
    idx_task = cellfun(@(C) contains(C, task(1:end-1), 'IgnoreCase', true), eye_names);

    ieeg_names = ieeg_names(idx_task);
    eye_names = eye_names(idx_task);
    
end