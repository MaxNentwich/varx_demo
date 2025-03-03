%% Load eyetracking data

function eye = load_eye(data_dir, patient, task, eye_name)

    if strcmp(task, 'Fixation')

        load(sprintf('%s/%s/Eyetracking/%s/%s', data_dir, patient, task, eye_name), 'Eye_fixation')

        eye.left = Eye_fixation.LEC_fv';
        eye.right = Eye_fixation.REC_fv'; 

        % Determine Eyetracking sampling rate
        eye.fs = Eye_fixation.setup.fixDur;

        % Process triggers
        eye.triggers.time = double(cell2mat(Eye_fixation.timing.ET_time(:,2)));
        eye.triggers.time_SDK = double(cell2mat(Eye_fixation.timing.SDK_time(:,2)))./1e6;
        
        eye.triggers.ID = Eye_fixation.timing.ET_time(:,1);
        eye.triggers.ID = cellfun(@(C) str2double(C(2:end)), eye.triggers.ID);

        % Time vector
        eye.time = double(Eye_fixation.TS_fv)';

        % Additional parameters
        eye.setup = Eye_fixation.setup;

    else

        load(sprintf('%s/%s/Eyetracking/%s/%s', data_dir, patient, task, eye_name), 'Eye_movie')

        eye.left = Eye_movie.leftEye';
        eye.right = Eye_movie.rightEye'; 

        % Determine Eyetracking sampling rate
        eye.fs = Eye_movie.setup.currentFrameRate;

        % Process triggers
        eye.triggers.time = double(cell2mat(Eye_movie.timing.ET_time(:,2)));
        eye.triggers.time_SDK = double(cell2mat(Eye_movie.timing.SDK_time(:,2)))./1e6;
        
        eye.triggers.ID = Eye_movie.timing.ET_time(:,1);
        eye.triggers.ID = cellfun(@(C) str2double(C(2:end)), eye.triggers.ID);

        % Time vector
        eye.time = double(Eye_movie.timeStamp)';

        % Additional parameters
        eye.n_fr = Eye_movie.fr;
        eye.frame_time_SDK = Eye_movie.vbl_aggr;
        eye.setup = Eye_movie.setup;

    end

    % Check if there are gaps in the eyetraking timestamps
    % This happens because Tobii can only hold 600 s in buffer. For longer
    % recordings some data is not saved 
    
    if max(diff(eye.time))/median(diff(eye.time)) > 2

        [~, idx_shift] = max(diff(eye.time));
        idx_shift = idx_shift + 1;

        % Insert NaN samples
        time_real = mean(diff(eye.time(idx_shift:end)));

        time_eye_gap = eye.time(idx_shift-1)+time_real : time_real : eye.time(idx_shift)-time_real;

        eye.time = [eye.time(1:idx_shift-1), time_eye_gap, eye.time(idx_shift:end)];

        eye.left  = [eye.left(:, 1:idx_shift-1), ...
            NaN([size(eye.left ,1), size(time_eye_gap,2)]), ...
            eye.left(:, idx_shift:end)];

        eye.right  = [eye.right(:, 1:idx_shift-1), ...
            NaN([size(eye.right ,1), size(time_eye_gap,2)]), ...
            eye.right(:, idx_shift:end)];

        warning('%4.1f ms worth of NaNs filled into gap in eyetracking data!', size(time_eye_gap,2)*time_real/1e3);

    end
    
end