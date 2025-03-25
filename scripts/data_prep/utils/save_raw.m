%% Save a downsampeled version of the iEEG data

function save_raw(options)

for pat = 1:length(options.patients)

    data_files = dir(sprintf('%s/%s/%s', options.data_dir, options.patients(pat).name, options.neural_prep_dir));
    data_files([data_files.isdir]) = [];
    
    data_files(cellfun(@(C) contains(C, 'artifacts'), extractfield(data_files, 'name'))) = [];
    
    out_dir = sprintf('%s/%s/%s', options.data_dir, options.patients(pat).name, options.raw_dir);

    for v = 1:length(data_files)
        
        out_file = sprintf('%s/%s', out_dir, data_files(v).name);
        
        if exist(out_file, 'file') ~= 0
            continue
        end
        
        fprintf('Processing %s, %s ...\n', options.patients(pat).name, data_files(v).name)
        
        % Load ieeg data
        load(sprintf('%s/%s', data_files(v).folder, data_files(v).name), 'ieeg')
        
        % Load eyetracking data
        if ~contains(data_files(v).name, 'Eyes_Closed_Rest')
            load(sprintf('%s/%s/%s/%s', options.data_dir, options.patients(pat).name, options.eye_dir, data_files(v).name), 'eye')
        else
            load(sprintf('%s/%s/%s/Despicable_Me_English.mat', options.data_dir, options.patients(pat).name, options.eye_dir), 'eye')
            eye.time = ieeg.time(1):1/eye.fs:ieeg.time(end);
        end
        
        % Trigger samples 
        if contains(data_files(v).name, 'Resting_fixation.mat')
            start_sample = eye.triggers.trigger_sample(1);
            end_sample = eye.triggers.trigger_sample(2);
        elseif ~contains(data_files(v).name, 'Eyes_Closed_Rest')
            start_sample = eye.triggers.trigger_sample(eye.triggers.trigger_ID == options.trigger_IDs.start_ID);
    
            end_sample = eye.triggers.trigger_sample(eye.triggers.trigger_ID == options.trigger_IDs.end_ID_1);
            if isempty(end_sample) 
                end_sample = eye.triggers.trigger_sample(eye.triggers.trigger_ID == options.trigger_IDs.end_ID_2);
            end
        end
        
        % Anti-aliasing filter
        for ch = 1:size(ieeg.data,2)
            ieeg.data(:,ch) = LpFilter(ieeg.data(:,ch), 5, eye.fs/2, ieeg.fs); 
        end
        
        % Downsample data to eyetracking data
        ieeg.data = interp1(ieeg.time, ieeg.data, eye.time);

        % Cut data at triggers
        if ~contains(data_files(v).name, 'Eyes_Closed_Rest')
            ieeg.data = ieeg.data(start_sample:end_sample, :);
            ieeg.time = eye.time(start_sample:end_sample);
        end

        ieeg.fs = eye.fs;
        
        % Save data
        if exist(out_dir, 'dir') == 0, mkdir(out_dir), end

        % Save
        save(out_file, 'ieeg', '-v7.3')

        clearvars start_sample end_sample

    end
    
end

end