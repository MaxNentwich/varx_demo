%% Extract different frequency bands and phase information by wavelet tranform 

function extract_bandpower_phase(options)

for pat = 1:length(options.patients)

    data_files = dir(sprintf('%s/%s/%s', options.data_dir, options.patients(pat).name, options.neural_prep_dir));
    data_files([data_files.isdir]) = []; 
    
    data_files(cellfun(@(C) contains(C, 'artifacts'), extractfield(data_files, 'name'))) = [];

    env_out_dir = sprintf('%s/%s/%s/%s', options.data_dir, options.patients(pat).name, options.env_dir);

    for v = 1:length(data_files)
        
        % Check if the output for any band is missing
        for f = 1:length(options.freq_bands)
            file_check(f) = exist(sprintf('%s/%s/%s', env_out_dir, options.band_names{f}, data_files(v).name), 'file') ~= 0;
        end
        
        if sum(file_check == 0) == 0
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
        end

        data = ieeg.data;
        ieeg.data = [];
        
        fs_ieeg = ieeg.fs;
        time_ieeg = ieeg.time;
        
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

        % Bandpass filter 
        for f = 1:length(options.freq_bands)
            
            % Define output directory and file name
            out_dir = sprintf('%s/%s', env_out_dir, options.band_names{f});
            if exist(out_dir, 'dir') == 0, mkdir(out_dir), end
            
            out_file = sprintf('%s/%s', out_dir, data_files(v).name);            
            if exist(out_file, 'file') ~= 0, continue, end
            
            fprintf('%s ... \n', options.band_names{f})
            
            ieeg.band = [];
            
            % Hilbert transform of bandpassed signal
            for ch = 1:size(data,2)
                ieeg.band(:,ch) = BpFilter(data(:,ch), 5, options.freq_bands{f}, fs_ieeg);
                ieeg.hilbert(:,ch) = hilbert(ieeg.band(:,ch)); 
            end
            
            % Power envelope 
            ieeg.envelope = abs(ieeg.hilbert);   
            
            % Instantaneous phase 
            ieeg.phase = angle(ieeg.hilbert); 
            
            ieeg.hilbert = [];
            
            % Lowpass filter envelope and bandsfor anti-aliasing 
            for ch = 1:size(ieeg.envelope,2)
                ieeg.envelope(:,ch) = LpFilter(ieeg.envelope(:,ch), 5, eye.fs/2, fs_ieeg);
                ieeg.band(:,ch) = LpFilter(ieeg.band(:,ch), 5, eye.fs/2, fs_ieeg);
            end
            
            % Downsample data to eyetracking data
            if contains(data_files(v).name, 'Eyes_Closed_Rest')
               eye.time = ieeg.time(1):1/eye.fs:ieeg.time(end);
            end

            ieeg.envelope = interp1(time_ieeg, ieeg.envelope, eye.time);
            ieeg.band = interp1(time_ieeg, ieeg.band, eye.time);
            ieeg.phase = interp1(time_ieeg, ieeg.phase, eye.time);

            % Cut data at triggers
            if contains(data_files(v).name, 'Eyes_Closed_Rest')
                
                ieeg.time = eye.time;
                ieeg.fs = eye.fs;

            else

                ieeg.envelope = ieeg.envelope(start_sample:end_sample, :);
                ieeg.band = ieeg.band(start_sample:end_sample, :);
                ieeg.phase = ieeg.phase(start_sample:end_sample, :);
                
                % Time and sampling rate
                ieeg.time = eye.time(start_sample:end_sample);
                ieeg.fs = eye.fs;

            end


            % Resample the envelopes to 60 Hz 
            for ch = 1:size(ieeg.envelope,2)
                ieeg.envelope(:,ch) = LpFilter(ieeg.envelope(:,ch), 5, ieeg.fs/options.dsf/2, ieeg.fs);
                ieeg.band(:,ch) = LpFilter(ieeg.band(:,ch), 5, ieeg.fs/options.dsf/2, ieeg.fs);
            end
            time_ds = linspace(ieeg.time(1), ieeg.time(end), round(length(ieeg.envelope)/options.dsf));
            envelope_ds = interp1(ieeg.time, ieeg.envelope, time_ds);
            fs_ds = eye.fs/options.dsf;    
            
            % Resample to 10 Hz
            time_10Hz = linspace(ieeg.time(1), ieeg.time(end), round(length(ieeg.envelope)/(ieeg.fs/10)));
            envelope_10Hz = interp1(ieeg.time, ieeg.envelope, time_10Hz);

            % Resample the band
            band_ds = interp1(ieeg.time, ieeg.band, time_ds);
            band_10Hz = interp1(ieeg.time, ieeg.envelope, time_10Hz);
            
            % Resample phase
            phase_ds = interp1(ieeg.time, ieeg.phase, time_ds);

            % Save data
            labels = ieeg.label;
            
            save(out_file, 'ieeg', 'envelope_ds', 'envelope_10Hz', 'band_ds', 'band_10Hz', 'time_ds', 'time_10Hz', ...
                'phase_ds', 'fs_ds', 'labels', '-v7.3')
 
        end
        
        clearvars start_sample end_sample data
                
    end
    
end

end