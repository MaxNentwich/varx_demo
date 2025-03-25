%% Load iEEG and eyetracking data, align and save as .mat file

function load_data(options)

% Load the file name correspondance
T = readtable('file_names.xlsx');
        
for p = 1:length(options.patients)
    
    % Load file name correspondence (iEEG and eyetracking)
    [ieeg_names, eye_names] = file_corresp(options.data_dir, options.patients(p).name, options.task);
    
    % Define the filenames of the output directory 
    out_dir = sprintf('%s/%s/matlab_data', options.data_dir, options.patients(p).name);
    
    % Loop over task blocks
    for n = 1:length(ieeg_names)
        
        % Define output file names
        file_name_eye = sprintf('%s/Eyetracking/%s.mat', out_dir, T.mat_file{cellfun(@(C) contains(eye_names{n}, C), ...
            T.eye_file)});
        file_name_ieeg = sprintf('%s/Neural/%s.mat', out_dir, T.mat_file{cellfun(@(C) contains(eye_names{n}, C), ...
            T.eye_file)});

        if sum(ismember(eye_names, 'eyes_closed_rest')) == 2 && exist(file_name_ieeg, 'file') ~= 0
            file_name_ieeg = strrep(file_name_ieeg, 'Eyes_Closed_Rest.mat', 'Eyes_Closed_Rest_2.mat');
            n = 2;
        end

        % Skip if both output files exist
        if strcmp(options.task, 'Eyes_Closed_Rest')
            if exist(file_name_ieeg, 'file') ~= 0, continue, end
        else
            if exist(file_name_eye, 'file') ~= 0 && exist(file_name_ieeg, 'file') ~= 0
                continue
            end
        end

        fprintf('Loading subject %s, %s \n', options.patients(p).name, strrep(eye_names{n}, '.mat', ''))
        
        %% Load data
        ieeg = load_NS(options.data_dir, options.patients(p).name, options.task, ieeg_names{n}); 
        if ~strcmp(options.task, 'Eyes_Closed_Rest')
            eye = load_eye(options.data_dir, options.patients(p).name, options.task, eye_names{n});
        end

        %% Manual edits of triggers 
       
        % NS131 has an extra trigger in iEEG data
        if strcmp(options.patients(p).name, 'NS131') && strcmp(eye_names{n}, 'Eye_movie_sublij131_Monkey1avi_rep1.mat')
            
            ieeg.triggers = ieeg.triggers(2:end);
        
        % NS151 is missing triggers in iEEG data at the end. The triggers are taken from eyetracking data. This might 
        % introduce smallerrors in the range of ~3-4 ms 
        elseif strcmp(options.patients(p).name, 'NS151') && contains(eye_names{n}, 'The_Present')
            
           time_diff = mean(eye.triggers.time(1:length(ieeg.triggers))/1e6 - ieeg.triggers);
           
           ieeg.triggers(end+1 : length(eye.triggers.time)) = eye.triggers.time(length(ieeg.triggers)+1 : end)/1e6 - ...
               time_diff;
            
        end
            
        % Check if all triggers are in both TDT and eyetracking data
        if ~strcmp(options.task, 'Eyes_Closed_Rest')

            if length(ieeg.triggers) ~= length(eye.triggers.time)
                warning('Number of triggers in iEEG and Eyetracking don''t match! Have been corrected!')
                 
                [indx_ieeg, indx_eye] = match_triggers(ieeg.triggers, eye.triggers.time);
                
                ieeg.triggers = ieeg.triggers(indx_ieeg);
                eye.triggers.ID = eye.triggers.ID(indx_eye);
                eye.triggers.time = eye.triggers.time(indx_eye);
                eye.triggers.time_SDK = eye.triggers.time_SDK(indx_eye);
                
            end

        end

        % NS155_02: Monkey2 TDT block is named fixation -> rename
        % ieeg_names{n} or some data is processed like a fixation recording
        if strcmp(options.patients(p).name, 'NS155_02') && strcmp(ieeg_names{n}, 'B26_Fixation')
            ieeg_names{n} = 'Bedit_Monkey2';
        end
        
        %% Find trigger samples and check accuracy

        % Find Video Start trigger 
        if ~strcmp(options.task, 'Eyes_Closed_Rest')

            if contains(ieeg_names{n}, 'Fixation')
    
                idx_start = 1;
                idx_end = 2;
    
            else
    
                idx_start = find(eye.triggers.ID == options.trigger_IDs.start_ID);
                idx_end = find(eye.triggers.ID  == options.trigger_IDs.end_ID_1);
        
                if isempty(idx_end)
                    idx_end = find(eye.triggers.ID == options.trigger_IDs.end_ID_2);
                end
    
            end
    
            % Align time vectors to start triggers
            indx_eye = round(interp1(eye.time, 1:size(eye.time,2), eye.triggers.time));
    
            eye_T0 = eye.time(indx_eye(idx_start));
            eye.time = eye.time - eye_T0;
            eye.time = eye.time./1e6;
            
            % NS151 has missing triggers and data, fill in the missing data with NaNs
            if max(ieeg.triggers) > max(ieeg.time)
               
                ieeg.time = [ieeg.time, ieeg.time(end)+1/ieeg.fs:1/ieeg.fs:ieeg.triggers(end)];
                ieeg.data = [ieeg.data, nan(size(ieeg.data,1), length(ieeg.time)-length(ieeg.data))];
                
                warning(sprintf(['%.3fs of data missing at the end of the iEEG recording! Filled with NaNs based on', ...
                'eyetracking tiggers... \n'], sum(isnan(ieeg.data(1,:)))/ieeg.fs))
                
            end
    
            indx_ieeg = round(interp1(ieeg.time, 1:size(ieeg.time,2), ieeg.triggers));
    
            ieeg_T0 = ieeg.time(indx_ieeg(idx_start));
            ieeg.time = ieeg.time - ieeg_T0;
    
            % Check if distance between triggers in TDT and eyetracking data 
            % is similar
            jit = max(abs(1e3*(ieeg.time(indx_ieeg(idx_start:idx_end)) - ...
                eye.time(indx_eye(idx_start:idx_end)))));
    
            if jit > 1e3*(1/eye.fs)
                warning('Jitter between EEG and Eyetracking triggers: %4.1f ms', jit)
            end
    
            % Find movie start and end trigger
            start_eye = eye.triggers.time(idx_start);
            end_eye = eye.triggers.time(idx_end);
    
            start_ieeg = ieeg.triggers(idx_start);
            end_ieeg = ieeg.triggers(idx_end);
    
            % Check if the time between start and end trigger of ECoG and 
            % eyetracking are the same
            se_diff = ((end_eye - start_eye)/1e6) - (end_ieeg - start_ieeg);
            if se_diff > (1/eye.fs)
                warning('Difference of movie start and end trigger between ECoG and Eyetracking: %4.1f ms', 1e3*se_diff)
            end
            
            %% Organize trigger data
            trigger_ID = eye.triggers.ID;
            
            % iEEG
            trigger_time_original = ieeg.triggers; 
            trigger_time = ieeg.triggers - ieeg_T0;
            trigger_sample = indx_ieeg;
            ieeg.triggers = table(trigger_ID, trigger_time, trigger_sample, trigger_time_original);
    
            % Eyetracking
            trigger_time_original = eye.triggers.time;
            trigger_time = (eye.triggers.time - eye_T0)./1e6;
            trigger_sample = indx_eye;
            trigger_time_SDK = eye.triggers.time_SDK;
    
            eye.triggers = table(trigger_ID, trigger_time, trigger_sample, trigger_time_original, ...
                trigger_time_SDK);
            
            %% Timing of frames 
    
            if ~contains(ieeg_names{n}, 'Fixation')
    
                % SDK time vector interpolalated from tigger timing (has a linear drift to eye.time)
                eye.time_SDK = interp1(eye.triggers.trigger_time, eye.triggers.trigger_time_SDK, eye.time);
        
                % Find samples of eyetracking data that correspond to each frame in the video
                sample_SDK = 1:length(eye.time_SDK);             
                eye.frame_sample = round(interp1(eye.time_SDK(~isnan(eye.time_SDK)), sample_SDK(~isnan(eye.time_SDK)), ...
                    eye.frame_time_SDK));
        
                eye.frame_time = eye.time(eye.frame_sample);
            
            end

        end

        %% Save the data 
        if ~strcmp(options.task, 'Eyes_Closed_Rest')

            if exist(sprintf('%s/Eyetracking', out_dir),'dir') == 0
                mkdir(sprintf('%s/Eyetracking', out_dir))
            end

            save(file_name_eye, 'eye')

        end
        
        if exist(sprintf('%s/Neural', out_dir), 'dir') == 0
            mkdir(sprintf('%s/Neural', out_dir))
        end
           
        save(file_name_ieeg, 'ieeg', '-v7.3')

    end
    
end

end