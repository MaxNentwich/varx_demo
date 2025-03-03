
function [eye, start_sample, end_sample] = load_et(data_dir, eye_dir, patient_name, file_name, trigger_IDs)

    load(sprintf('%s/%s/%s/%s', data_dir, patient_name, eye_dir, file_name), 'eye')

    % Trigger samples
    if contains(file_name, 'fixation')

        start_sample = eye.triggers.trigger_sample(1);
        end_sample = eye.triggers.trigger_sample(2);

    else

        start_sample = eye.triggers.trigger_sample(eye.triggers.trigger_ID == trigger_IDs.start_ID);
    
        end_sample = eye.triggers.trigger_sample(eye.triggers.trigger_ID == trigger_IDs.end_ID_1);
        if isempty(end_sample) 
            end_sample = eye.triggers.trigger_sample(eye.triggers.trigger_ID == trigger_IDs.end_ID_2);
        end

    end

    eye.left = eye.left(:, start_sample:end_sample);
    eye.right = eye.right(:, start_sample:end_sample);
    eye.time = eye.time(start_sample:end_sample);
    
end