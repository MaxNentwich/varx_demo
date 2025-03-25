%% Reorganize data and select data relevant for this analysis

% Source directories
data_dir = '/media/max/e61df479-8f57-4855-b852-04ccdfb12a6c/ECoGData/Tobii/Patients';
saccade_dir = '/media/max/Workspace/Data/saccade_data';

% Output directory
sample_dir = '/media/max/Workspace/Data/movies_legacy_sample';

% List of patients 
patient_list = readtable('../../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

for i = 1:length(patients)

    data_pat = sprintf('%s/%s/matlab_data/Envelope_phase/BHA', data_dir, patients{i});
    sample_pat = sprintf('%s/%s/BHA', sample_dir, patients{i});
    sample_pat_sim = sprintf('%s/%s/Stimuli', sample_dir, patients{i});
    
    vids = dir(data_pat);
    vids([vids.isdir]) = [];

    % LFP
    data_lfp_pat = sprintf('%s/%s/matlab_data/raw', data_dir, patients{i});
    sample_lfp_pat = sprintf('%s/%s/LFP', sample_dir, patients{i});

    if exist(sample_pat, 'dir') == 0, mkdir(sample_pat), end
    if exist(sample_pat_sim, 'dir') == 0, mkdir(sample_pat_sim), end
    if exist(sample_lfp_pat, 'dir') == 0, mkdir(sample_lfp_pat), end

    for v = 1:length(vids)
    
        % HFA
        hfa_file = sprintf('%s/%s', sample_pat, vids(v).name);

        load(sprintf('%s/%s', data_pat, vids(v).name), 'envelope_ds', 'fs_ds', 'time_ds')
        envelope = envelope_ds;
        fs = fs_ds;
    
        if exist(hfa_file, 'file') == 0
            save(hfa_file, 'envelope', 'fs')
        end
    
        % Saccades and fixation
        if strcmp(vids(v).name, 'Eyes_Closed_Rest.mat') || strcmp(vids(v).name, 'Eyes_Closed_Rest_2.mat')

            saccades = [];
            fixations = [];
            pupil = [];

        else

            load(sprintf('%s/%s_%s', saccade_dir, patients{i}, vids(v).name), 'saccade_onset', 'fixation_onset', 'pupil', 'eye')
            eye_time = eye.time;
        
            fix_sample = round(interp1(time_ds, 1:length(time_ds), eye_time(fixation_onset == 1)));
            fixations = zeros(length(time_ds),1);
            fixations(fix_sample) = 1;
        
            sac_sample = round(interp1(time_ds, 1:length(time_ds), eye_time(saccade_onset == 1)));
            saccades = zeros(length(time_ds),1);
            saccades(sac_sample) = 1;
    
            % Resample Pupil
            pupil = interp1(eye_time, pupil, time_ds);

        end
    
        % Frame corresponding to saccade
        if ~strcmp(vids(v).name, 'Resting_fixation.mat') && ~strcmp(vids(v).name, 'Eyes_Closed_Rest.mat')

            idx_frame_time = 1:length(eye.frame_time);

            idx_duplicate = find(diff(eye.frame_time) == 0);

            eye.frame_time(idx_duplicate) = [];
            idx_frame_time(idx_duplicate) = [];

            frame_saccade = round((interp1(eye.frame_time, idx_frame_time, eye_time(saccade_onset == 1))));
            frame_fixation = round((interp1(eye.frame_time, idx_frame_time, eye_time(fixation_onset == 1))));

        else
            frame_saccade = []; frame_fixation = [];
        end

        save(sprintf('%s/%s', sample_pat_sim, vids(v).name), 'saccades', 'fixations', 'frame_saccade', 'frame_fixation', 'pupil', 'fs')
    
        % LFP
        lfp_file = sprintf('%s/%s', sample_lfp_pat, vids(v).name);

        if exist(lfp_file, 'file') == 0

            if exist(sprintf('%s/%s', data_lfp_pat, vids(v).name), 'file') ~= 0
    
                load(sprintf('%s/%s', data_lfp_pat, vids(v).name), 'ieeg')
                lfp = ieeg.data;
                fs_lfp = ieeg.fs;
    
                % Resample 
                lfp = resample(lfp, 1e4, round(1e4*length(lfp)/length(envelope)));
        
                save(lfp_file, 'lfp', 'fs')
    
            end

        end

    end

end