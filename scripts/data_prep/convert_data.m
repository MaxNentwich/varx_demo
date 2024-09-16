%% Reorganize data and select data relevant for this analysis

% Source directories
data_dir = '/media/DATA/ECoGData/Tobii/Patients';
saccade_dir = '/media/DATA/ieeg_plot/Data/saccade_data';

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
        load(sprintf('%s/%s', data_pat, vids(v).name), 'envelope_ds', 'fs_ds', 'time_ds')
        envelope = envelope_ds;
        fs = fs_ds;
    
        save(sprintf('%s/%s', sample_pat, vids(v).name), 'envelope', 'fs')
    
        load(sprintf('%s/%s_%s', saccade_dir, patients{i}, vids(v).name), 'saccade_onset', 'fixation_onset', 'eye')
        eye_time = eye.time;
    
        fix_sample = round(interp1(time_ds, 1:length(time_ds), eye_time(fixation_onset == 1)));
        fixations = zeros(length(time_ds),1);
        fixations(fix_sample) = 1;
    
        sac_sample = round(interp1(time_ds, 1:length(time_ds), eye_time(saccade_onset == 1)));
        saccades = zeros(length(time_ds),1);
        saccades(sac_sample) = 1;
    
        save(sprintf('%s/%s', sample_pat_sim, vids(v).name), 'saccades', 'fixations', 'fs')
    
        % LFP
        if exist(sprintf('%s/%s', data_lfp_pat, vids(v).name), 'file') ~= 0

            load(sprintf('%s/%s', data_lfp_pat, vids(v).name), 'ieeg')
            lfp = ieeg.data;
            fs_lfp = ieeg.fs;

            % Resample 
            lfp = resample(lfp, 1e4, round(1e4*length(lfp)/length(envelope)));
    
            save(sprintf('%s/%s', sample_lfp_pat, vids(v).name), 'lfp', 'fs')

        end

    end

end