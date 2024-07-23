%% Select a patient -> add this variable everywhere
pat_select = 'NS135';

%% Types of recordings to preprocess
% inscapes, rest
data_type = 'inscapes';

%% Directories 
options.w_dir = pwd; 
options.drive_dir = '/media/DATA/ieeg_plot';                                        % Path to organized metadata
options.data_drive_dir = '/media/DATA/ECoGData';                                    % Path to raw data on hard drive
options.sample_dir = '/media/max/Workspace/Data/movies_legacy_sample';

options.data_dir = sprintf('%s/Tobii/Patients', options.data_drive_dir);            % Data directory
options.neural_prep_dir = 'matlab_data/Neural_prep';                                % Preprocessed neural data 
options.env_dir = 'matlab_data/Envelope_phase';                                     % Envelope and phase
options.eye_dir = 'matlab_data/Eyetracking';                                        % Eytracking data 
options.raw_dir = 'matlab_data/raw';                                                % LFP
options.im_data_dir = sprintf('%s/Data', options.drive_dir);                        % Data from intermediary analysis steps

% Add some necessary paths
addpath(genpath(sprintf('%s/Functions', options.w_dir)))                            % Functions
addpath(genpath(sprintf('%s/Organize', options.drive_dir)))                         % Tables and variables for organization
addpath(genpath(sprintf('%s/External', options.w_dir)))                             % External packages

%% Task (Movies, Freeviewing, Fixation, RhythmicSaccade) 
if strcmp(data_type, 'rest')
    options.task = 'Fixation';
elseif strcmp(data_type, 'inscapes')
    options.task = 'Movies';
end

%% Set flags for visualization and videos with eye movements or annotations
options.visualize_preproc = true;                                                   % Plot spectrograms and sample channels 
options.visualize_trfs = true;                                                      % Plot stimuli for trf analysis 

%% List of patients
options.patients = dir(options.sample_dir);
options.patients = options.patients(cellfun(@(C) strcmp(C, pat_select), {options.patients.name}));

%% Definition of frequency bands
options.freq_bands = [70,150];                                                      % Band limits in Hz
options.band_names = {'BHA'};                                                       % Names
options.dsf = 5;                                                                    % Downsampling factor for all frequency bands
options.band_select = {'BHA'};                                                      % Frequency band {'raw', 'Theta', 'Alpha', 'Beta', 'BHA'}

if strcmp(data_type, 'rest')
    options.vid_names = {'Resting_fixation'};                                       % Videos {'Monkey', 'Despicable_Me_English', 'Despicable_Me_Hungarian', 'The_Present_Rep_1', 'The_Present_Rep_2'}
elseif strcmp(data_type, 'inscapes')
    options.vid_names = {'Inscapes'};  
end

%% Saccade detection
options.vel_th = 2;                                                                 % Threshold for saccade detection (standard deviations)
options.peri_saccade_window = [0.033, 0.12];                                        % Window around saccades used to detect amplitude and speed, in seconds
options.min_diff = 0.11;                                                            % Minimum distance required between saccades (in seconds)
options.t_after = 1;                                                                % Time after cuts to look for saccades (in seconds)
options.t_before = 0.5;                                                             % Time before cuts to look for saccades (in seconds)
options.t_around = 1;                                                               % Time window around cuts to exclude from 'other' saccades
options.saccade_selection = 'first_last';                                           % include 'all' or the first and last ('first_last') saccade before and after cuts

options.screen_size = [1920, 1080];                                                 % Screen size of the Tobii eyetracker [px]
options.screen_dimension = [509.2 286.4];                                           % Screen size of the Tobii eyetracker [mm]

%% Triggers in eyetracking data
options.trigger_IDs.start_ID = 11;
options.trigger_IDs.end_ID_1 = 12;
options.trigger_IDs.end_ID_2 = 13;

%% Load the data
load_data(options)

%% Check if the labels make sense
label_check(options)

%% Preprocess the iEEG data
preprocess_ieeg(options)

%% Extract the data in different frequency bands
extract_bandpower_phase(options)

%% Save the downsampled 'raw' LFP 
save_raw(options)

%% Save saccade data
save_saccade_data(options)

%% Save only necessary data

% BHA
if strcmp(data_type, 'rest')

    load(sprintf('/media/DATA/ECoGData/Tobii/Patients/%s/matlab_data/Envelope_phase/BHA/Resting_fixation.mat', pat_select), 'envelope_ds', 'fs_ds', 'labels')
    load(sprintf('/media/DATA/ieeg_plot/Data/saccade_data/%s_Resting_fixation.mat', pat_select), 'saccade_onset', 'fixation_onset', 'eye')
    
    saccades = resample_peaks(saccade_onset, eye.fs, eye.fs/fs_ds);
    fixations = resample_peaks(fixation_onset, eye.fs, eye.fs/fs_ds);
    
    envelope = envelope_ds;
    fs = fs_ds;
    
    save(sprintf('%s/%s/BHA/Resting_fixation.mat', options.sample_dir, pat_select), 'envelope', 'fs')
    save(sprintf('%s/%s/Stimuli/Resting_fixation.mat', options.sample_dir, pat_select), 'saccades', 'fixations', 'fs')
    save(sprintf('%s/%s/ieeg_channel_labels.mat', options.sample_dir, pat_select), 'labels')
end

% Also add inscapes
if strcmp(data_type, 'inscapes')

    load(sprintf('%s/%s/matlab_data/Envelope_phase/BHA/Inscapes.mat', options.data_dir, pat_select), 'envelope_ds', 'fs_ds')
    load(sprintf('/media/DATA/ieeg_plot/Data/saccade_data/%s_Inscapes.mat', pat_select), 'saccade_onset', 'fixation_onset', 'eye')

    saccades = resample_peaks(saccade_onset, eye.fs, eye.fs/fs_ds);
    fixations = resample_peaks(fixation_onset, eye.fs, eye.fs/fs_ds);

    envelope = envelope_ds;
    fs = fs_ds;
    save(sprintf('%s/%s/BHA/Inscapes.mat', options.sample_dir, pat_select), 'envelope', 'fs')
    save(sprintf('%s/%s/Stimuli/Inscapes.mat', options.sample_dir, pat_select), 'saccades', 'fixations', 'fs')

end