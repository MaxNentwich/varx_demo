
%% Example of input filters 

% Select a patient
pat = 'NS127_02';

% Select signal type
signal_type = 'LFP';

% Name of inputs
input_names = {'Fixation Onset', 'Film Cuts', 'Audio Envelope'}; 
input_filenames = {'fixtion_onset', 'film_cuts', 'audio_envelope'}; 

% Font size
fig_font = 12;

% Figure directory
fig_dir = '../results/figures';

% Colos
cols_plot = ['r', 'g', 'b'];

%% Load data
if strcmp(signal_type, 'LFP')
    load(sprintf('../results/models/%s_varx_models.mat', pat), 'm_varx', 'labels', 'fs_neural')
elseif strcmp(signal_type, 'HFA')
    load(sprintf('../results/models/%s_varx_models.mat', pat), 'm_varx_hfa', 'labels', 'fs_neural')
    m_varx = m_varx_hfa;
end

% Find example channels (biggest effect size)
[~, ch_max] = max(m_varx.B_Rvalue);

% Time 
time = (1:size(m_varx.B,1)) / fs_neural;

n_feature = size(m_varx.B,3);

for f = 1:n_feature

    if f == n_feature
        figure('Units','inches', 'Position', [10,9,3,1.72])
    else
        figure('Units','inches', 'Position', [10,9,3,1.354])
    end

    plot(time, m_varx.B(:,ch_max(f),f), 'LineWidth', 1, 'Color', cols_plot(f))
    
    if f == n_feature
        xlabel('Delay to input [s]')
    else
        xticks([])
    end
    ylabel('Weight')
    title(input_names{f})

    fontsize(gcf, fig_font, 'points')

    exportgraphics(gcf, sprintf('%s/fig3_addon_filters_B_%s_%s.png', fig_dir, signal_type, input_filenames{f}), 'Resolution', 600)

end

%% Example of signals

% Define paths and parameters
data_dir = '/media/max/Workspace/Data/movies_legacy_sample/';
et_dir = '/media/max/e61df479-8f57-4855-b852-04ccdfb12a6c/ECoGData/Tobii/Patients';
cut_dir = '../data/scene_cuts_frames';
speech_dir = '../data/speech_files';
sample_data_dir = '/media/max/Workspace/Data/varx_data';

addpath(genpath('../src'))

stim_features = {'fixations'};
str = {'fixations'};

video_select = 'Despicable_Me_English.mat';

%% Compute all Models
data_pat_dir = sprintf('%s/%s', sample_data_dir, pat);

if exist(data_pat_dir, 'dir') == 0
    mkdir(data_pat_dir)
end

video = video_select;
[~,vid_name_legacy] = fileparts(video);

vid_file = [data_dir pat '/LFP/' video];

label_dir = sprintf('/media/max/e61df479-8f57-4855-b852-04ccdfb12a6c/ECoGData/Tobii/Patients/%s/matlab_data/Envelope_phase/BHA', pat);

load([data_dir pat '/LFP/' video],'lfp','fs');
load([data_dir pat '/BHA/' video],'envelope');
load([data_dir pat '/Stimuli/' video], stim_features{:});
load(sprintf('%s/%s', label_dir, video), 'labels')

x = [eval(stim_features{1})];
y = lfp; 
fs_neural = fs;
hfa = envelope;

time = (0:length(x)-1) / fs_neural;

%% Load film cuts
n_samples = length(time);

vid_name_legacy = strrep(vid_name_legacy, '_Rep_1', '');
vid_name_legacy = strrep(vid_name_legacy, '_Rep_2', '');

load(sprintf('%s/%s/matlab_data/Eyetracking/%s', et_dir, pat, video), 'eye')      
scene_table = readtable(sprintf('%s/%s_scenes.xlsx', cut_dir, vid_name_legacy));

if strfind(vid_name_legacy, 'Monkey') == 1
    scene_frame = scene_table.Frame * 2;
    scene_frame = scene_frame(scene_frame < length(eye.frame_time));
else
    scene_frame = scene_table.Frame;
end

scenes_time = eye.frame_time(scene_frame);
scenes_samp = round(interp1(time-time(1), 1:n_samples, scenes_time));

scenes_vec = zeros(1, n_samples);
scenes_vec(scenes_samp) = 1;

x = [x, scenes_vec'];

clearvars eye

%% Load audio
speech_vid = sprintf('%s/%s.mat', speech_dir, vid_name_legacy);

load(sprintf('%s/%s.mat', speech_dir, vid_name_legacy), 'audio', 'fs')
time_audio = (0:size(audio, 1)-1) / fs;

audio_env = abs(hilbert(audio));

audio_env = audio_env(time_audio < time(end));
time_audio = time_audio(time_audio < time(end));

audio_ds = resample(audio_env, fs_neural, fs);

if length(audio_ds) < length(x)
    audio_ds = [audio_ds, zeros(1, length(x) - length(audio_ds))];
elseif length(audio_ds) > length(x)
    audio_ds = audio_ds(1:length(audio_ds)-length(x));
end

x = [x, audio_ds'];

%% Plot
fig_font = 14;

% Fixation onset
figure('Units', 'inches', 'Position', [7,6,5.5,1.64])
plot(time(1:10*fs_neural), x(1:10*fs_neural,1), 'Color', cols_plot(1))
ylim([-0.1, 1.1])
xticks({})
yticks([0, 1])
title('Fixation onset')
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_addon_fixation_onset.png', fig_dir), 'Resolution', 600)

% Film cuts
figure('Units', 'inches', 'Position', [7,6,5.5,1.64])
plot(time(1:10*fs_neural), x(1:10*fs_neural,2), 'Color', cols_plot(2))
ylim([-0.1, 1.1])
xticks({})
yticks([0, 1])
title('Film cuts')
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_addon_film_cuts.png', fig_dir), 'Resolution', 600)

% Audio envelope
figure('Units', 'inches', 'Position', [7,6,5.5,2])
plot(time(1:10*fs_neural), x(1:10*fs_neural,3), 'Color', cols_plot(3))
xlabel('Time [s]')
title('Audio envelope')
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_addon_audio_envelope.png', fig_dir), 'Resolution', 600)

% Neural signal
y = y./std(y);

figure('Units', 'inches', 'Position', [7,6,5.5,1.59])
plot(time(1:10*fs_neural), y(1:10*fs_neural,15), 'Color', 'k')
xticks([])
title('LFP')
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_addon_LFP.png', fig_dir), 'Resolution', 600)

% HFA
hfa = hfa./std(hfa);

figure('Units', 'inches', 'Position', [7,6,5.5,2])
plot(time(1:10*fs_neural), hfa(1:10*fs_neural,15), 'Color', 'k')
xlabel('Time [s]')
title('BHA')
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_addon_HFA.png', fig_dir), 'Resolution', 600)
