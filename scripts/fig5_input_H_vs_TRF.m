%% Difference between VARX model input filter H and TRF

%% Define patients
model_dir = '../results/models_revision_1';

example_pat = 'NS127_02';

signal_type = 'HFA';

% Stimulus 'fixations', 'film_cuts', 'audio_env'
stim = 'fixations';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

fig_font = 20;

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

addpath('../src/Violinplot-Matlab-master')
addpath('../src')

load('../data/movie_subs_table.mat', 'movie_subs_table');

% Poster settings
poster_size = false;

if poster_size
    fig_font = 26;
end

%% Load data

% Initialize array to collect data
model_varx = cell(1, length(patients));
H_filter = cell(1, length(patients));
labels_all = cell(1, length(patients));

for pat = 1:length(patients)
    
    if strcmp(signal_type, 'LFP')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx', 'H', 'labels', 'fs_neural', 'stim_features')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx_hfa', 'H_hfa', 'labels', 'fs_neural', 'stim_features')
        m_varx = m_varx_hfa;
        H = H_hfa;
    end

    model_varx{pat} = m_varx;
    H_filter{pat} = H;
    labels_all{pat} = cellfun(@(C) sprintf('%s_%s', patients{pat}, C), labels, 'UniformOutput', false);

end

idx_stim = find(ismember(stim_features, stim));

%% Summarize data

% Find all channels with significant auditory response
B_pval = [];
labels_cat = [];
H_cat = [];
B_cat = [];

for pat = 1:length(patients)
    B_pval = [B_pval; model_varx{pat}.B_pval(:, idx_stim)];
    labels_cat = [labels_cat; labels_all{pat}];
    H_cat = [H_cat; H_filter{pat}(:,:,idx_stim)'];
    B_cat = [B_cat; model_varx{pat}.B(:,:,idx_stim)'];
end

% Find significant channels after FDR correction 
fdr = mafdr(B_pval, 'BHFDR', true);

idx_sig = fdr < 0.05;

labels_cat = labels_cat(idx_sig);
H_cat = H_cat(idx_sig, :);
B_cat = B_cat(idx_sig, :);
B_pval = B_pval(idx_sig);

%% Plot an example 
idx_example = ismember(patients, example_pat);

labels_sig_example = labels_cat(cellfun(@(C) contains(C, example_pat), labels_cat));
idx_sig_example = ismember(labels_all{idx_example}, labels_sig_example);

% if strcmp(example_pat, 'NS127_02') && idx_stim == 1 && strcmp(signal_type, 'LFP')
%     idx_sig_example(setdiff(1:length(idx_sig_example), [20:24, 87:92])) = 0;
% end

% Audio envelope
x_label_str = 'Delay [s]';
title_str = '';
out_file = sprintf('%s/fig5_trf_versus_B_audio_env_example_%s_%s_stim_%d.png', fig_dir, example_pat, signal_type, idx_stim);

plot_trf_b(model_varx{idx_example}, H_filter{idx_example}, fs_neural, idx_sig_example, idx_stim, x_label_str, title_str, poster_size, fig_font, out_file)

%% Measure length of responses
H_avg_pow = mean(H_cat.^2, 2);
B_avg_pow = mean(B_cat.^2, 2);

fwhm_H = resp_length(H_cat, fs_neural, 10, 5);
fwhm_B = resp_length(B_cat, fs_neural, 10, 5);

% Get patient labels
labels_pat = cell(1, length(labels_cat));

for l = 1:length(labels_cat)
    
    label_parts = strsplit(labels_cat{l}, '_');

    if length(label_parts) == 3
        labels_pat{l} = sprintf('%s_%s', label_parts{1}, label_parts{2});
    elseif length(label_parts) == 2
        labels_pat{l} = label_parts{1};
    end

end

% Compute power and length of responses across channels for each patient 
H_pow_pat = nan(1, length(patients));
B_pow_pat = nan(1, length(patients));
H_len_pat = nan(1, length(patients));
B_len_pat = nan(1, length(patients));

for p = 1:length(patients)

    idx_pat = ismember(labels_pat, patients{p});

    H_pow_pat(p) = mean(H_avg_pow(idx_pat));
    B_pow_pat(p) = mean(B_avg_pow(idx_pat));

    H_len_pat(p) = mean(fwhm_H(idx_pat)) / fs_neural * 1e3;
    B_len_pat(p) = mean(fwhm_B(idx_pat)) / fs_neural * 1e3;

end

%% Stats
B_pow_pat = B_pow_pat(~isnan(B_pow_pat));
H_pow_pat = H_pow_pat(~isnan(H_pow_pat));
B_len_pat = B_len_pat(~isnan(B_len_pat));
H_len_pat = H_len_pat(~isnan(H_len_pat));

p_pow = signrank(B_pow_pat, H_pow_pat);
p_len = signrank(B_len_pat, H_len_pat);

med_pow = median(B_pow_pat-H_pow_pat);
med_len = median(B_len_pat-H_len_pat);

fprintf('Power: median%c=%1.1e, p=%1.5f, N=%d\n', 916, med_pow, p_pow, length(B_pow_pat))
fprintf('Length: median%c=%2.2fms, p=%1.5f, N=%d\n', 916, med_len, p_len, length(B_len_pat))

%% Figures for difference of power and length of responses in each patient

if poster_size 
    figure('Units', 'inches', 'Position', [1,1,3.57,5]);
else
    figure('Position', [500,275,300,420])
end

hold on 

for i = 1:length(H_pow_pat)
    
    plot([1,2], [B_pow_pat(i), H_pow_pat(i)], 'k', 'LineWidth', 1.2)
    plot(1, B_pow_pat(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7)
    plot(2, H_pow_pat(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7)

end

xlim([0.75, 2.25])
xticks([1,2])
xticklabels({'B', 'H'})

ylabel('Average Power [a.u.]')
title('Power')

set(gca,'YAxisLocation','right')

grid on

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig5_trf_versus_B_power_across_patients_%s_stim_%d.png', ...
    fig_dir, signal_type, idx_stim), 'Resolution', 300)

%% Length of responses
if poster_size 
    figure('Units', 'inches', 'Position', [1,1,3.57,5]);
else
    figure('Position', [500,275,300,420])
end
hold on 

for i = 1:length(H_len_pat)
    
    plot([1,2], [B_len_pat(i), H_len_pat(i)], 'k', 'LineWidth', 1.2)
    plot(1, B_len_pat(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7)
    plot(2, H_len_pat(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7)

end

xlim([0.75, 2.25])
xticks([1,2])
xticklabels({'B', 'H'})

ylabel('Length of responses [ms]')
title('Length')

set(gca,'YAxisLocation','right')

grid on

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig5_trf_versus_B_fwhm_across_patients_%s_stim_%d.png', ...
    fig_dir, signal_type, idx_stim), 'Resolution', 300)

[~,p_fwhm,~,stats_fwhm] = ttest(B_len_pat, H_len_pat);
