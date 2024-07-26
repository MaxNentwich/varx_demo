%% Difference between VARX model input filter H and TRF

%% Define patients
example_pat = 'NS127_02';

signal_type = 'LFP';

% Stimulus 3 -> audio; 1 -> fixations; 2 -> cuts
idx_stim = 3;

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

fig_font = 16;

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

addpath('../src/Violinplot-Matlab-master')
addpath('../src')

load('../data/movie_subs_table.mat', 'movie_subs_table');

%% Load data

% Initialize array to collect data
model_varx = cell(1, length(patients));
H_filter = cell(1, length(patients));
labels_all = cell(1, length(patients));

for pat = 1:length(patients)
    
    if strcmp(signal_type, 'LFP')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx', 'H', 'labels', 'fs_neural')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx_hfa', 'H_hfa', 'labels', 'fs_neural')
        m_varx = m_varx_hfa;
        H = H_hfa;
    end

    model_varx{pat} = m_varx;
    H_filter{pat} = H;
    labels_all{pat} = cellfun(@(C) sprintf('%s_%s', patients{pat}, C), labels, 'UniformOutput', false);

end

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

% Audio envelope
x_label_str = 'Delay [s]';
title_str = '';
out_file = sprintf('%s/fig5_trf_versus_B_audio_env_example_%s_%s_stim_%d.png', fig_dir, example_pat, signal_type, idx_stim);

plot_trf_b(model_varx{idx_example}, H_filter{idx_example}, fs_neural, idx_sig_example, idx_stim, x_label_str, title_str, out_file)

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
[~, p_pow, ~, stat_pow] = ttest(B_pow_pat, H_pow_pat);
[~, p_len, ~, stat_len] = ttest(B_len_pat, H_len_pat);

fprintf('Power: t(%d)=%1.2f, p=%1.5f\n', stat_pow.df, stat_pow.tstat, p_pow)
fprintf('length: t(%d)=%1.2f, p=%1.5f\n', stat_len.df, stat_len.tstat, p_len)

%% Figures for difference of power and length of responses in each patient

figure('Position', [500,275,300,420])
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

grid on

fontsize(16, 'points')

exportgraphics(gcf, sprintf('%s/fig5_trf_versus_B_power_across_patients_%s_stim_%d.png', ...
    fig_dir, signal_type, idx_stim), 'Resolution', 300)

%% Length of responses
figure('Position', [500,275,300,420])
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

grid on

fontsize(16, 'points')

exportgraphics(gcf, sprintf('%s/fig5_trf_versus_B_fwhm_across_patients_%s_stim_%d.png', ...
    fig_dir, signal_type, idx_stim), 'Resolution', 300)

[~,p_fwhm,~,stats_fwhm] = ttest(B_len_pat, H_len_pat);
