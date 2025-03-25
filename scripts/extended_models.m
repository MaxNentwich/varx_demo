
model_dir = '/media/max/Workspace/Code/varx_demo/results/models_revision_1';
fig_dir = '/media/max/Workspace/Code/varx_demo/results/figures';

signal_type = 'HFA';

addpath(genpath('./data_prep'))

m_files = dir(model_dir);
m_files([m_files.isdir]) = [];

%% Effect of inividual regressors
fd = cell(length(m_files),1);

for m = 1:length(m_files)

    if strcmp(signal_type, 'LFP')
        load(sprintf('%s/%s', model_dir, m_files(m).name), 'm_varx', 'm_varx_features', 'stim_features')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('%s/%s', model_dir, m_files(m).name), 'm_varx', 'm_varx_features_hfa', 'stim_features')
        m_varx_features = m_varx_features_hfa;
    end

    fd{m} = zeros(1,length(stim_features)); 
    
    for f = 1:length(stim_features)
        fd{m}(f) = mean(mean((m_varx_features{f}.A_Rvalue(m_varx.B_pval(:,f) < 0.05, :) - m_varx.A_Rvalue(m_varx.B_pval(:,f) < 0.05, :)))); 
    end

end

fd = cell2mat(fd);

% Stats
p_fd = zeros(1, size(fd,2));

for jj = 1:size(fd, 2)
    p_fd(jj) = signrank(fd(:,jj));
end

p_fd = mafdr(p_fd, 'BHFDR', true);

% Summarize
mean_fd = nanmean(fd);
err_fd = nanstd(fd) / sqrt(size(fd,1));

[~, idx_sort] = sort(mean_fd, 'descend');

mean_fd = mean_fd(idx_sort);
err_fd = err_fd(idx_sort);
stim_features = stim_features(idx_sort);
p_fd = p_fd(idx_sort);

% Figure
max_lim = max(mean_fd + err_fd);

x_axis = 1:length(stim_features);
x_vector = [x_axis, fliplr(x_axis)];

figure
patch = fill(x_vector, [mean_fd+err_fd,fliplr(mean_fd-err_fd)], 'k');

set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.2);
hold on;
plot(x_axis, mean_fd, 'k.-', 'MarkerSize', 20, 'LineWidth', 2);

for i = 1:length(p_fd)
    if p_fd(i) < 0.05
        plot(i, 1.3*max_lim, 'k*', 'MarkerSize', 10)
    end
    fprintf('Removing %s: DeltaR = %1.2e, p = %1.5f\n', stim_features{i}, mean_fd(:,i), p_fd(i))
end

xticks(1:length(stim_features))
xticklabels(cellfun(@(C) strrep(C, '_', ' '), stim_features, 'UniformOutput', false))
xtickangle(60)

xlim([0.8, length(stim_features)+0.2])

ylabel('\DeltaR')

title('Removing features')

fontsize(20, 'points')

grid on

exportgraphics(gcf, sprintf('%s/revision_fig3_removing_features.png', fig_dir))

%% Decrease of R when adding features progressively 
R_diff = cell(length(m_files),1);

for m = 1:length(m_files)

    if strcmp(signal_type, 'LFP')
        load(sprintf('%s/%s', model_dir, m_files(m).name), 'm_varx_additive', 'features_sort')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('%s/%s', model_dir, m_files(m).name), 'm_varx_additive_hfa', 'features_sort')
        m_varx_additive = m_varx_additive_hfa;
    end

    R_zero = m_varx_additive{1}.A_Rvalue;
    R_zero = abs(R_zero(eye(size(R_zero)) ~= 1));
    
    R_diff{m} = zeros(1,length(features_sort)); 
    
    for f = 1:length(features_sort)
        R_feat = m_varx_additive{f+1}.A_Rvalue;
        R_feat = abs(R_feat(eye(size(R_feat)) ~= 1));    
        R_diff{m}(f) = mean(R_feat - R_zero);
    end

end

R_diff = cell2mat(R_diff);

% Test which features change compared to the previous ones
p_diff = zeros(1, size(R_diff,2));
R_step = zeros(1, size(R_diff,2));

p_diff(1) = signrank(R_diff(:,1));

for j = 1:size(R_diff,2)-1
    R_step(j+1) = mean(R_diff(:,j+1)) - mean(R_diff(:,j));
    p_diff(j+1) = ranksum(R_diff(:,j), R_diff(:,j+1));
end

p_diff_corr = mafdr(p_diff, 'BHFDR', true);

% Summarize
mean_RD = mean(R_diff);
err_RD = std(R_diff) / sqrt(size(R_diff,1));

% Plot
min_lim = round(min(mean_RD - err_RD), 6);

x_axis = 1:length(features_sort);
x_vector = [x_axis, fliplr(x_axis)];

figure
patch = fill(x_vector, [mean_RD+err_RD,fliplr(mean_RD-err_RD)], 'k');

set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', 0.2);
hold on;
plot(x_axis, mean_RD, 'k.-', 'MarkerSize', 20, 'LineWidth', 2);

for s = 1:length(p_diff_corr)
    if p_diff_corr(s) < 0.05
        plot(s, 1.1*min_lim, 'k*', 'MarkerSize', 10)
        fprintf('Adding %s: DeltaR = %1.2e, R step = %1.2e, p = %1.5f\n', features_sort{s}, mean_RD(:,s), R_step(s), p_diff_corr(s))
    end
end

xticks(1:length(features_sort))
xticklabels(cellfun(@(C) strrep(C, '_', ' '), features_sort, 'UniformOutput', false))
xtickangle(60)

xlim([0.8, length(features_sort)+0.2])

ylabel('\DeltaR')

title('Adding features')

fontsize(20, 'points')

grid on

exportgraphics(gcf, sprintf('%s/revision_fig3_adding_features.png', fig_dir))


%% Effect of removing regressor on B
load(sprintf('%s/%s', model_dir, m_files(7).name), ...
    'm_varx_ch_cut', 'ch_keep', 'm_varx', 'm_varx_example', ...
    'example_features', 'stim_features', 'fs_neural')

idx_env = ismember(stim_features, 'audio_env');
idx_edge = ismember(stim_features, 'acoustic_edges');
idx_fix = ismember(stim_features, 'fixations');

label_1 = strjoin(cellfun(@(C) strrep(C, '_', ' '), example_features{1}, 'UniformOutput', false), ' + ');
label_2 = strjoin(cellfun(@(C) strrep(C, '_', ' '), example_features{2}, 'UniformOutput', false), ' + ');
label_3 = strjoin(cellfun(@(C) strrep(C, '_', ' '), example_features{3}, 'UniformOutput', false), ' + ');

idx_sig = m_varx_example{3}.B_pval(:,idx_env) < 0.0001;

B_ex1 = m_varx_example{1}.B_coeff(:, idx_sig, idx_env)';
B_ex2 = m_varx_example{2}.B_coeff(:, idx_sig, idx_env)';
B_ex3 = m_varx_example{3}.B_coeff(:, idx_sig, idx_env)';

B_ex3_sac = m_varx_example{3}.B_coeff(:, idx_sig, idx_fix)';
B_ex3_edge = m_varx_example{3}.B_coeff(:, idx_sig, idx_edge)';

time = (0:size(B_ex1,2)-1) / fs_neural;

% Channel in model with removed channels
ch_id = 1:length(idx_sig);
idx_sig_cut = ismember(ch_keep, ch_id(idx_sig));

B_ex_chns_all = m_varx.B_coeff(:, idx_sig, idx_env)';
B_ex_chns_cut = m_varx_ch_cut.B_coeff(:, idx_sig_cut, idx_env)';

%% Line plot of an example
ch_example = 12;

% Edge response in a model with all channels compared to a model with some 
% channels that were removed
figure('Position', [500,750,475,300])
hold on

plot(time, LpFilter(B_ex_chns_all(ch_example,:), 5, 10, 60), 'Color', [0.9, 0.5, 0], 'LineWidth', 2)
plot(time, LpFilter(B_ex_chns_cut(ch_example,:), 5, 10, 60), 'Color', [0, 0.7, 0.8], 'LineStyle', '--', 'LineWidth', 2)

ax = gca();
ax.FontSize = 16;

ylabel('B - Envelope')
xlabel('Time [s]')
title('Removing Channels')
grid on 

legend({'All Channels', 'Channels removed'}, 'Position', [0.62,0.25,0.3,0.2], ...
    'FontSize', 12)

exportgraphics(gcf, sprintf('%s/revision_supp1_remove_channels.png', fig_dir))

% Edge respone for models with correlated and uncorrelated input features
figure('Position', [1250,750,550,300])
hold on

plot(time, LpFilter(B_ex1(ch_example,:), 5, 10, 60), 'Color', [0.9, 0.01, 0.15], 'LineWidth', 2)
plot(time, LpFilter(B_ex2(ch_example,:), 5, 10, 60), 'Color', [0.1, 0.1, 0.9], 'LineWidth', 2)
plot(time, LpFilter(B_ex3(ch_example,:), 5, 10, 60), 'Color', [0.02, 0.8, 0.05], 'LineStyle', '--', 'LineWidth', 2)

ylabel('B - Envelope')
xlabel('Time [s]')
title('Adding Inputs')
grid on 

ax = gca();
ax.OuterPosition = [0 0 0.8 1];
ax.FontSize = 16;

legend({'Envelope', 'Envelope + Edges', 'Envelope + Edges + Fixations'}, ...
    'Position', [0.47,0.25,0.5,0.2], 'FontSize', 12)

exportgraphics(gcf, sprintf('%s/revision_supp1_correlated_features.png', fig_dir))

% All responses in the full model
figure('Position', [1700,750,475,300])
hold on

plot(time, LpFilter(B_ex3(ch_example,:), 5, 10, 60), 'k', 'LineWidth', 2)
plot(time, LpFilter(B_ex3_edge(ch_example,:), 5, 10, 60) - 0.02, 'k', 'LineWidth', 2)
plot(time, LpFilter(B_ex3_sac(ch_example,:), 5, 10, 60) - 0.04, 'k', 'LineWidth', 2)

yticks([-0.04, -0.02, 0])
ylim([-0.05, 0.01])
yticklabels({'fixation onset', 'acoustic edges', 'auditory envelope'})
ylabel('Feature')
xlabel('Time [s]')
title('Responses for full model')
grid on 

fontsize(16, 'points')

exportgraphics(gcf, sprintf('%s/revision_supp1_correlated_features_full_model.png', fig_dir))
