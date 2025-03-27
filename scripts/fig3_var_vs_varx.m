%% Difference between VAR and VARX models

%% Define patients
model_dir = '../results/models_revision_1';

example_pat = 'NS127_02';

signal_type = 'HFA';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

fig_font = 20;

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

addpath('../src/Violinplot-Matlab-master')

p_thresh = 0.0001;

select_sig = false;

% Plot p-values ('pval') or R-values ('rval') for example
example_val = 'pval';

% Poster settings
poster_size = false;

if poster_size
    fig_font = 26;
end

%% Load data

% Initialize array to collect data
model_var = cell(1, length(patients));
model_varx = cell(1, length(patients));
var_Rvalue = cell(1, length(patients));
varx_Rvalue = cell(1, length(patients));
n_sig_var = nan(1, length(patients));
n_sig_varx = nan(1, length(patients));
r_diff_mean = nan(1, length(patients));
R_diff_features = cell(1, length(patients));

for pat = 1:length(patients)
    
    if strcmp(signal_type, 'LFP')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx', 'm_varx_additive', 'stim_features')
        idx_none = ismember(stim_features, 'None');
        m_var = m_varx_additive{1};
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx_hfa', 'm_varx_additive_hfa', 'stim_features')
        idx_none = ismember(stim_features, 'None');
        m_var = m_varx_additive_hfa{1};
        m_varx = m_varx_hfa;
        m_varx_features = m_varx_additive_hfa;
    end

    model_var{pat} = m_var;
    model_varx{pat} = m_varx;

    model_var{pat}.A_pval(model_var{pat}.A_pval == 0) = 1e-10;
    model_varx{pat}.A_pval(model_varx{pat}.A_pval == 0) = 1e-10;

    % Number of significant connections
    n_ch = size(model_var{pat}.A_pval,1);
    n_conn = n_ch^2 - n_ch;

    n_sig_var(pat) = (sum(model_var{pat}.A_pval(:) < p_thresh) - n_ch) / n_conn;
    n_sig_varx(pat) = (sum(model_varx{pat}.A_pval(:) < p_thresh) - n_ch) / n_conn;

    %% Effect size
    var_Rvalue{pat} = sqrt(1-exp(-model_var{pat}.A_Deviance/model_var{pat}.T));
    varx_Rvalue{pat} = sqrt(1-exp(-model_varx{pat}.A_Deviance/model_varx{pat}.T));

    % R values of significant connections
    if select_sig
        idx_sig = (model_var{pat}.A_pval < p_thresh) & (model_varx{pat}.A_pval < p_thresh);
        idx_sig = idx_sig - eye(size(idx_sig));
    else
        idx_sig = eye(size(var_Rvalue{pat})) == 0;
    end

    r_diff_mean(pat) = mean(varx_Rvalue{pat}(idx_sig==1) - var_Rvalue{pat}(idx_sig==1));


end

%% Difference of significant cannels and effect size across patients
n_sig_diff = n_sig_varx - n_sig_var;

p_chans = signrank(n_sig_diff);
p_rval = signrank(real(r_diff_mean));

fprintf('Ratio of significant channels: median=%1.1e p=%1.5f, N=%d\n', ...
    median(n_sig_diff), p_chans, length(n_sig_diff))
fprintf('Effect size: median=%1.1e p=%1.5f, N=%d\n', ...
    median(real(r_diff_mean)), p_rval, length(r_diff_mean))

%% Plot an example
idx_example = ismember(patients, example_pat);

if strcmp(example_val, 'pval')
    mat_var = log10(model_var{idx_example}.A_pval);
    mat_varx = log10(model_varx{idx_example}.A_pval);
    cmap = 'summer';
elseif strcmp(example_val, 'rval')
    mat_var = var_Rvalue{idx_example};
    mat_varx = varx_Rvalue{idx_example};
    cmap = 'Reds';
end

plot_range_ar = [min(mat_var(:)), max(mat_var(:))];
plot_range = [min(mat_varx(:)), max(mat_varx(:))];
plot_range = min([plot_range_ar; plot_range]);

if strcmp(example_val, 'rval')
    plot_range(1) = 0;
    plot_range(2) = plot_range(2)*0.1;
end

plot_diff = max(abs(mat_var(:) - mat_varx(:)));
plot_range_diff = [-plot_diff, plot_diff];

if poster_size 
    fig1 = figure('Units', 'inches', 'Position', [1,1,20,5]);
else
    fig1 = figure('Position', [3000,350,1450,450]);
end

t = tiledlayout(1,3);

ax1 = nexttile;

imagesc(mat_var)

clim(ax1, plot_range)
axis square
colormap(ax1, slanCM(cmap))
xlabel('Channels')
ylabel('Channels')
title('No Inputs')

ax2 = nexttile;

imagesc(mat_varx)

clim(ax2, plot_range)
axis square
colormap(ax2, slanCM(cmap))
cb = colorbar(); 
if strcmp(example_val, 'pval')
    ylabel(cb,'log p-value','Rotation',90)
elseif strcmp(example_val, 'rval')
    ylabel(cb,'R','Rotation',90)
end
xlabel('Channels')
yticks([])
title('Inputs')

ax3 = nexttile;

imagesc(mat_var - mat_varx)

clim(ax3, 0.05*plot_range_diff)

axis square
colormap(ax3, slanCM('bwr'))
cb = colorbar(); 
ylabel(cb,['\Delta log p-value' newline '(No Input - Input)'] ,'Rotation',90)
xlabel('Channels')
yticks([])
title('Difference')

fontsize(fig1, fig_font, 'points')

exportgraphics(fig1, sprintf('%s/fig3_varx_full_none_example_%s_%s.png', ...
    fig_dir, example_pat, signal_type), 'Resolution', 600)

%% Difference of significant connections and R values over all paitents
R_varx = cellfun(@(C) abs(mean(C(:))), varx_Rvalue);
R_var = cellfun(@(C) abs(mean(C(:))), var_Rvalue);

R_plot = [R_var; R_varx];
ch_plot = [n_sig_var; n_sig_varx];

R_plot = R_plot - R_plot(1,:);
ch_plot = ch_plot - ch_plot(1,:);

figure('Position', [1000,750,750,450])

% Significant Channels
tiledlayout(1,10);

nexttile(1,[1,4])
hold on

plot(ch_plot(:, diff(ch_plot) > 0), '.-', 'Color', [0.75, 0.3, 0], 'LineWidth', 2, 'MarkerSize', 20)
plot(ch_plot(:, diff(ch_plot) <= 0), '.-', 'Color', [0, 0.3, 0.75], 'LineWidth', 2, 'MarkerSize', 20)

xticks([1, 2])
xticklabels({'No Input', 'Input'})
xtickangle(45)
xlim([0.75, 2.25])

ylabel('Fraction of Channels')
set(gca, 'YAxisLocation', 'right')
ylim([1.1*min(ch_plot(:)), 1.1*max(ch_plot(:))])

fontsize(fig_font, 'Points')
title(['Significant' newline 'Connections'])

legend(['Increase', repmat({''}, 1, length(ch_plot)-2), 'Decrease'], ...
    'Position', [0.77,0.02,0.22,0.15]);

grid on

% Effect size
nexttile(5,[1,4])
hold on 

plot(R_plot(:, diff(R_plot) > 0), '.-', 'Color', [0.75, 0.3, 0], 'LineWidth', 2, 'MarkerSize', 20)
plot(R_plot(:, diff(R_plot) < 0), '.-', 'Color', [0, 0.3, 0.75], 'LineWidth', 2, 'MarkerSize', 20)

xticks([1, 2])
xticklabels({'No Input', 'Input'})
xtickangle(45)
xlim([0.75, 2.25])

ylabel('R')
set(gca, 'YAxisLocation', 'right')
ylim([1.1*min(R_plot(:)), 1.05*max(R_plot(:))])

fontsize(fig_font, 'Points')
title('Effect Size')

grid on

exportgraphics(gcf, sprintf('%s/fig3_varx_full_none_sig_ch_diff_R_%s_sig_ch_%d.png', ...
    fig_dir, signal_type, select_sig), 'Resolution', 600)

%% 
% if poster_size 
%     figure('Units', 'inches', 'Position', [1,1,8.5,5]);
% else
%     figure('Position', [400,300,725,450])
% end
% 
% tiledlayout(1,9);
% 
% % Violin plot of differences for one paitent
% nexttile(1,[1,4])
% hold on
% 
% scatter(0.1*randn(1, length(n_sig_diff)), n_sig_diff, 'k', 'filled')
% 
% plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
% plot([-0.25, 0.25], median(n_sig_diff)*ones(2,1), 'k--')
% 
% grid on
% box on
% xticks([])
% xlim([-0.3, 0.3])
% ylim_abs =  1.2 * max(abs(n_sig_diff));
% ylim([-ylim_abs, ylim_abs])
% set(gca, 'YAxisLocation', 'right')
% ylabel(['\DeltaRatio' newline '(Input - No Input)'])
% 
% title(['Significant' newline 'Connections'])
% 
% ax = ancestor(gca, 'axes');
% ax.YAxis.Exponent = 0;
% if strcmp(signal_type, 'LFP')
%     ytickformat('%0.3f')
% elseif strcmp(signal_type, 'HFA')
%     ytickformat('%0.4f')
% end
% 
% %% Effect size on same panel 
% nexttile(6,[1,4])
% hold on 
% 
% scatter(0.1*randn(1, length(r_diff_mean)), r_diff_mean, 'k', 'filled')
% plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
% plot([-0.25, 0.25], median(r_diff_mean)*ones(2,1), 'k--')
% ylim_abs =  1.2 * max(abs(r_diff_mean));
% 
% grid on
% box on
% xticks([])
% xlim([-0.3, 0.3])
% ylim([-ylim_abs, ylim_abs])
% set(gca, 'YAxisLocation', 'right')
% ylabel(['Mean \DeltaR' newline '(Input - No Input)'])
% 
% title(['Effect' newline 'size'])
% 
% fontsize(gcf, fig_font, 'points')
