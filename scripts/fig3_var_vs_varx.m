%% Difference between VAR and VARX models

%% Define patients
example_pat = 'NS127_02';

signal_type = 'LFP';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

fig_font = 16;

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

addpath('../src/Violinplot-Matlab-master')

p_thresh = 0.0001;

select_sig = false;

%% Load data

% Initialize array to collect data
model_var = cell(1, length(patients));
model_varx = cell(1, length(patients));
n_sig_var = nan(1, length(patients));
n_sig_varx = nan(1, length(patients));
r_diff_mean = nan(1, length(patients));
R_diff_features = cell(1, length(patients));

for pat = 1:length(patients)
    
    if strcmp(signal_type, 'LFP')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx', 'm_varx_features', 'feature_combos')
        idx_none = ismember(feature_combos, 'None');
        m_var = m_varx_features{idx_none};
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx_hfa', 'm_varx_hfa_features', 'feature_combos')
        idx_none = ismember(feature_combos, 'None');
        m_var = m_varx_hfa_features{idx_none};
        m_varx = m_varx_hfa;
        m_varx_features = m_varx_hfa_features;
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
    var_Rvalue = sqrt(1-exp(-model_var{pat}.A_Deviance/model_var{pat}.T));
    varx_Rvalue = sqrt(1-exp(-model_varx{pat}.A_Deviance/model_varx{pat}.T));

    % R values of significant connections
    if select_sig
        idx_sig = (model_var{pat}.A_pval < p_thresh) & (model_varx{pat}.A_pval < p_thresh);
        idx_sig = idx_sig - eye(size(idx_sig));
    else
        idx_sig = eye(size(var_Rvalue)) == 0;
    end

    r_diff_mean(pat) = mean(var_Rvalue(idx_sig==1) - varx_Rvalue(idx_sig==1));

    %% Coordinates

    %% Compute R values for models with different features
    if length(feature_combos) == 0, continue, end

    R_varx_feat = NaN(size(varx_Rvalue,1), size(varx_Rvalue,2), length(feature_combos));
    for feat = 1:length(feature_combos)
        R_varx_feat(:,:,feat) = abs(sqrt(1-exp(-m_varx_features{feat}.A_Deviance/m_varx_features{feat}.T)));  
    end
    
    if select_sig

        idx_sig_varx = m_varx.A_pval < p_thresh;
        idx_sig_varx = idx_sig_varx - (eye(size(idx_sig_varx)) .* diag(idx_sig_varx));
    
        idx_sig_feat = NaN(size(varx_Rvalue,1), size(varx_Rvalue,2), length(feature_combos));
        for feat = 1:length(feature_combos)
            idx_sig_feat(:,:,feat) = m_varx_features{feat}.A_pval < p_thresh;
            idx_sig_feat(:,:,feat) = idx_sig_feat(:,:,feat) - (eye(size(idx_sig_feat(:,:,feat))) .* diag(idx_sig_feat(:,:,feat)));
        end
    
        idx_sig_feat = (sum(idx_sig_feat,3) == length(feature_combos)) & idx_sig_varx;

    else
        idx_sig_feat = eye(size(varx_Rvalue)) == 0;
    end

    idx_none = ismember(feature_combos, 'None');
    R_varx_none = R_varx_feat(:,:,idx_none);
    
    R_varx_feat(:,:,idx_none) = [];
    R_varx_feat = cat(3, varx_Rvalue, R_varx_feat);

    feature_combos(idx_none) = [];
    C =  {{'Fixation, Cuts and Sound'}, feature_combos};
    feature_combos = cat(2, C{:});

    R_diff_vec = NaN(sum(idx_sig_feat(:)), length(feature_combos));
    for feat = 1:length(feature_combos)
        R_feat = R_varx_feat(:,:,feat);
        R_diff_vec(:,feat) = R_varx_none(idx_sig_feat) - R_feat(idx_sig_feat);     
    end

    R_diff_features{pat} = mean(R_diff_vec);

end

%% Difference of significant cannels and effect size across patients
n_sig_diff = n_sig_var - n_sig_varx;
[~,p_chans,~,stats_chans] = ttest(n_sig_diff);

[~,p_rval,~,stats_rval] = ttest(abs(r_diff_mean));

% Models with different features
R_feat_mat = abs(cell2mat(R_diff_features));
R_feat_mat = reshape(R_feat_mat, [length(feature_combos), length(R_diff_features)]);

ci = zeros(size(R_feat_mat,1),2);
for i = 1:size(R_feat_mat,1)
    [~,~,ci(i,:),~] = ttest(R_feat_mat(i,:));
end

% Test difference of individual features
[~,p_fix_cuts,~,stats_fix_cuts] = ttest(R_feat_mat(ismember(feature_combos, 'Fixations'), :),  ...
    R_feat_mat(ismember(feature_combos, 'Cuts'), :));

[~,p_fix_sound,~,stats_fix_sound] = ttest(R_feat_mat(ismember(feature_combos, 'Fixations'), :),  ...
    R_feat_mat(ismember(feature_combos, 'Sound'), :));

[~,p_cuts_sound,~,stats_cuts_sound] = ttest(R_feat_mat(ismember(feature_combos, 'Cuts'), :),  ...
    R_feat_mat(ismember(feature_combos, 'Sound'), :));

[~,p_all_sound,~,stats_all_sound] = ttest(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :),  ...
    R_feat_mat(ismember(feature_combos, 'Sound'), :));

[~,p_all_cuts,~,stats_all_cuts] = ttest(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :),  ...
    R_feat_mat(ismember(feature_combos, 'Cuts'), :));

[~,p_all_fix,~,stats_all_fix] = ttest(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :),  ...
    R_feat_mat(ismember(feature_combos, 'Fixations'), :));

%% Plot an example
idx_example = ismember(patients, example_pat);

plot_range_ar = [min(log10(model_var{idx_example}.A_pval(:))), max(log10(model_var{idx_example}.A_pval(:)))];
plot_range = [min(log10(model_varx{idx_example}.A_pval(:))), max(log10(model_varx{idx_example}.A_pval(:)))];
plot_range = min([plot_range_ar; plot_range]);

plot_diff = max(abs(log10(model_var{idx_example}.A_pval(:)) - log10(model_varx{idx_example}.A_pval(:))));
plot_range_diff = [-plot_diff, plot_diff];

fig1 = figure('Position', [3000,350,1450,450]);

t = tiledlayout(1,3);

ax1 = nexttile;

imagesc(log10(model_var{idx_example}.A_pval))

clim(ax1, plot_range)
axis square
colormap(ax1, slanCM('summer'))
xlabel('Channels')
ylabel('Channels')
title('No Inputs')

ax2 = nexttile;

imagesc(log10(model_varx{idx_example}.A_pval))

clim(ax2, plot_range)
axis square
colormap(ax2, slanCM('summer'))
cb = colorbar(); 
ylabel(cb,'log p-value','Rotation',90)
xlabel('Channels')
yticks([])
title('Inputs')

ax3 = nexttile;

imagesc(log10(model_var{idx_example}.A_pval) - log10(model_varx{idx_example}.A_pval))

clim(ax3, 0.05*plot_range_diff)

axis square
colormap(ax3, slanCM('bwr'))
cb = colorbar(); 
ylabel(cb,'\Delta log p-value (No Input - Input)' ,'Rotation',90)
xlabel('Channels')
yticks([])
title('Difference')

fontsize(fig1, fig_font, 'points')

exportgraphics(fig1, sprintf('%s/fig3_varx_full_none_example_%s_%s.png', ...
    fig_dir, example_pat, signal_type), 'Resolution', 300)

%% Difference of significant connections and R values over all paitents

figure('Position', [400,300,725,450])
tiledlayout(1,9);

% Violin plot of differences for one paitent
nexttile(1,[1,4])
hold on

scatter(0.1*randn(1, length(n_sig_diff)), n_sig_diff, 'k', 'filled')

plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
plot([-0.25, 0.25], median(n_sig_diff)*ones(2,1), 'k--')

grid on
box on
xticks([])
xlim([-0.3, 0.3])
ylim_abs =  1.2 * max(abs(n_sig_diff));
ylim([-ylim_abs, ylim_abs])
set(gca, 'YAxisLocation', 'right')
ylabel('\Delta Ratio (No Input - Input)')

title('Sig. Connections')

ax = ancestor(gca, 'axes');
ax.YAxis.Exponent = 0;
if strcmp(signal_type, 'LFP')
    ytickformat('%0.3f')
elseif strcmp(signal_type, 'HFA')
    ytickformat('%0.4f')
end

%% Effect size on same panel 
nexttile(6,[1,4])
hold on 

scatter(0.1*randn(1, length(r_diff_mean)), r_diff_mean, 'k', 'filled')
plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
plot([-0.25, 0.25], median(r_diff_mean)*ones(2,1), 'k--')
ylim_abs =  1.2 * max(abs(r_diff_mean));

grid on
box on
xticks([])
xlim([-0.3, 0.3])
ylim([-ylim_abs, ylim_abs])
set(gca, 'YAxisLocation', 'right')
ylabel('Mean \DeltaR (No Input - Input)')

if strcmp(signal_type, 'LFP')
    ax = ancestor(gca, 'axes');
    ax.YAxis.Exponent = 0;
    ytickformat('%0.4f')
end

title('Effect size')

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_varx_full_none_sig_ch_diff_R_%s_sig_ch_%d.png', ...
    fig_dir, signal_type, select_sig), 'Resolution', 300)

%% Features
figure('Position', [725,700,700,525])
hold on

plot(flipud(mean(R_feat_mat,2)), 'r*-', 'LineWidth', 2)
plot(flipud(ci(:,1)), 'k*-')
plot(flipud(ci(:,2)), 'k*-')

xticklabels(fliplr(feature_combos))
xtickangle(45)

ylabel('Mean \Delta R-value')
title('No Input - Input')
grid('on')

fontsize(18, 'points')

exportgraphics(gcf, sprintf('%s/fig3_R_features_%s_sig_ch_%d.png', ...
    fig_dir, signal_type, select_sig), 'Resolution', 300)
