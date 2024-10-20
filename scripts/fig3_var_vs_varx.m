%% Difference between VAR and VARX models

%% Define patients
example_pat = 'NS127_02';

signal_type = 'HFA';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

fig_font = 16;

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

    %% Coordinates

    %% Compute R values for models with different features
    if length(feature_combos) == 0, continue, end

    R_varx_feat = NaN(size(varx_Rvalue{pat},1), size(varx_Rvalue{pat},2), length(feature_combos));
    for feat = 1:length(feature_combos)
        R_varx_feat(:,:,feat) = abs(sqrt(1-exp(-m_varx_features{feat}.A_Deviance/m_varx_features{feat}.T)));  
    end
    
    if select_sig

        idx_sig_varx = m_varx.A_pval < p_thresh;
        idx_sig_varx = idx_sig_varx - (eye(size(idx_sig_varx)) .* diag(idx_sig_varx));
    
        idx_sig_feat = NaN(size(varx_Rvalue{pat},1), size(varx_Rvalue{pat},2), length(feature_combos));
        for feat = 1:length(feature_combos)
            idx_sig_feat(:,:,feat) = m_varx_features{feat}.A_pval < p_thresh;
            idx_sig_feat(:,:,feat) = idx_sig_feat(:,:,feat) - (eye(size(idx_sig_feat(:,:,feat))) .* diag(idx_sig_feat(:,:,feat)));
        end
    
        idx_sig_feat = (sum(idx_sig_feat,3) == length(feature_combos)) & idx_sig_varx;

    else
        idx_sig_feat = eye(size(varx_Rvalue{pat})) == 0;
    end

    idx_none = ismember(feature_combos, 'None');
    R_varx_none = R_varx_feat(:,:,idx_none);
    
    R_varx_feat(:,:,idx_none) = [];
    R_varx_feat = cat(3, varx_Rvalue{pat}, R_varx_feat);

    feature_combos(idx_none) = [];
    C =  {{'Fixation, Cuts and Sound'}, feature_combos};
    feature_combos = cat(2, C{:});

    R_diff_vec = NaN(sum(idx_sig_feat(:)), length(feature_combos));
    for feat = 1:length(feature_combos)
        R_feat = R_varx_feat(:,:,feat);
        R_diff_vec(:,feat) = R_feat(idx_sig_feat) - R_varx_none(idx_sig_feat);     
    end

    R_diff_features{pat} = mean(R_diff_vec);

end

%% Difference of significant cannels and effect size across patients
n_sig_diff = n_sig_varx - n_sig_var;

p_chans = signrank(n_sig_diff);
p_rval = signrank(real(r_diff_mean));

fprintf('Ratio of significant channels: median=%1.1e p=%1.5f, N=%d\n', ...
    median(n_sig_diff), p_chans, length(n_sig_diff))
fprintf('Effect size: median=%1.1e p=%1.5f, N=%d\n', ...
    median(real(r_diff_mean)), p_rval, length(r_diff_mean))

% Models with different features
R_feat_mat = real(cell2mat(R_diff_features));
R_feat_mat = reshape(R_feat_mat, [length(feature_combos), length(R_diff_features)]);

ci = zeros(size(R_feat_mat,1),2);
for i = 1:size(R_feat_mat,1)
    [~,~,ci(i,:),~] = ttest(R_feat_mat(i,:));
end

% Test difference of individual features
p_fix_cuts = signrank(R_feat_mat(ismember(feature_combos, 'Fixations'), :), ...
    R_feat_mat(ismember(feature_combos, 'Cuts'), :));
med_fix_cuts = median(R_feat_mat(ismember(feature_combos, 'Fixations'), :) ...
    - R_feat_mat(ismember(feature_combos, 'Cuts'), :));

p_fix_sound = signrank(R_feat_mat(ismember(feature_combos, 'Fixations'), :), ...
    R_feat_mat(ismember(feature_combos, 'Sound'), :));
med_fix_sound = median(R_feat_mat(ismember(feature_combos, 'Fixations'), :) ...
    - R_feat_mat(ismember(feature_combos, 'Sound'), :));

p_cuts_sound = signrank(R_feat_mat(ismember(feature_combos, 'Cuts'), :), ...
    R_feat_mat(ismember(feature_combos, 'Sound'), :));
med_cuts_sound = median(R_feat_mat(ismember(feature_combos, 'Cuts'), :) ...
    - R_feat_mat(ismember(feature_combos, 'Sound'), :));

p_all_sound = signrank(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :), ...
    R_feat_mat(ismember(feature_combos, 'Sound'), :));
med_all_sound = median(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :) ...
    - R_feat_mat(ismember(feature_combos, 'Sound'), :));

p_all_cuts = signrank(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :), ...
    R_feat_mat(ismember(feature_combos, 'Cuts'), :));
med_all_cuts = median(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :) ...
    - R_feat_mat(ismember(feature_combos, 'Cuts'), :));

p_all_fix = signrank(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :), ...
    R_feat_mat(ismember(feature_combos, 'Fixations'), :));
med_all_fix = median(R_feat_mat(ismember(feature_combos, 'Fixation, Cuts and Sound'), :) ...
    - R_feat_mat(ismember(feature_combos, 'Fixations'), :));

fprintf('Fixations vs cuts: median%cR=%1.1e, p=%1.5f, N=%d\n', ...
    916, med_fix_cuts, p_fix_cuts, size(R_feat_mat,2))
fprintf('Fixations vs sound: median%cR=%1.1e, p=%1.5f, N=%d\n', ...
    916, med_fix_sound, p_fix_sound, size(R_feat_mat,2))
fprintf('Cuts vs sound: median%cR=%1.1e, p=%1.5f, N=%d\n', ...
    916, med_cuts_sound, p_cuts_sound, size(R_feat_mat,2))
fprintf('All vs sound: median%cR=%1.1e, p=%1.5f, N=%d\n', ...
    916, med_all_sound, p_all_sound, size(R_feat_mat,2))
fprintf('All vs cuts: median%cR=%1.1e, p=%1.5f, N=%d\n', ...
    916, med_all_cuts, p_all_cuts, size(R_feat_mat,2))
fprintf('All vs fixations: median%cR=%1.1e, p=%1.5f, N=%d\n', ...
    916, med_all_fix, p_all_fix, size(R_feat_mat,2))

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

imagesc(mat_varx)

clim(ax1, plot_range)
axis square
colormap(ax1, slanCM(cmap))
xlabel('Channels')
ylabel('Channels')
title('Inputs')

ax2 = nexttile;

imagesc(mat_var)

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
title('No Inputs')

ax3 = nexttile;

imagesc(mat_varx - mat_var)

clim(ax3, 0.05*plot_range_diff)

axis square
colormap(ax3, slanCM('bwr'))
cb = colorbar(); 
ylabel(cb,['\Delta log p-value' newline '(Input - No Input)'] ,'Rotation',90)
xlabel('Channels')
yticks([])
title('Difference')

fontsize(fig1, fig_font, 'points')

exportgraphics(fig1, sprintf('%s/fig3_varx_full_none_example_%s_%s.png', ...
    fig_dir, example_pat, signal_type), 'Resolution', 600)

%% Difference of significant connections and R values over all paitents

if poster_size 
    figure('Units', 'inches', 'Position', [1,1,8.5,5]);
else
    figure('Position', [400,300,725,450])
end

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
ylabel(['\DeltaRatio' newline '(Input - No Input)'])

title(['Significant' newline 'Connections'])

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
ylabel(['Mean \DeltaR' newline '(Input - No Input)'])

title(['Effect' newline 'size'])

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_varx_full_none_sig_ch_diff_R_%s_sig_ch_%d.png', ...
    fig_dir, signal_type, select_sig), 'Resolution', 600)

%% Features
if poster_size 
    figure('Units', 'inches', 'Position', [1,1,7.5,5]);
    feature_combos = cellfun(@(C) strrep(C, 'and', '&'), feature_combos, 'UniformOutput', false);
    feature_combos = cellfun(@(C) strrep(C, 'Fixations', 'Fix.'), feature_combos, 'UniformOutput', false);
    feature_combos = cellfun(@(C) strrep(C, 'Fixation', 'Fix.'), feature_combos, 'UniformOutput', false);
else
    figure('Position', [725,700,700,525])
end

hold on

plot(flipud(mean(R_feat_mat,2)), 'r*-', 'LineWidth', 2)
plot(flipud(ci(:,1)), 'k*-')
plot(flipud(ci(:,2)), 'k*-')

xticklabels(fliplr(feature_combos))
xtickangle(45)

xlim([0.5, length(ci)+0.5])

ylabel('Mean \DeltaR')
title('Input - No Input')
grid('on')

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig3_R_features_%s_sig_ch_%d.png', ...
    fig_dir, signal_type, select_sig), 'Resolution', 600)
