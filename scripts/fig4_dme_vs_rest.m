%% Figure 1 Difference between Despicable Me English and resting fixation 

%% Figure settings
fig_font = 16;

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

coord_dir = '../results/spatial_plots';
if exist(coord_dir, 'dir') == 0, mkdir(coord_dir), end

addpath('../src/Violinplot-Matlab-master')
addpath('../src')

%% Define patients
example_pat = 'NS127_02';

signal_type = 'HFA';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient(table2array(sum(patient_list(:,2:end),2) == 2));

% Movie table for coordinates
load('../data/movie_subs_table.mat', 'movie_subs_table');

% Matrix to convert fsaverage to MNI152
% (https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems)
fsaverage2mni152 = [0.9975, -0.0073, 0.0176, -0.0429; ...
    0.0146, 1.0009, -0.0024, 1.5496; ...
    -0.0130, -0.0093, 0.9971, 1.1840];

% Select significant channels or not
sig_channels = false;
p_thresh = 0.001;

% Initialize array to collect data
varX_Rvalue_dme = cell(1, length(patients));
varX_Rvalue_rest = cell(1, length(patients));
varx_R_diff = cell(1, length(patients));
varx_diff_mean = nan(1, length(patients));
varx_diff_dme_mean = nan(1, length(patients));
n_sig_diff = nan(1, length(patients));
coords = cell(1, length(patients));

for pat = 1:length(patients)
    
    % Load data
    if strcmp(signal_type, 'LFP')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx_mov', 'vid_recs', 'labels')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx_mov_hfa', 'vid_recs', 'labels')
        m_varx_mov = m_varx_mov_hfa;
    end

    % Find the 5 minute version of Despicable Me and Resting state data
    idx_dme_5 = ismember(vid_recs, 'Despicable_Me_English_5min');
    idx_dme_shift = ismember(vid_recs, 'Despicable_Me_English_5min_shift');
    idx_rest = ismember(vid_recs, 'Resting_fixation');

    if sum(idx_dme_5) == 0, continue, end

    % Compute R values
    varX_Rvalue_dme{pat} = abs(sqrt(1-exp(-m_varx_mov{idx_dme_5}.A_Deviance/m_varx_mov{idx_dme_5}.T)));  
    varX_Rvalue_dme_shift = abs(sqrt(1-exp(-m_varx_mov{idx_dme_shift}.A_Deviance/m_varx_mov{idx_dme_shift}.T)));  
    varX_Rvalue_rest{pat} = abs(sqrt(1-exp(-m_varx_mov{idx_rest}.A_Deviance/m_varx_mov{idx_rest}.T)));  

    %% Difference of significant connections
    if sig_channels

        idx_sig = m_varx_mov{idx_dme_5}.A_pval < p_thresh & m_varx_mov{idx_rest}.A_pval < p_thresh;
        idx_sig = idx_sig - (eye(size(idx_sig)) .* diag(idx_sig));

        idx_sig_shift = m_varx_mov{idx_dme_5}.A_pval < p_thresh & m_varx_mov{idx_dme_shift}.A_pval < p_thresh;
        idx_sig_shift = idx_sig_shift - (eye(size(idx_sig_shift)) .* diag(idx_sig_shift));

    else
        idx_sig = eye(size(varX_Rvalue_dme{pat})) ~= 1;
        idx_sig_shift = idx_sig;
    end

    varX_Rvalue_dme{pat}(idx_sig ~= 1) = 0;
    varX_Rvalue_rest{pat}(idx_sig ~= 1) = 0;
    varX_Rvalue_dme_shift(idx_sig_shift ~= 1) = 0;

    % Compute difference of R values
    varx_R_diff{pat} = varX_Rvalue_dme{pat} - varX_Rvalue_rest{pat};
    varx_R_diff_dme = varX_Rvalue_dme{pat} - varX_Rvalue_dme_shift;

    % Vectorize and compute median
    varx_R_diff_vec = varx_R_diff{pat}(idx_sig == 1);
    varx_diff_mean(pat) = mean(varx_R_diff_vec);

    varx_R_diff_dme_vec = varx_R_diff_dme(idx_sig_shift == 1);
    varx_diff_dme_mean(pat) = mean(varx_R_diff_dme_vec);

    %% Significant connetions
    n_ch = size(m_varx_mov{idx_dme_5}.A_pval,1);
    n_conn = n_ch^2 - n_ch;

    n_sig_dme5 = (sum(m_varx_mov{idx_dme_5}.A_pval(:) < p_thresh) - n_ch) / n_conn;
    n_sig_rest = (sum(m_varx_mov{idx_rest}.A_pval(:) < p_thresh) - n_ch) / n_conn;

    n_sig_diff(pat) = n_sig_dme5 - n_sig_rest;

    %% Get coordinates for spatial plots
    if strcmp(patients{pat}, 'NS174_02')
        pat_coord = 'NS174';
    else
        pat_coord = patients{pat};
    end

    movie_table_pat = movie_subs_table(ismember(movie_subs_table.SubID, pat_coord), :);

    labels_coord =  cellfun(@(C) sprintf('%s_%s', pat_coord, C), labels, 'UniformOutput', false);
    [~, coords{pat}, ~] = localize_elecs_bipolar(labels_coord, 'dk');

    % Convert to MNI152
    coords{pat} = mean(coords{pat}, 3);
    coords{pat} = fsaverage2mni152 * [coords{pat}, zeros(size(coords{pat},1), 1)]';
    coords{pat} = coords{pat}';

end

%% Remove empty patients
idx_empty = isnan(varx_diff_mean);
varx_diff_mean(idx_empty) = [];
varx_diff_dme_mean(idx_empty) = [];

n_sig_diff(idx_empty) = [];

%% Compute differences and stats
[~, p_n_sig, ~, stat_n_sig] = ttest(n_sig_diff);
[~, p_effect_size, ~, stat_effect_size] = ttest(varx_diff_mean);
[~, p_movie, ~, stat_movie] = ttest(varx_diff_dme_mean);

%% Save the data for a spatial plot
idx_example = ismember(patients, example_pat);

varX_Rvalue_rest = varX_Rvalue_rest{idx_example};
varX_Rvalue_dme = varX_Rvalue_dme{idx_example};

coords = coords{idx_example};
save(sprintf('%s/fig4_matrices_rest_movie_coords.mat', coord_dir), ...
    'varX_Rvalue_rest', 'varX_Rvalue_dme', 'coords')

%% Plot an example
plot_max = max([max(varX_Rvalue_dme(:)), max(varX_Rvalue_rest(:))]);
plot_range = [0, plot_max];

figure('Position', [3000,350,1450,450]);

tiledlayout(1,3);

% Connectitivy of Despicable Me
ax1 = nexttile;

imagesc(varX_Rvalue_dme)

clim(0.2*plot_range)
axis square
colormap(ax1, slanCM('Reds'))
xlabel('Channels')
ylabel('Channels')
title('Movie')

% Connectitivy of Resting State
ax2 = nexttile;

imagesc(varX_Rvalue_rest)

clim(0.2*plot_range)
axis square
colormap(ax2, slanCM('Reds'))
cb = colorbar(); 
ylabel(cb,'R' ,'Rotation',90)
yticks([])
xlabel('Channels')
title('Rest')

% Difference
ax3 = nexttile;

imagesc(varx_R_diff{idx_example})

clim([-0.1*plot_range(2), 0.1*plot_range(2)])

axis square
colormap(ax3, slanCM('bwr'))
cb = colorbar(); 
ylabel(cb,'\DeltaR (Movie - Rest)' ,'Rotation',90)
yticks([])
xlabel('Channels')
title('Difference')

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_example_%s_%s.png', ...
    fig_dir, example_pat, signal_type), 'Resolution', 300)

%% Difference of all connection for one patient and summary for all patients
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
ylabel('\Delta Ratio (Movie - Rest)')

ax = ancestor(gca, 'axes');
ax.YAxis.Exponent = 0;
if strcmp(signal_type, 'LFP')
    ytickformat('%0.2f')
elseif strcmp(signal_type, 'HFA')
    ytickformat('%0.3f')
end

title('Sig. Connections')

% Effect size
nexttile(6,[1,4])
hold on 

scatter(0.1*randn(1, length(varx_diff_mean)), varx_diff_mean, 'k', 'filled')

plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
plot([-0.25, 0.25], median(varx_diff_mean)*ones(2,1), 'k--')

grid on
box on
xticks([])
xlim([-0.3, 0.3])
ylim_abs =  1.2 * max(abs(varx_diff_mean));
ylim([-ylim_abs, ylim_abs])
set(gca, 'YAxisLocation', 'right')
ylabel('\DeltaR (Movie - Rest)')

ax = ancestor(gca, 'axes');
ax.YAxis.Exponent = 0;
ytickformat('%0.4f')

title('Effect Size')

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_summary_%s.png', ...
    fig_dir, signal_type), 'Resolution', 300)

%% Positive control showing movie data does differ when more features are included
figure('Position', [400,300,350,450])
hold on 

scatter(0.1*randn(1, length(varx_diff_dme_mean)), varx_diff_dme_mean, 'k', 'filled')

plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
plot([-0.25, 0.25], median(varx_diff_dme_mean)*ones(2,1), 'k--')

grid on
box on
xticks([])
xlim([-0.3, 0.3])
ylim_abs =  1.2 * max(abs(varx_diff_dme_mean));
ylim([-ylim_abs, ylim_abs])
set(gca, 'YAxisLocation', 'right')
ylabel('\DeltaR (Inputs - No Inputs)')

ax = ancestor(gca, 'axes');
ax.YAxis.Exponent = 0;
ytickformat('%0.5f')

title('Movie')

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_features_summary_%s.png', ...
    fig_dir, signal_type), 'Resolution', 300)

%% Run python file to plot connections
pyenv('Version','/home/max/anaconda3/envs/varx_plots/bin/python3.8')
pyrunfile('connectivity_plot_movie_vs_rest.py')
