%% Figure 1 Difference between Despicable Me English and resting fixation 

%% Figure settings
fig_font = 20;

model_dir = '../results/models_revision_1';

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

coord_dir = '../results/spatial_plots';
if exist(coord_dir, 'dir') == 0, mkdir(coord_dir), end

addpath('../src/Violinplot-Matlab-master')
addpath('../src')

%% Define patients
example_pat = 'NS137';

signal_type = 'LFP';

% Movie segment ['Despicable_Me_English_5min', 'Despicable_Me_English_last_5min', 'Inscapes_5min', 'Inscapes_last_5min', 'Monkey_5min', 'Despicable_Me_English_5min_shift']
movie_select = 'Despicable_Me_English_5min';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient(table2array(sum(patient_list(:,2:3),2) == 2));

if strcmp(movie_select, 'Eyes_Closed_Rest')
    example_pat = 'NS136';
end

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

% Poster settings
poster_size = true;

if poster_size
    fig_font = 26;
end

% Initialize array to collect data
varX_Rvalue_dme = cell(1, length(patients));
varX_Rvalue_rest = cell(1, length(patients));
varX_pvalue_dme = cell(1, length(patients));
varX_Rvalue_dme_shift = cell(1, length(patients));
varX_pvalue_rest = cell(1, length(patients));
varx_R_diff = cell(1, length(patients));
varx_diff_mean = nan(1, length(patients));
varx_diff_dme_mean = nan(1, length(patients));
n_sig_diff = nan(1, length(patients));
coords = cell(1, length(patients));

for pat = 1:length(patients)
    
    % Load data
    if strcmp(signal_type, 'LFP')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx_mov', 'vid_recs', 'labels')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx_mov_hfa', 'vid_recs', 'labels')
        m_varx_mov = m_varx_mov_hfa;
    end

    % Find the 5 minute version of Despicable Me and Resting state data
    idx_dme_5 = ismember(vid_recs, movie_select);
    idx_dme_shift = ismember(vid_recs, 'Despicable_Me_English_5min_shift');
    idx_rest = ismember(vid_recs, 'Resting_fixation');

    if sum(idx_dme_5) == 0, continue, end

    % Compute R values
    varX_Rvalue_dme{pat} = abs(sqrt(1-exp(-m_varx_mov{idx_dme_5}.A_Deviance/m_varx_mov{idx_dme_5}.T)));  
    varX_Rvalue_dme_shift{pat} = abs(sqrt(1-exp(-m_varx_mov{idx_dme_shift}.A_Deviance/m_varx_mov{idx_dme_shift}.T)));  
    varX_Rvalue_rest{pat} = abs(sqrt(1-exp(-m_varx_mov{idx_rest}.A_Deviance/m_varx_mov{idx_rest}.T)));  

    %% Difference of significant connections
    varX_pvalue_dme{pat} = m_varx_mov{idx_dme_5}.A_pval;
    varX_pvalue_rest{pat} = m_varx_mov{idx_rest}.A_pval;

    if sig_channels

        idx_sig = varX_pvalue_dme{pat} < p_thresh & varX_pvalue_rest{pat} < p_thresh;
        idx_sig = idx_sig - (eye(size(idx_sig)) .* diag(idx_sig));

        idx_sig_shift = varX_pvalue_dme{pat} < p_thresh & m_varx_mov{idx_dme_shift}.A_pval < p_thresh;
        idx_sig_shift = idx_sig_shift - (eye(size(idx_sig_shift)) .* diag(idx_sig_shift));

    else
        idx_sig = eye(size(varX_Rvalue_dme{pat})) ~= 1;
        idx_sig_shift = idx_sig;
    end

    varX_Rvalue_dme{pat}(idx_sig ~= 1) = 0;
    varX_Rvalue_rest{pat}(idx_sig ~= 1) = 0;
    varX_Rvalue_dme_shift{pat}(idx_sig_shift ~= 1) = 0;

    % Compute difference of R values
    varx_R_diff{pat} = varX_Rvalue_dme{pat} - varX_Rvalue_rest{pat};
    varx_R_diff_dme = varX_Rvalue_dme{pat} - varX_Rvalue_dme_shift{pat};

    % Vectorize and compute median
    varx_R_diff_vec = varx_R_diff{pat}(idx_sig == 1);
    varx_diff_mean(pat) = mean(varx_R_diff_vec);

    varx_R_diff_dme_vec = varx_R_diff_dme(idx_sig_shift == 1);
    varx_diff_dme_mean(pat) = mean(varx_R_diff_dme_vec);

    %% Significant connetions
    n_ch = size(varX_pvalue_dme{pat},1);
    n_conn = n_ch^2 - n_ch;

    n_sig_dme5 = (sum(varX_pvalue_dme{pat}(:) < p_thresh) - n_ch) / n_conn;
    n_sig_rest = (sum(varX_pvalue_rest{pat}(:) < p_thresh) - n_ch) / n_conn;

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

varX_Rvalue_dme(idx_empty) = [];
varX_Rvalue_dme_shift(idx_empty) = [];

n_sig_diff(idx_empty) = [];

%% Compute differences and stats
p_n_sig = signrank(n_sig_diff);
p_effect_size = signrank(varx_diff_mean);
p_movie = signrank(varx_diff_dme_mean);

fprintf('Ratio of significant channels: median=%1.1e, p=%1.2f, N=%d\n', ...
    median(n_sig_diff), p_n_sig, length(n_sig_diff))
fprintf('Effect size: median=%1.1e, p=%1.2f, N=%d\n', ...
    median(varx_diff_mean), p_effect_size, length(varx_diff_mean))
fprintf('Movie features median=%1.1e, p=%1.5f, N=%d\n', ...
    median(varx_diff_dme_mean), p_movie, length(varx_diff_dme_mean))

%% Save the data for a spatial plot
R_dme = cellfun(@(C) mean(C(:)), varX_Rvalue_dme);
R_dme_shift = cellfun(@(C) mean(C(:)), varX_Rvalue_dme_shift);

idx_example = ismember(patients, example_pat);

varX_Rvalue_rest = varX_Rvalue_rest{idx_example};
varX_Rvalue_dme = varX_Rvalue_dme{idx_example};

varX_pvalue_rest = varX_pvalue_rest{idx_example};
varX_pvalue_dme = varX_pvalue_dme{idx_example};

coords = coords{idx_example};

plot_max = max([max(varX_Rvalue_dme(:)), max(varX_Rvalue_rest(:))]);
plot_range = [0, plot_max];

save(sprintf('%s/fig4_matrices_rest_movie_coords.mat', coord_dir), ...
    'varX_Rvalue_rest', 'varX_Rvalue_dme', 'varX_pvalue_rest', 'varX_pvalue_dme', ...
    'coords', 'plot_range')

%% Plot an example

if poster_size 
    figure('Units', 'inches', 'Position', [1,1,20,5]);
else
    figure('Position', [3000,350,1450,450]);
end

tiledlayout(1,3);

% Connectitivy of Despicable Me
ax1 = nexttile;

imagesc(varX_Rvalue_rest)

clim(plot_range)
axis square
colormap(ax1, flipud(slanCM('amber')))
colormap(ax1, slanCM('Reds'))
ylabel('Channels')
xlabel('Channels')
title('Rest')


% Connectitivy of Resting State
ax2 = nexttile;

imagesc(varX_Rvalue_dme)

clim(plot_range)
axis square
colormap(ax2, flipud(slanCM('amber')))
colormap(ax2, slanCM('Reds'))
cb = colorbar(); 
ylabel(cb,'R' ,'Rotation',90)
yticks([])
xlabel('Channels')
title('Movie')

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

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_example_%s_%s_%s.png', ...
    fig_dir, example_pat, signal_type, movie_select), 'Resolution', 600)

%% Difference of all connection for one patient and summary for all patients
if poster_size 
    figure('Units', 'inches', 'Position', [1,1,8.5,5]);
else
    figure('Position', [400,300,675,450])
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
ylabel(['\DeltaRatio' newline '(Movie - Rest)'])

ax = ancestor(gca, 'axes');
ax.YAxis.Exponent = 0;
ytickformat('%0.3f')

title(['Significant' newline 'Connections'])

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
ylabel(['\DeltaR' newline '(Movie - Rest)'])

title(['Effect' newline 'Size'])

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_summary_%s_%s.png', ...
    fig_dir, signal_type, movie_select), 'Resolution', 600)

%% Positive control showing movie data does differ when more features are included
R_plot = [R_dme_shift; R_dme];
R_plot = R_plot - R_plot(1,:);

figure('Position', [450,200,300,450])
hold on 

plot(R_plot(:, diff(R_plot) > 0), '.-', 'Color', [0.75, 0.3, 0], 'LineWidth', 2, 'MarkerSize', 20)
plot(R_plot(:, diff(R_plot) < 0), '.-', 'Color', [0, 0.3, 0.75], 'LineWidth', 2, 'MarkerSize', 20)

xticks([1, 2])
xticklabels({'No Input', 'Input'})
xtickangle(45)
xlim([0.75, 2.25])

ylabel('R')
set(gca, 'YAxisLocation', 'right')
ylim([1.05*min(R_plot(:)), 1.05*max(R_plot(:))])

fontsize(fig_font, 'Points')
title('Effect Size')

grid on

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_features_summary_%s.png', ...
    fig_dir, signal_type), 'Resolution', 300)

% figure('Position', [450,200,275,300])
% hold on 
% 
% scatter(0.1*randn(1, length(varx_diff_dme_mean)), varx_diff_dme_mean, 'k', 'filled')
% 
% plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
% plot([-0.25, 0.25], median(varx_diff_dme_mean)*ones(2,1), 'k--')
% 
% grid on
% box on
% xticks([])
% xlim([-0.3, 0.3])
% ylim_abs =  1.2 * max(abs(varx_diff_dme_mean));
% ylim([-ylim_abs, ylim_abs])
% set(gca, 'YAxisLocation', 'right')
% ylabel(['\DeltaR' newline '(Inputs - No Inputs)'])
% 
% title('Movie')
% 
% fontsize(gcf, fig_font, 'points')


%% Run python file to plot connections

% pyenv('Version','/home/max/anaconda3/envs/varx_plots/bin/python3.8')
% pyrunfile('connectivity_plot_movie_vs_rest.py')
