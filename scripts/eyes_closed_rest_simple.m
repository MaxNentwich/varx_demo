
%%
model_dir = '../results/models_revision_1';

fig_dir = '/media/max/Workspace/Code/varx_demo/results/figures';

pats_rest_closed = {'NS136', 'NS138', 'NS140', 'NS142', 'NS144', 'NS155_02', 'NS174_03', 'NS178'};

p_thresh = 0.0001;

fig_font = 20;

%% 
mean_vec_closed = nan(1, length(pats_rest_closed));
d_sig_closed = nan(1, length(pats_rest_closed));

R_movie = nan(1, length(pats_rest_closed));
R_rest = nan(1, length(pats_rest_closed));

n_sig_movie = nan(1, length(pats_rest_closed));
n_sig_rest = nan(1, length(pats_rest_closed));

for pat = 1:length(pats_rest_closed)

    % Load data
    load(sprintf('%s/%s_varx_models.mat', model_dir, pats_rest_closed{pat}), 'm_varx_mov', 'vid_recs', 'labels')


    % Find the 5 minute version of Despicable Me and Resting state data
    idx_dme_5 = ismember(vid_recs, 'Despicable_Me_English_5min');
    idx_rest_closed = ismember(vid_recs, 'Eyes_Closed_Rest');
    
    % Compute R values
    R_dme = abs(m_varx_mov{idx_dme_5}.A_Rvalue);  
    R_rest_closed = abs(m_varx_mov{idx_rest_closed}.A_Rvalue); 
    
    % p-values
    p_dme = m_varx_mov{idx_dme_5}.A_pval;
    p_closed = m_varx_mov{idx_rest_closed}.A_pval;
    
    idx_sig = p_closed < p_thresh & p_dme < p_thresh;
    idx_sig = idx_sig - (eye(size(idx_sig)) .* diag(idx_sig));
    
    %% Plot matrices
    if strcmp(pats_rest_closed{pat}, 'NS155_02')
      
        example_mat = {R_dme, R_rest_closed};
        example_titles = {'Movie', 'Rest (eyes-closed)'}; 
        example_short = {'R_dme', 'R_rest_closed'};

        for e = 1:length(example_mat)

            figure
            imagesc(example_mat{e})
            colormap(slanCM('Reds'))
            clim([0, 0.15])
            cb = colorbar;
            cb.Label.String = 'R';
            axis square
            xlabel('Channels')
            ylabel('Channels')
            title(example_titles{e})
            fontsize(fig_font, 'Points')
            exportgraphics(gcf, ...
                sprintf('%s/revision_S2_%s_mat.png', fig_dir, example_short{e}), ...
                'Resolution', 300)
    
        end

    end
    
    %% Difference of significant connections
    R_dme(idx_sig ~= 1) = 0;
    R_rest_closed(idx_sig ~= 1) = 0;
    
    % Compute difference of R values
    R_diff_closed = R_dme - R_rest_closed;

    R_movie(pat) = mean(R_dme(:));
    R_rest(pat) = mean(R_rest_closed(:));

    Rd_vec_closed = R_diff_closed(idx_sig == 1);
    mean_vec_closed(pat) = mean(Rd_vec_closed);
    
    %% Significant connetions
    n_ch = size(p_dme,1);
    n_conn = n_ch^2 - n_ch;
    
    n_sig_movie(pat) = (sum(p_dme(:) < p_thresh) - n_ch) / n_conn;
    n_sig_rest(pat) = (sum(p_closed(:) < p_thresh) - n_ch) / n_conn;
    
    d_sig_closed(pat) = n_sig_movie(pat) - n_sig_rest(pat);

end

%% Stats
fprintf('Eyes Closed Rest VARX, Delta R = %1.4f, p = %1.2f, N = %d\n', median(mean_vec_closed), signrank(mean_vec_closed), length(mean_vec_closed))
fprintf('Eyes Closed Rest VARX, Sig. Chans = %1.4f, p = %1.2f, N = %d\n', median(d_sig_closed), signrank(d_sig_closed), length(d_sig_closed))

%% Figures
R_plot = [R_rest; R_movie];
ch_plot = [n_sig_rest; n_sig_movie];

figure('Position', [1000,750,750,450])

% Significant Channels
tiledlayout(1,10);

nexttile(1,[1,4])
hold on

plot(ch_plot(:, diff(ch_plot) > 0), '.-', 'Color', [0.75, 0.3, 0], 'LineWidth', 2, 'MarkerSize', 20)
plot(ch_plot(:, diff(ch_plot) < 0), '.-', 'Color', [0, 0.3, 0.75], 'LineWidth', 2, 'MarkerSize', 20)

xticks([1, 2])
xticklabels({'Rest', 'Movie'})
xtickangle(45)
xlim([0.75, 2.25])

ylabel('Fraction of Channels')
set(gca, 'YAxisLocation', 'right')
ylim([0.9*min(ch_plot(:)), 1.1*max(ch_plot(:))])

legend(['Increase', repmat({''}, 1, length(R_plot)-2), 'Decrease'], ...
    'Position', [0.77,0.02,0.22,0.15]);

fontsize(fig_font, 'Points')
title(['Significant' newline 'Connections'])

grid on

% Effect size
nexttile(5,[1,4])
hold on 

plot(R_plot(:, diff(R_plot) > 0), '.-', 'Color', [0.75, 0.3, 0], 'LineWidth', 2, 'MarkerSize', 20)
plot(R_plot(:, diff(R_plot) < 0), '.-', 'Color', [0, 0.3, 0.75], 'LineWidth', 2, 'MarkerSize', 20)

xticks([1, 2])
xticklabels({'Movies', 'Rest'})
xtickangle(45)
xlim([0.75, 2.25])

ylabel('R')
set(gca, 'YAxisLocation', 'right')
ylim([0.95*min(R_plot(:)), 1.05*max(R_plot(:))])

fontsize(fig_font, 'Points')
title('Effect Size')

grid on

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_summary_LFP_eyes_closed_rest.png', ...
    fig_dir), 'Resolution', 600)

% tiledlayout(1,9);
% 
% % Violin plot of differences for one paitent
% nexttile(1,[1,4])
% hold on
% 
% scatter(0.1*randn(1, length(d_sig_closed)), d_sig_closed, 'k', 'filled')
% 
% plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
% plot([-0.25, 0.25], median(d_sig_closed)*ones(2,1), 'k--')
% 
% grid on
% box on
% xticks([])
% xlim([-0.3, 0.3])
% ylim_abs =  1.2 * max(abs(d_sig_closed));
% ylim([-ylim_abs, ylim_abs])
% set(gca, 'YAxisLocation', 'right')
% ylabel(['\DeltaRatio' newline '(Movie - Rest)'])
% 
% ax = ancestor(gca, 'axes');
% ax.YAxis.Exponent = 0;
% ytickformat('%0.3f')
% 
% title(['Significant' newline 'Connections'])
% 
% % Effect size
% nexttile(6,[1,4])
% hold on 
% 
% scatter(0.1*randn(1, length(mean_vec_closed)), mean_vec_closed, 'k', 'filled')
% 
% plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
% plot([-0.25, 0.25], median(mean_vec_closed)*ones(2,1), 'k--')
% 
% grid on
% box on
% xticks([])
% xlim([-0.3, 0.3])
% ylim_abs =  1.2 * max(abs(mean_vec_closed));
% ylim([-ylim_abs, ylim_abs])
% set(gca, 'YAxisLocation', 'right')
% ylabel(['\DeltaR' newline '(Movie - Rest)'])
% 
% title(['Effect' newline 'Size'])
% 
% fontsize(gcf, fig_font, 'points')
% 
% exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_summary_LFP_eyes_closed_rest.png', ...
%     fig_dir), 'Resolution', 600)