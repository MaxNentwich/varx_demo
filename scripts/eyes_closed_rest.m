
%%
model_dir = '/media/max/Workspace/Code/varx_demo/results/models_revision_1';
fig_dir = '/media/max/Workspace/Code/varx_demo/results/figures';

pats_rest_closed = {'NS136', 'NS138', 'NS140', 'NS142', 'NS144', 'NS145', 'NS155_02', 'NS174_03', 'NS178'};

example_pat = 'NS174_03';

p_thresh = 0.0001;

%% 

mean_vec_rest = nan(1, length(pats_rest_closed));
mean_vec_rest_shift = nan(1, length(pats_rest_closed));
mean_vec_closed = nan(1, length(pats_rest_closed));
mean_vec_dme = nan(1, length(pats_rest_closed));

d_sig_rest = nan(1, length(pats_rest_closed));
d_sig_rest_shift = nan(1, length(pats_rest_closed));
d_sig_closed = nan(1, length(pats_rest_closed));
d_sig_dme_shift = nan(1, length(pats_rest_closed));

for pat = 1:length(pats_rest_closed)

    % Load data
    load(sprintf('%s/%s_varx_models.mat', model_dir, pats_rest_closed{pat}), 'm_varx_mov', 'vid_recs', 'labels')
    
    % Find the 5 minute version of Despicable Me and Resting state data
    idx_dme_5 = ismember(vid_recs, 'Despicable_Me_English_5min');
    idx_dme_shift = ismember(vid_recs, 'Despicable_Me_English_5min_shift');
    idx_rest = ismember(vid_recs, 'Resting_fixation');
    idx_rest_shift = ismember(vid_recs, 'Resting_fixation_shift');
    idx_rest_closed = ismember(vid_recs, 'Eyes_Closed_Rest');
    
    % Compute R values
    R_dme = abs(m_varx_mov{idx_dme_5}.A_Rvalue);  
    R_dme_shift = abs(m_varx_mov{idx_dme_shift}.A_Rvalue);  
    R_rest = abs(m_varx_mov{idx_rest}.A_Rvalue);  
    R_rest_shift = abs(m_varx_mov{idx_rest_shift}.A_Rvalue);
    R_rest_closed = abs(m_varx_mov{idx_rest_closed}.A_Rvalue); 
    
    % p-values
    p_dme = m_varx_mov{idx_dme_5}.A_pval;
    p_dme_shift = m_varx_mov{idx_dme_shift}.A_pval;
    p_rest = m_varx_mov{idx_rest}.A_pval;
    p_rest_shift = m_varx_mov{idx_rest_shift}.A_pval;
    p_closed = m_varx_mov{idx_rest_closed}.A_pval;
    
    idx_sig = p_closed < p_thresh & p_dme < p_thresh & p_rest < p_thresh & p_rest_shift < p_thresh & p_dme_shift < p_thresh;
    idx_sig = idx_sig - (eye(size(idx_sig)) .* diag(idx_sig));
    
    %% Plot matrices
    if strcmp(pats_rest_closed{pat}, 'NS155_02')
      
        example_mat = {R_dme, R_rest, R_rest_shift, R_rest_closed, R_dme_shift};
        example_titles = {'Despicable Me (VARX)', 'Resting Fixation (VARX)', 'Resting Fixation (VAR)', 'Eyes Closed (VAR)', 'Despicable Me (VAR)'}; 
        example_short = {'R_dme', 'R_rest', 'R_rest_shift', 'R_rest_closed', 'R_dme_shift'};

        for e = 1:length(example_mat)

            figure
            imagesc(example_mat{e})
            colormap(slanCM('Reds'))
            clim([0, 0.15])
            colorbar
            axis square
            xlabel('Channels')
            ylabel('Channels')
            title(example_titles{e})
            fontsize(16, 'Points')
            exportgraphics(gcf, ...
                sprintf('%s/revision_S2_%s_mat.png', fig_dir, example_short{e}), ...
                'Resolution', 300)
    
        end

    end
    
    %% Difference of significant connections
    R_dme(idx_sig ~= 1) = 0;
    R_rest(idx_sig ~= 1) = 0;
    R_dme_shift(idx_sig ~= 1) = 0;
    R_rest_shift(idx_sig ~= 1) = 0;
    R_rest_closed(idx_sig ~= 1) = 0;
    
    % Compute difference of R values
    R_diff_rest = R_dme - R_rest;
    R_diff_rest_shift = R_dme - R_rest_shift;
    R_diff_closed = R_dme - R_rest_closed;
    R_diff_dme = R_dme - R_dme_shift;
    
    % Vectorize and compute median
    Rd_vec_rest = R_diff_rest(idx_sig == 1);
    mean_vec_rest(pat) = mean(Rd_vec_rest);
    
    Rd_vec_rest_shift = R_diff_rest_shift(idx_sig == 1);
    mean_vec_rest_shift(pat) = mean(Rd_vec_rest_shift);
    
    Rd_vec_closed = R_diff_closed(idx_sig == 1);
    mean_vec_closed(pat) = mean(Rd_vec_closed);
    
    Rd_vec_dme = R_diff_dme(idx_sig == 1);
    mean_vec_dme(pat) = mean(Rd_vec_dme);
    
    %% Significant connetions
    n_ch = size(p_dme,1);
    n_conn = n_ch^2 - n_ch;
    
    n_sig_dme5 = (sum(p_dme(:) < p_thresh) - n_ch) / n_conn;
    n_sig_dme_shift = (sum(p_dme_shift(:) < p_thresh) - n_ch) / n_conn;
    n_sig_rest = (sum(p_rest(:) < p_thresh) - n_ch) / n_conn;
    n_sig_rest_shift = (sum(p_rest_shift(:) < p_thresh) - n_ch) / n_conn;
    n_sig_closed = (sum(p_closed(:) < p_thresh) - n_ch) / n_conn;
    
    d_sig_rest(pat) = n_sig_dme5 - n_sig_rest;
    d_sig_rest_shift(pat) = n_sig_dme5 - n_sig_rest_shift;
    d_sig_closed(pat) = n_sig_dme5 - n_sig_closed;
    d_sig_dme_shift(pat) = n_sig_dme5 - n_sig_dme_shift;

end

%% Figures
figure
hold on
bar([mean(mean_vec_rest), mean(mean_vec_rest_shift), mean(mean_vec_closed), mean(mean_vec_dme)], ...
    'FaceColor',[0.5 .5 .5],'EdgeColor',[0 0 0]);
plot([mean_vec_rest; mean_vec_rest_shift; mean_vec_closed; mean_vec_dme], '.-', 'MarkerSize', 10)
xticks(1:4)
xticklabels({'DME - Rest', 'DME - Rest (VAR)', 'DME - Eyes Closed', 'DME - DME (VAR)'})
ylabel('\DeltaR')
title('R-values')
fontsize(16, 'Points')
exportgraphics(gcf, sprintf('%s/revision_S2_R_vals_bar.png', fig_dir), 'Resolution', 300)

figure
hold on
bar([mean(d_sig_rest), mean(d_sig_rest_shift), mean(d_sig_closed), mean(d_sig_dme_shift)], ...
    'FaceColor',[0.5 .5 .5],'EdgeColor',[0 0 0]);
plot([d_sig_rest; d_sig_rest_shift; d_sig_closed; d_sig_dme_shift], '.-', 'MarkerSize', 10)
xticks(1:4)
xticklabels({'DME - Rest', 'DME - Rest (VAR)', 'DME - Eyes Closed', 'DME - DME (VAR)'})
ylabel('\Delta Sig. Chns.')
title('Significant Channels')
fontsize(16, 'Points')
exportgraphics(gcf, sprintf('%s/revision_S2_p_vals_bar.png', fig_dir), 'Resolution', 300)

%% Stats
fprintf('Resting Fixation VARX, Delta R = %1.4f, p = %1.2f\n', mean(mean_vec_rest), signrank(mean_vec_rest))
fprintf('Resting Fixation VAR, Delta R = %1.4f, p = %1.2f\n', mean(mean_vec_rest_shift), signrank(mean_vec_rest_shift))
fprintf('Eyes Closed Rest VARX, Delta R = %1.4f, p = %1.2f\n', mean(mean_vec_closed), signrank(mean_vec_closed))
fprintf('DME, Delta R = %1.4f, p = %1.2f\n', mean(mean_vec_dme), signrank(mean_vec_dme))

fprintf('Resting Fixation VARX, Sig. Chans = %1.4f, p = %1.2f\n', mean(d_sig_rest), signrank(d_sig_rest))
fprintf('Resting Fixation VAR, Sig. Chans = %1.4f, p = %1.2f\n', mean(d_sig_rest_shift), signrank(d_sig_rest_shift))
fprintf('Eyes Closed Rest VARX, Sig. Chans = %1.4f, p = %1.2f\n', mean(d_sig_closed), signrank(d_sig_closed))
fprintf('DME, Sig. Chans = %1.4f, p = %1.2f\n', mean(d_sig_dme_shift), signrank(d_sig_dme_shift))