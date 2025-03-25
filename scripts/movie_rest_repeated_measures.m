%% Figure 1 Difference between Despicable Me English and resting fixation 

%% Figure settings
fig_font = 20;

model_dir = '../results/models_revision_1';

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

%% Define patients
signal_type = 'LFP';

% Movie segment 
conditions = {'Resting_fixation', 'Despicable_Me_English_5min', 'Despicable_Me_English_last_5min', 'Inscapes_5min', 'Monkey_5min'};

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient(table2array(sum(patient_list(:,2:3),2) == 2));

% Movie table for coordinates
load('../data/movie_subs_table.mat', 'movie_subs_table');

% Select significant channels or not
sig_channels = false;
p_thresh = 0.001;

if exist('movie_rest_stats.mat', 'file') == 0

    % Initialize array to collect data
    R_values = nan(length(patients), length(conditions));
    sig_chans = nan(length(patients), length(conditions));
    
    for pat = 1:length(patients)
        
        % Load data
        if strcmp(signal_type, 'LFP')
            load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx_mov', 'vid_recs', 'labels')
        elseif strcmp(signal_type, 'HFA')
            load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx_mov_hfa', 'vid_recs', 'labels')
            m_varx_mov = m_varx_mov_hfa;
        end
    
        for c = 1:length(conditions)
    
            idx_c = ismember(vid_recs, conditions{c});
    
            if sum(idx_c) == 0, continue, end
    
            R_val = abs(m_varx_mov{idx_c}.A_Rvalue);
            p_val = m_varx_mov{idx_c}.A_pval;
    
            idx_sig = eye(size(R_val)) ~= 1;
            R_values(pat, c) = mean(R_val(idx_sig == 1));
    
            % Significant connections
            n_ch = size(R_val,1);
            n_conn = n_ch^2 - n_ch;
        
            sig_chans(pat, c) = (sum(p_val(:) < p_thresh) - n_ch) / n_conn;
       
        end
    
    end
    
    %% Remove empty patients
    idx_empty = isnan(sum(R_values, 2));
    
    R_values(idx_empty, :) = [];
    sig_chans(idx_empty, :) = [];
    patients(idx_empty) = [];
    
    %% Compute differences and stats
    R_mat = [];
    chn_mat = [];
    
    for m = 2:size(R_values, 2)
        R_mat = [R_mat; [R_values(:,1), R_values(:,m)]];
        chn_mat = [chn_mat; [sig_chans(:,1), sig_chans(:,m)]];
    end
    
    conditions_vec = repmat(conditions(2:end), length(patients), 1);
    conditions_vec = conditions_vec(:);
    
    T_R = table(conditions_vec, R_mat(:,1), R_mat(:,2), ...
        'VariableNames', {'movie_type', 'rest', 'movies'});
    T_chn = table(conditions_vec, chn_mat(:,1), chn_mat(:,2), ...
        'VariableNames', {'movie_type', 'rest', 'movies'});
    
    save('movie_rest_stats.mat', 'T_R', 'T_chn', 'patients')

else
    load('movie_rest_stats.mat', 'T_R', 'T_chn', 'patients')
end

Cond = table([1 2]','VariableNames',{'Conditions'});



%% Tests
% 
% % Repeated measures ANOVA (parametric)
% rm_R = fitrm(T_R, 'rest-movies~movie_type', 'WithinDesign', Cond);
% ranovatbl_R = ranova(rm_R);
% 
% p_R_param = ranovatbl_R.pValue(1);
% 
% % Friedman test (non-parametric)
% p_R_non_param = friedman(table2array(T_R(:, 2:end), length(patients)));
% 
% %% Significant channels
% 
% % Repeated measures ANOVA (parametric)
% rm_chn = fitrm(T_chn, 'rest-movies~movie_type', 'WithinDesign', Cond);
% ranovatbl_chn = ranova(rm_chn);
% 
% p_chn_param = ranovatbl_chn.pValue(1);
% 
% % Friedman test (non-parametric)
% p_chn_non_param = friedman(table2array(T_chn(:, 2:end), length(patients)));

%% Reorganize again for a plot
movies = unique(T_R.movie_type);

R_mat = nan(length(movies), size(T_R,1) / length(movies));
ch_mat = nan(length(movies), size(T_chn,1) / length(movies));

for m = 1:length(movies)

    % R values
    R_mat(m,:) = T_R(ismember(T_R.movie_type, movies{m}), :).movies;
    R_rest = T_R(ismember(T_R.movie_type, movies{m}), :).rest';

    % Significant channels
    ch_mat(m,:) = T_chn(ismember(T_chn.movie_type, movies{m}), :).movies;
    ch_rest = T_chn(ismember(T_chn.movie_type, movies{m}), :).rest';

end

R_mat = mean(R_mat);
ch_mat = mean(ch_mat);

R_diff = R_mat - R_rest;
ch_diff = ch_mat - ch_rest;

%% Figure
figure('Position', [400,300,675,450])

tiledlayout(1,9);

nexttile(1,[1,4])
hold on

scatter(0.1*randn(1, length(ch_diff)), ch_diff, 'k', 'filled')

plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
plot([-0.25, 0.25], median(ch_diff)*ones(2,1), 'k--')

grid on
box on
xticks([])
xlim([-0.3, 0.3])
ylim_abs =  1.2 * max(abs(ch_diff));
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

scatter(0.1*randn(1, length(R_diff)), R_diff, 'k', 'filled')

plot([-0.3, 0.3], zeros(2,1), 'k', 'LineWidth', 2)
plot([-0.25, 0.25], median(R_diff)*ones(2,1), 'k--')

grid on
box on
xticks([])
xlim([-0.3, 0.3])
ylim_abs =  1.2 * max(abs(R_diff));
ylim([-ylim_abs, ylim_abs])
set(gca, 'YAxisLocation', 'right')
ylabel(['\DeltaR' newline '(Movie - Rest)'])

title(['Effect' newline 'Size'])

fontsize(gcf, 20, 'points')

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_summary_%s_movie_avg.png', ...
    fig_dir, signal_type), 'Resolution', 600)