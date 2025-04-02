%% Figure 1 Difference between Despicable Me English and resting fixation 

%% Figure settings
fig_font = 20;

model_dir = '../results/models_revision_1';

fig_dir = '../results/figures';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

font_size = 20;

%% Define patients
signal_type = 'HFA';

% Movie segment 
conditions = {'Resting_fixation', 'Despicable_Me_English_5min', 'Despicable_Me_English_last_5min', 'Inscapes_5min', 'Monkey_5min'};
% conditions = {'Resting_fixation_shift', 'Despicable_Me_English_5min_shift', ...
%     'Despicable_Me_English_last_5min_shift', 'Inscapes_5min_shift', 'Monkey_5min_shift'};

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient(table2array(sum(patient_list(:,2:3),2) == 2));

% Movie table for coordinates
load('../data/movie_subs_table.mat', 'movie_subs_table');

% Select significant channels or not
sig_channels = false;
p_thresh = 0.001;

res_file = sprintf('../data/movie_rest_stats_%s.mat');

if exist(res_file, 'file') == 0

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
    
    save(res_file, 'T_R', 'T_chn', 'patients')

else
    load(res_file, 'T_R', 'T_chn', 'patients')
end

Cond = table([1 2]','VariableNames',{'Conditions'});



%% Tests

%% Linear mixed effect model
% Remove outlier
N = length(patients);

R = [T_R.movies; T_R.rest(1:length(patients))];
chn = [T_chn.movies; T_chn.rest(1:length(patients))];
type = [T_R.movie_type; repmat({'Rest'}, length(patients), 1)];
stim = [ones(height(T_R),1); zeros(length(patients),1)];
subjects = repmat(patients, length(R) / length(patients), 1);

T = table(type,R,chn,subjects,stim);

fitlme(T,'R ~ stim  + (1|subjects)')
fitlme(T,'chn ~ stim  + (1|subjects)')

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

R_avg = mean(R_mat);
ch_avg = mean(ch_mat);

R_diff = R_mat - R_rest;
ch_diff = ch_avg - ch_rest;

%% Figure
R_plot = [R_rest; R_avg];
ch_plot = [ch_rest; ch_avg];

figure('Position', [1000,750,750,450])

% Significant Channels
tiledlayout(1,10);

nexttile(1,[1,4])
hold on

plot(ch_plot(:, diff(ch_plot) > 0), '.-', 'Color', [0.75, 0.3, 0], 'LineWidth', 2, 'MarkerSize', 20)
plot(ch_plot(:, diff(ch_plot) < 0), '.-', 'Color', [0, 0.3, 0.75], 'LineWidth', 2, 'MarkerSize', 20)

xticks([1, 2])
xticklabels({'Rest', 'Movies'})
xtickangle(45)
xlim([0.75, 2.25])

ylabel('Fraction of Channels')
set(gca, 'YAxisLocation', 'right')
ylim([0.9*min(ch_plot(:)), 1.1*max(ch_plot(:))])

fontsize(font_size, 'Points')
title(['Significant' newline 'Connections'])

grid on

% Effect size
nexttile(5,[1,4])
hold on 

plot(R_plot(:, diff(R_plot) > 0), '.-', 'Color', [0.75, 0.3, 0], 'LineWidth', 2, 'MarkerSize', 20)
plot(R_plot(:, diff(R_plot) < 0), '.-', 'Color', [0, 0.3, 0.75], 'LineWidth', 2, 'MarkerSize', 20)

xticks([1, 2])
xticklabels({'Rest', 'Movies'})
xtickangle(45)
xlim([0.75, 2.25])

ylabel('R')
set(gca, 'YAxisLocation', 'right')
ylim([0.95*min(R_plot(:)), 1.05*max(R_plot(:))])

fontsize(font_size, 'Points')
title('Effect Size')

grid on

legend(['Increase', repmat({''}, 1, length(R_plot)-2), 'Decrease'], ...
    'Position', [0.77,0.02,0.21,0.15]);

exportgraphics(gcf, sprintf('%s/fig4_connectivity_dme_rest_summary_%s_movie_avg.png', ...
    fig_dir, signal_type), 'Resolution', 600)