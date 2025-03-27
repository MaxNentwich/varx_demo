%% Compare eye movements between movies and rest

% Figure settings
fig_font = 16;
fig_dir = '../results/figures';

% Source directories
data_dir = '/media/max/e61df479-8f57-4855-b852-04ccdfb12a6c/ECoGData/Tobii/Patients';
saccade_dir = '/media/max/Workspace/Data/saccade_data';

% List of patients 
patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

% Video file list
vids_select = {'Despicable_Me_English.mat', 'Monkey1_Rep_1.mat', 'Monkey2_Rep_1.mat', 'Inscapes.mat', 'Resting_fixation.mat'};

saccade_interval = cell(length(patients), length(vids_select));
saccade_amp = cell(length(patients), length(vids_select));
fixation_position = cell(length(patients), length(vids_select));
rec_time = cell(length(patients), length(vids_select));

% 'Screens' for eye position heatmap
screen_size = [1080, 1920];

screen_dme = zeros(screen_size);
screen_monkey = zeros(screen_size);
screen_inscapes = zeros(screen_size);
screen_rest = zeros(screen_size);

for i = 1:length(patients)

    data_pat = sprintf('%s/%s/matlab_data/Envelope_phase/BHA', data_dir, patients{i});
    
    vids = dir(data_pat);
    vids([vids.isdir]) = [];
    
    idx_vids = ismember({vids.name}, vids_select);
    vids = vids(idx_vids);

    for v = 1:length(vids)
    
        % Load eyetracking data
        load(sprintf('%s/%s_%s', saccade_dir, patients{i}, vids(v).name), 'eye', 'saccade_onset', 'saccade_amplitude', 'pos_post')

        idx_vid = ismember(vids_select, vids(v).name);

        saccade_interval{i,idx_vid} = diff(eye.time(find(saccade_onset)))';
        saccade_amp{i,idx_vid} = saccade_amplitude;
        fixation_position{i,idx_vid} = pos_post;
        rec_time{i,idx_vid} = eye.time(end);

        for p = 1:size(pos_post,1)
            
            pos_pix = round(pos_post(p,:) .* fliplr(screen_size));

            if strcmp('Despicable_Me_English.mat', vids(v).name)
                screen_dme(pos_pix(2), pos_pix(1)) = 1;
            elseif strcmp('Monkey1_Rep_1.mat', vids(v).name) || strcmp('Monkey2_Rep_1.mat', vids(v).name) 
                screen_monkey(pos_pix(2), pos_pix(1)) = 1;
            elseif strcmp('Inscapes.mat', vids(v).name)
                screen_inscapes(pos_pix(2), pos_pix(1)) = 1;
            elseif strcmp('Resting_fixation.mat', vids(v).name)
                screen_rest(pos_pix(2), pos_pix(1)) = 1;
            end

        end

    end

end

% Organize
sacc_int_dme = cell2mat(saccade_interval(:, ismember(vids_select, 'Despicable_Me_English.mat')));
sacc_int_monkey = [cell2mat(saccade_interval(:,ismember(vids_select, 'Monkey1_Rep_1.mat'))); ...
    cell2mat(saccade_interval(:,ismember(vids_select, 'Monkey2_Rep_1.mat')))];
sacc_int_inscapes = cell2mat(saccade_interval(:, ismember(vids_select, 'Inscapes.mat')));
sacc_int_rest = cell2mat(saccade_interval(:, ismember(vids_select, 'Resting_fixation.mat')));

sacc_amp_dme = cell2mat(saccade_amp(:, ismember(vids_select, 'Despicable_Me_English.mat')));
sacc_amp_monkey = [cell2mat(saccade_amp(:,ismember(vids_select, 'Monkey1_Rep_1.mat'))); ...
    cell2mat(saccade_amp(:,ismember(vids_select, 'Monkey2_Rep_1.mat')))];
sacc_amp_inscapes = cell2mat(saccade_amp(:, ismember(vids_select, 'Inscapes.mat')));
sacc_amp_rest = cell2mat(saccade_amp(:, ismember(vids_select, 'Resting_fixation.mat')));

% Overall saccade frequency and number
n_dme = length(sacc_amp_dme);
n_monkey = length(sacc_amp_monkey);
n_inscapes = length(sacc_amp_inscapes);
n_rest = length(sacc_amp_rest);

%% Inter-saccade-interval
bin_scale = 30;

figure('Position', [500,325,675,350])
hold on

histogram(1e3*sacc_int_dme, int64(bin_scale*max(sacc_int_dme)), 'Normalization', 'probability')
histogram(1e3*sacc_int_monkey, int64(bin_scale*max(sacc_int_monkey)), 'Normalization', 'probability')
histogram(1e3*sacc_int_inscapes, int64(bin_scale*max(sacc_int_inscapes)), 'Normalization', 'probability')
histogram(1e3*sacc_int_rest, int64(bin_scale*max(sacc_int_rest)), 'Normalization', 'probability')

xlim([0,1000])
grid on
xlabel('Inter-saccade-interval [ms]')
ylabel('Probability')
legend({'Despicable Me English', 'Monkey', 'Inscapes', 'Rest'})
fontsize(gcf(), fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig_saccade_params_isi.png', fig_dir), 'Resolution', 300)

%% Saccade amplitude
bin_scale = 4;

figure('Position', [500,375,550,350])
hold on

histogram(sacc_amp_dme, int64(bin_scale*max(sacc_amp_dme)), 'Normalization', 'probability')
histogram(sacc_amp_monkey, int64(bin_scale*max(sacc_amp_monkey)), 'Normalization', 'probability')
histogram(sacc_amp_inscapes, int64(bin_scale*max(sacc_amp_inscapes)), 'Normalization', 'probability')
histogram(sacc_amp_rest, int64(bin_scale*max(sacc_amp_rest)), 'Normalization', 'probability')

xlim([0,15])
grid on
xlabel('Saccade amplitude [DVA]')
ylabel('Probability')
legend({'Despicable Me English', 'Monkey', 'Inscapes', 'Rest'})
fontsize(gcf(), fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig_saccade_params_amplitude.png', fig_dir), 'Resolution', 300)

%% Heatmaps of fixation position 

% Despicable Me
figure('Position', [525,350,550,350])
imagesc(imgaussfilt(screen_dme, 20))
xticks([])
yticks([])
title('Despicable Me English')
fontsize(gcf(), fig_font, 'points')
exportgraphics(gcf, sprintf('%s/fig_saccade_params_heatmap_dme.png', fig_dir), 'Resolution', 300)

% Monkeys
figure('Position', [525,350,550,350])
imagesc(imgaussfilt(screen_monkey, 20))
xticks([])
yticks([])
title('Monkey')
fontsize(gcf(), fig_font, 'points')
exportgraphics(gcf, sprintf('%s/fig_saccade_params_heatmap_monkey.png', fig_dir), 'Resolution', 300)

% Inscaples
figure('Position', [525,350,550,350])
imagesc(imgaussfilt(screen_inscapes, 20))
xticks([])
yticks([])
title('Inscapes')
fontsize(gcf(), fig_font, 'points')
exportgraphics(gcf, sprintf('%s/fig_saccade_params_heatmap_inscapes.png', fig_dir), 'Resolution', 300)

% Resting State
figure('Position', [525,350,550,350])
imagesc(imgaussfilt(screen_rest, 20))
xticks([])
yticks([])
title('Resting State')
fontsize(gcf(), fig_font, 'points')
exportgraphics(gcf, sprintf('%s/fig_saccade_params_heatmap_rest.png', fig_dir), 'Resolution', 300)
