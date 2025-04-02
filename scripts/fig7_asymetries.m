%% Difference between VAR and VARX models

%% Define patients
example_pat = 'NS127_02';
signal_type = 'HFA';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

fig_font = 24;

c_scale = 0.4;

fig_dir = '../results/figures';
model_dir = '../results/models_extend';
% model_dir = '../results/models_revision_1';

load('../data/movie_subs_table.mat', 'movie_subs_table');

addpath('../src')

% Choose annotation ['myelin', 'timescale', 'gradient_1', 'gradient_2']
annot_name = 'myelin';

% Normalize asymetry by areas
norm_areas = false;

% Poster settings
poster_size = false;

if poster_size
    fig_font = 26;
end

%% Run python script to compute T1w/T2w ratio (myelination) in each parcel of the Desikan-Killiany atlas
% Tricky setup -> ran file in spyder
% pyenv('Version','/home/max/anaconda3/envs/varx_plots/bin/python3.8')
% pyrunfile('hierarchy_dk_atlas.py')

%% Load data

% Initialize array to collect data
varx_R_asym = cell(1, length(patients));
loc = cell(1, length(patients));
direction_contact = cell(1, length(patients));
pat_name = cell(1, length(patients));

for pat = 1:length(patients)
    
    if exist(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat})) == 0
        continue
    end

    if strcmp(signal_type, 'LFP')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx', 'labels')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('%s/%s_varx_models.mat', model_dir, patients{pat}), 'm_varx_hfa', 'labels')
        m_varx = m_varx_hfa;
    end

    varx_Rvalue = sqrt(1-exp(-m_varx.A_Deviance/m_varx.T));
    varx_R_asym{pat} = abs(varx_Rvalue) - abs(varx_Rvalue)';

    % Get electrode labels and coordinates
    if strcmp(patients{pat}, 'NS174_02')
        pat_coord = 'NS174';
    else
        pat_coord = patients{pat};
    end

    movie_table_pat = movie_subs_table(ismember(movie_subs_table.SubID, pat_coord), :);

    labels_coord =  cellfun(@(C) sprintf('%s_%s', pat_coord, C), labels, 'UniformOutput', false);
    [loc{pat}, ~, ~] = localize_elecs_bipolar(labels_coord, 'dk');
    loc{pat} = loc{pat}';

    pat_name{pat} = repmat({patients{pat}}, 1, length(labels_coord));

    % Normalize number of electrodes
    if norm_areas

        locs_unique = unique(loc{pat});
        varx_R_asym_norm = zeros(length(locs_unique), size(varx_R_asym{pat},2));
    
        for l = 1:length(locs_unique)
            varx_R_asym_norm(l,:) = mean(varx_R_asym{pat}(ismember(loc{pat}, locs_unique{l}), :));
        end

        direction_contact{pat} = mean(varx_R_asym_norm);

    else
        direction_contact{pat} = mean(varx_R_asym{pat});
    end

end

%% Direction by brain area
loc_all = cat(2, loc{:});
direction_contact = cell2mat(direction_contact);
pat_name = [pat_name{:}];

areas = unique(loc_all);
areas(ismember(areas, {'accumbens', ...
    'brainstem', ...
    'caudate', ...
    'choroidplexus', ...
    'unknown', ...
    'ventricle', ...
    'whitematter'})) = [];

direction_area = cell(length(areas), 1);
median_area = nan(length(areas), 1);
patient_area = cell(length(areas), 1);
pat_u = cell(length(areas), 1);
direction_area_pat = cell(length(areas), 1);

for a = 1:length(areas)

    direction_area{a} = direction_contact(ismember(loc_all, areas{a}));
    median_area(a) = median(direction_area{a});

    patient_area{a} = pat_name(ismember(loc_all, areas{a}));

    pat_u{a} = unique(patient_area{a});
    direction_area_pat{a} = zeros(1, length(pat_u{a}));

    for u = 1:length(pat_u{a})
        direction_area_pat{a}(u) = median(direction_area{a}(ismember(patient_area{a}, pat_u{a}{u})));
    end

end

[~,idx_sort] = sort(median_area, 'descend');
areas = areas(idx_sort);
median_area = median_area(idx_sort);
patient_area = patient_area(idx_sort);
pat_u = pat_u(idx_sort);
direction_area = direction_area(idx_sort);
direction_area_pat = direction_area_pat(idx_sort);

%% Correlate with myelination (cortical hierarchy)
annot_data = readtable(sprintf('../data/%s_lh_parcels_aparc.xlsx', annot_name));

map_areas = table2array(annot_data(:,2));
map_vals = table2array(annot_data(:,3));

map_areas = cellfun(@(C) strrep(C, 'b''', ''), map_areas, 'UniformOutput', false);
map_areas = cellfun(@(C) strrep(C, '''', ''), map_areas, 'UniformOutput', false);

idx_areas = cellfun(@(C) find(ismember(areas, C)), map_areas, 'UniformOutput', false);

idx_empty = cellfun(@(C) isempty(C), idx_areas);
map_areas = map_areas(~idx_empty);
map_vals = map_vals(~idx_empty);
idx_areas = cell2mat(idx_areas(~idx_empty));

areas = areas(idx_areas);
median_area = median_area(idx_areas);
patient_area = patient_area(idx_areas);
pat_u = pat_u(idx_areas);
direction_area = direction_area(idx_areas);
direction_area_pat = direction_area_pat(idx_areas);

error_area = cellfun(@(C) std(C)/sqrt(length(C)), direction_area_pat);

% Save
varx_dir = table(areas', median_area, 'VariableNames', {'parcel_name', 'varx_dir'});
writetable(varx_dir, '../data/varx_dir_aparc.xlsx')

%% Plot
X = [ones(length(map_vals),1) map_vals];
b = X\median_area;
median_area_h = X*b;

if poster_size 
    figure('Units', 'inches', 'Position', [1,1,3.5,3.5])
else
    figure()
end

hold on 

error_color = 0.5*ones(1,3);

for m = 1:length(map_vals)
    plot(map_vals(m)*[1,1], median_area(m) + error_area(m)*[-1,1], ...
        'Color', error_color, 'LineWidth', 1)
    plot(map_vals(m)*[0.995, 1.005], (median_area(m)-error_area(m))*[1,1], ...
        'Color', error_color, 'LineWidth', 1)
    plot(map_vals(m)*[0.995, 1.005], (median_area(m)+error_area(m))*[1,1], ...
        'Color', error_color, 'LineWidth', 1)
end

scatter(map_vals, median_area, 'filled', 'k')
plot(map_vals, median_area_h, 'r')

xlabel('T1w/T2w ratio')
ylabel('Mean R-R^T')
xlim([min(map_vals)-0.01, max(map_vals)+0.01])
fontsize(gcf, fig_font, 'points')
set(gca, 'XDir', 'reverse')

exportgraphics(gcf, sprintf('%s/fig7_direction_vs_%s_%s.png', fig_dir, annot_name, signal_type), 'Resolution', 300)

[r_area, p_area] = corr(map_vals, median_area);

fprintf('r=%1.3f, p=%1.3f\n', r_area, p_area)

% % Without oultier
% idx_out = find(ismember(areas, 'transversetemporal'));
% idx_keep = setdiff(1:length(areas), idx_out);
% 
% [r_area, p_area] = corr(map_vals(idx_keep), median_area(idx_keep));
% fprintf('r=%1.3f, p=%1.3f\n', r_area, p_area)

% Any patients outliers?
% figure('Position', [1000,750,675,450])
% plot(direction_area_pat, 'o-')
% xticklabels(cellfun(@(C) strrep(C, '_', ' '), pat_u{idx_out}, 'UniformOutput', false))
% ylabel('Mean R-R^T')
% xtickangle(45)
% title(areas{idx_out})
% grid on
% fontsize(gcf, fig_font, 'points')
% 
% exportgraphics(gcf, sprintf('%s/fig7_direction_vs_%s_%s_outlier_per_patient.png', fig_dir, annot_name, signal_type), 'Resolution', 300)

%% Linear mixed effect model
area_labels = {};
map_all = [];
for a = 1:length(areas)
    area_labels = [area_labels; repmat({areas{a}}, length(direction_area_pat{a}), 1)];
    map_all = [map_all; repmat(map_vals(a), length(direction_area_pat{a}), 1)];
end

direction_area_pat = cell2mat(direction_area_pat')';
patient_labels = [pat_u{:}]';

T = table(patient_labels, area_labels, direction_area_pat, map_all, ...
    'VariableNames', {'patient', 'area', 'direction', 'gradient'});

fitlme(T, 'direction ~ gradient  + (1|patient)')

save(sprintf('../results/stats/asymetry_stats_%s.mat', signal_type), 'T')

%% Plot example of matrix for one patient -> order by hierarchy

% Get data from a sample patient
id_ex = ismember(patients, example_pat);

varx_asym_ex = varx_R_asym{id_ex};
asym_max = max(abs(varx_asym_ex(:)));

loc_ex = loc{id_ex};

% Sort areas by hierachy
[~, idx_myelin] = sort(map_vals, 'descend');
areas_sort = map_areas(idx_myelin);

idx_sort_roi = cellfun(@(C) find(ismember(areas_sort, C)), loc_ex, 'UniformOutput', false);
idx_empty = cellfun(@(C) isempty(C), idx_sort_roi);

varx_asym_ex = varx_asym_ex(~idx_empty, ~idx_empty);
loc_ex = loc_ex(~idx_empty);
idx_sort_roi = idx_sort_roi(~idx_empty);

idx_sort_roi = cell2mat(idx_sort_roi);
[~,idx_sort] = sort(idx_sort_roi);

varx_asym_ex = varx_asym_ex(idx_sort, idx_sort);
loc_ex = loc_ex(idx_sort);

bnd = find(diff(idx_sort_roi(idx_sort)) ~= 0);

% Plot
if poster_size 
    figure('Units', 'inches', 'Position', [1,1,7,4])
else
    figure('Position', [350,250,675,450])
end

imagesc(varx_asym_ex)
hold on

y_lim = ylim();
x_lim = xlim();

for b = 1:length(bnd)
    plot([bnd(b)+0.5, bnd(b)+0.5], y_lim, 'k--')
end

clim([-c_scale*asym_max, c_scale*asym_max])
axis square
colormap(slanCM('bwr'))
xlabel('Cause y(t-1)')
ylabel('Effect on y(t)')
xticks([])
yticks([])
cb = colorbar(); 
ylabel(cb,'R - R^T' ,'Rotation',270)
axis tight
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig7_example_connectivity_%s.png', fig_dir, signal_type), 'Resolution', 300)

%% Plot spatial maps
% Also didn't run through matlab
% pyrunfile('plot_varx_dir_myelin_dk.py')
