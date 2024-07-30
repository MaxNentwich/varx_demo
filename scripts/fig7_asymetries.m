%% Difference between VAR and VARX models

%% Define patients
example_pat = 'NS127_02';
signal_type = 'LFP';

patient_list = readtable('../data/varx_patient_list.xlsx');
patients = patient_list.Patient;

fig_font = 16;

c_scale = 0.4;

fig_dir = '../results/figures';

load('../data/movie_subs_table.mat', 'movie_subs_table');

addpath('../src')

% Poster settings
poster_size = true;

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

for pat = 1:length(patients)
    
    if strcmp(signal_type, 'LFP')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx', 'labels')
    elseif strcmp(signal_type, 'HFA')
        load(sprintf('../results/models/%s_varx_models.mat', patients{pat}), 'm_varx_hfa', 'labels')
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

    direction_contact{pat} = mean(varx_R_asym{pat});

end

%% Direction by brain area
loc_all = cat(2, loc{:});
direction_contact = cell2mat(direction_contact);

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

for a = 1:length(areas)
    direction_area{a} = direction_contact(ismember(loc_all, areas{a}));
    median_area(a) = median(direction_area{a});
end

[~,idx_sort] = sort(median_area, 'descend');
areas = areas(idx_sort);
median_area = median_area(idx_sort);

%% Correlate with myelination (cortical hierarchy)
myelin = readtable('../data/myelin_lh_parcels_aparc.xlsx');

myelin_areas = table2array(myelin(:,2));
myelin_vals = table2array(myelin(:,3));

myelin_areas = cellfun(@(C) strrep(C, 'b''', ''), myelin_areas, 'UniformOutput', false);
myelin_areas = cellfun(@(C) strrep(C, '''', ''), myelin_areas, 'UniformOutput', false);

idx_areas = cellfun(@(C) find(ismember(areas, C)), myelin_areas, 'UniformOutput', false);

idx_empty = cellfun(@(C) isempty(C), idx_areas);
myelin_areas = myelin_areas(~idx_empty);
myelin_vals = myelin_vals(~idx_empty);
idx_areas = cell2mat(idx_areas(~idx_empty));

areas = areas(idx_areas);
median_area = median_area(idx_areas);

% Save
varx_dir = table(areas', median_area, 'VariableNames', {'parcel_name', 'varx_dir'});
writetable(varx_dir, '../data/varx_dir_aparc.xlsx')

%% Plot
X = [ones(length(myelin_vals),1) myelin_vals];
b = X\median_area;
median_area_h = X*b;

if poster_size 
    figure('Units', 'inches', 'Position', [1,1,3.5,3.5])
else
    figure()
end

hold on 
scatter(myelin_vals, median_area, 'filled')
plot(myelin_vals,median_area_h,'k')

xlabel('T1w/T2w ratio')
ylabel('Mean R-R^T')
grid on
fontsize(gcf, fig_font, 'points')
set(gca, 'XDir', 'reverse')

exportgraphics(gcf, sprintf('%s/fig7_direction_vs_t1wt2w_%s.png', fig_dir, signal_type), 'Resolution', 300)

[r_area, p_area] = corr(myelin_vals, median_area);

fprintf('r=%1.3f, p=%1.3f\n', r_area, p_area)

%% Plot example of matrix for one patient -> order by hierarchy

% Get data from a sample patient
id_ex = ismember(patients, example_pat);

varx_asym_ex = varx_R_asym{id_ex};
asym_max = max(abs(varx_asym_ex(:)));

loc_ex = loc{id_ex};

% Sort areas by hierachy
[~, idx_myelin] = sort(myelin_vals, 'descend');
areas_sort = myelin_areas(idx_myelin);

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
    figure()
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

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/fig7_example_connectivity_%s.png', fig_dir, signal_type), 'Resolution', 300)

%% Plot spatial maps
% Also didn't run through matlab
% pyrunfile('plot_varx_dir_myelin_dk.py')
