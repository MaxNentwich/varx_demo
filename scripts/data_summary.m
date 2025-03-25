%% Data summary 

data_dir = '/media/max/Workspace/Data/varx_data';

patients = dir(data_dir);
patients = patients(3:end);

rec_names = cell(length(patients),1);
length_rec = NaN(length(patients),1);
n_chns = NaN(length(patients),1);

duplicate_names = {'The_Present.mat', 'Monkey1.mat', 'Monkey2.mat', 'Monkey5.mat', 'Despicable_Me_English_5min.mat'};

%% Load data
for pat = 1:length(patients)

    recordings = dir(sprintf('%s/%s', data_dir, patients(pat).name));
    recordings([recordings.isdir]) = [];

    recordings(ismember({recordings.name}, duplicate_names)) = [];

    rec_names{pat} = {recordings.name};
    length_rec_pat = NaN(1,length(recordings));

    for rec = 1:length(recordings)

        load(sprintf('%s/%s/%s', data_dir, patients(pat).name, recordings(rec).name), 'lfp', 'fs');

        length_rec_pat(rec) = size(lfp,1) / fs ;
        n_chns(pat) = size(lfp,2);

    end

    length_rec(pat) = sum(length_rec_pat) / 60;

end

%% Load demographics
demographics = readtable('../data/demographics_ns.xlsx');
demographics(~cellfun(@(C) contains(C, 'NS'), demographics.PatientID), :) = [];

% Align with patient names
age = NaN(length(patients),1);
sex = cell(length(patients),1);
for pat = 1:length(patients)
    age(pat) = demographics.Age(cellfun(@(C) strcmp(C, patients(pat).name), demographics.PatientID));
    sex{pat} = demographics.Sex(cellfun(@(C) strcmp(C, patients(pat).name), demographics.PatientID));
end

%% Organize

% Find names of all recordings
rec_names_all = [];
for pat = 1:length(patients)
    rec_names_all = [rec_names_all, rec_names{pat}];
end

rec_names_unique = unique(rec_names_all);

% Don't need to consider 5 min version of Despicable Me here
rec_names_unique(cellfun(@(C) contains(C, '5min'), rec_names_unique)) = [];

%% Collect list of movies recorded for each patient, total duration, and number of channels 
rec_mat = zeros(length(patients), length(rec_names_unique));

for pat = 1:length(patients)
    rec_mat(pat,:) = ismember(rec_names_unique, rec_names{pat});
end

%% Reformat patients ID
pat_name = {patients.name};

pat_name_unique = cell(length(pat_name),1);
for pat = 1:length(pat_name)
    name_split = strsplit(pat_name{pat}, '_');
    pat_name_unique{pat} = name_split{1};
end

pat_name_unique = unique(pat_name_unique);

% Create new id's and removing 3 and 4 to match previous paper
pat_id = cell(length(pat_name_unique)+2,1);
for id = 1:length(pat_name_unique)+2
    pat_id{id} = sprintf('Pat_%d', id);
end
pat_id(3:4) = [];

patient_id = cell(length(pat_name),1);
for pat = 1:length(pat_name)
    
    patient_id{pat} = pat_id{cellfun(@(C) contains(pat_name{pat}, C), pat_name_unique)};

    name_split = strsplit(pat_name{pat}, '_');
    if length(name_split) == 2 && pat > 2
        patient_id{pat} = sprintf('%s_%s', patient_id{pat}, name_split{2});
    end

end

%% Create a table
rec_names_unique = cellfun(@(C) strrep(C, '.mat', ''), rec_names_unique, 'UniformOutput', false);
rec_names_unique = cellfun(@(C) strrep(C, '_', ' '), rec_names_unique, 'UniformOutput', false);

rec_mat = [rec_mat, length_rec, n_chns];

rec_table = array2table(rec_mat, 'VariableNames', [rec_names_unique, 'Total Length of Recordings [min]', 'Number of Channels']);
rec_table.Age = age;
rec_table.Sex = sex;
rec_table.patient_id = patient_id;

rec_table = renamevars(rec_table, rec_table.Properties.VariableNames, [rec_table.Properties.VariableNames(1:end-1), 'Patient ID']);

% Reorder columns
idx_id = find(ismember(rec_table.Properties.VariableNames, 'Patient ID'));
idx_age = find(ismember(rec_table.Properties.VariableNames, 'Age'));
idx_sex = find(ismember(rec_table.Properties.VariableNames, 'Sex'));

idx_reorder = [idx_id, idx_age, idx_sex, setdiff(1:length(rec_table.Properties.VariableNames), [idx_id, idx_age, idx_sex])];

rec_table = rec_table(:,idx_reorder);

% Save the table in CSV
writetable(rec_table, '../data/recording_summary.csv')  
