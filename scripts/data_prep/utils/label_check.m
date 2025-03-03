%% Check if the channel labels look right and find patients where they may need to be corrected manually

function label_check(options)

for p = 1:length(options.patients)
    
    % Output file
    out_dir = sprintf('%s/%s/matlab_data/Bad_channels', options.data_dir, options.patients(p).name);
    out_file = sprintf('%s/channel_labels.mat', out_dir);
    if exist(out_file, 'file') ~= 0, continue, end

    % Load the shortest file
    files = dir(sprintf('%s/%s/%s', options.data_dir, options.patients(p).name, options.neural_dir));
    files([files.isdir]) = [];
    
    [~, idx_small] = min([files.bytes]);
    load(sprintf('%s/%s/%s/%s', options.data_dir, options.patients(p).name, options.neural_dir, files(idx_small).name), 'ieeg')
    
    if length(ieeg.label) > size(ieeg.data,1)
        
        labels_unassigend = ieeg.label(size(ieeg.data,1)+1:end);
        elecs_unassigned = unique(cellfun(@(C) C(regexp(C, '\D')), labels_unassigend, 'UniformOutput', false));
        
        for i = 1:length(elecs_unassigned)
            fprintf('Electrode not found in data: %s ... removed\n', elecs_unassigned{i})
        end
        
        ieeg.label = ieeg.label(1:size(ieeg.data,1));
        
    end
    
    %% Filter the data
    ieeg.data(isnan(ieeg.data)) = 0;
    fs = ieeg.fs;
    
    %% Notch filter
    [z1,p1,k1] = butter(5, [59.5, 60.5]/(fs/2), 'stop');
    [z2,p2,k2] = butter(5, [119.5, 120.5]/(fs/2), 'stop');
    [z3,p3,k3] = butter(5, [179.5, 180.5]/(fs/2), 'stop');

    Z = [z1;z2;z3];
    P = [p1;p2;p3];
    K = k1*k2*k3;

    [sos,g] = zp2sos(Z,P,K);

    for ch = 1:size(ieeg.data)
        ieeg.data(ch,:) = filtfilt(sos, g, double(ieeg.data(ch,:))')';
        ieeg.data(ch,:) = HpFilter(ieeg.data(ch,:), 5, 0.5, fs)';
    end

    ieeg.data = ieeg.data - mean(ieeg.data);
    
    ieeg.data(abs(ieeg.data) > 5*iqr(ieeg.data(:))) = 0;
    
    corr_mat = corr(ieeg.data');
    
    figure
    imagesc(corr_mat)
    hold on

    % Find the labels of all electrode shafts
    elec_labels = cellfun(@(C) C(regexp(C, '\D')), ieeg.label, 'UniformOutput', false);
    elec_names = unique(elec_labels, 'stable');

    % Label each electrode
    idx_elec = zeros(size(ieeg.label));

    for e = 1:length(elec_names)
        idx_elec(cellfun(@(C) strcmp(C, elec_names{e}), elec_labels)) = e;
    end

    % 'Boundaries' between electrodes
    elec_bnd = find(diff(idx_elec) ~= 0);
    elec_bnd = [elec_bnd; length(idx_elec)];

    for b = 1:length(elec_bnd)
        plot([elec_bnd(b) elec_bnd(b)]+0.5, ylim, 'r', 'LineWidth', 2)
        plot(xlim, [elec_bnd(b) elec_bnd(b)]+0.5, 'r', 'LineWidth', 2)
    end

    % Add the names of electrodes
    elec_bnd  = [0; elec_bnd];
    n_elec = elec_bnd(2:end) - elec_bnd(1:end-1);
    label_pos = elec_bnd(1:end-1) + (n_elec/2) + 0.5;

    xticks(label_pos)
    yticks(label_pos)

    xticklabels(elec_names)
    yticklabels(elec_names)
    
    title(strrep(options.patients(p).name, '_', ' '))
  
    mean_corr = [];
    
    for i = 1:length(elec_bnd)-1
        for j = 1:length(elec_bnd)-1
            mean_corr(i,j) = mean2(corr_mat(elec_bnd(i)+1:elec_bnd(i+1), elec_bnd(j)+1:elec_bnd(j+1)));
        end
    end
    
    figure 
    imagesc(mean_corr)
    caxis([0 0.5])
    
    xticks(1:length(elec_names))
    yticks(1:length(elec_names))
    
    xticklabels(elec_names)
    yticklabels(elec_names)
    
    title(strrep(options.patients(p).name, '_', ' '))
    
    figure
    hold on
    [~, idx_max] = max(mean_corr);
    plot(1:size(mean_corr,1))
    plot(idx_max, '*')
    
    xticks(1:length(elec_names))
    yticks(1:length(elec_names))
    
    xticklabels(elec_names)
    yticklabels(elec_names)
    
    title(strrep(options.patients(p).name, '_', ' '))

    %% Manual input if data looks alright
    correct_input = 0;
    
    while ~correct_input               
        check_var = input('Do the channel labels look reasonable? (yes=1/no=0) ');
        if ~ismember(check_var, [0, 1])
            warning('Only 0 (no) and 1(yes) are allowed!\n')
        else
            correct_input = 1;
        end
    end
        
    if check_var 
        if exist(out_dir, 'dir') == 0, mkdir(out_dir), end
        save(out_file, 'check_var')
    elseif ~check_var
        error('Labels might not be correcly matched to data! Double check the anatomy folder!\n')
    end
    
    close all

end

end