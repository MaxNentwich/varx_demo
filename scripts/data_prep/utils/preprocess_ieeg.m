%% Preprocess data

function preprocess_ieeg(options)

if options.visualize_preproc  
    plot_ch = 1;
end

for pat = 1:length(options.patients)

    data_files = dir(sprintf('%s/%s/matlab_data/Neural', options.data_dir, options.patients(pat).name));
    data_files([data_files.isdir]) = [];

    for v = 1:length(data_files)
        
        %% Check if the file was already preprocessed
        out_file = sprintf('%s/%s/%s/%s', options.data_dir, options.patients(pat).name, options.neural_prep_dir, data_files(v).name);
        
        % if exist(out_file, 'file') ~= 0
        %     continue
        % end
        
        fprintf('Processing subject %s, %s \n', options.patients(pat).name, data_files(v).name)

        %% Load ECoG data
        load(sprintf('%s/%s', data_files(v).folder, data_files(v).name), 'ieeg')

        fs = ieeg.fs;
        
        % Subject NS151 is missing data at the end for 'The Present'; pad this data        
        if sum(isnan(ieeg.data(:))) ~= 0 && strcmp(options.patients(pat).name, 'NS151')

            idx_nan = isnan(ieeg.data(1,:));

            ieeg.data(:, idx_nan) = repmat(ieeg.data(:, find(idx_nan, 1, 'first') - 1), 1, sum(isnan(ieeg.data(1,:))));

        end

        %% Adjust labels
        
        % If there are more labels than data files it can be assumed that there are reference channels not included 
        % in the data 
        if length(ieeg.label) > size(ieeg.data,1)  

            ieeg.label = ieeg.label(1:size(ieeg.data,1));   

            warning('Unnacounted labels without matching data channels! \n')

        % Otherwise, cut the last data channels that are not labeled 
        elseif length(ieeg.label) < size(ieeg.data,1)     

            ieeg.data = ieeg.data(1:length(ieeg.label), :);

            warning('Unnacounted data channels without matching labels! \n')

        end

        %% Patient NS178 has some single sample artifacts in the data -> interpolates

        if strcmp(options.patients(pat).name, 'NS178')
    
            samples_data = 1:size(ieeg.data,2);
    
            data_diff = mean(abs(diff(ieeg.data,[],2)));
    
            [~, idx_gap] = findpeaks(double(isoutlier(data_diff)));
            idx_gap = idx_gap+1;
            
            idx_good = setdiff(samples_data, idx_gap);
    
            for ch = 1:size(ieeg.data,1)
                ieeg.data(ch,idx_gap) = interp1(samples_data(idx_good), ieeg.data(ch, idx_good), idx_gap);
            end
    
        end
            
        % Remove bad channels
        load(sprintf('%s/%s/matlab_data/Bad_channels/bad_channels.mat', options.data_dir, options.patients(pat).name), ...
            'bad_chns_manual')
        
        bad_chns_manual(cellfun(@(C) strcmp(C, '[]'), bad_chns_manual)) = [];
        
        idx_bad = cellfun(@(C) find(ismember(ieeg.label, C)), bad_chns_manual, 'UniformOutput', false);
        idx_bad = cell2mat(idx_bad(cellfun(@(C) ~isempty(C), idx_bad)));
        
        ieeg.data(idx_bad,:) = [];
        ieeg.label(idx_bad) = [];
        
        ieeg.data(cellfun(@(C) strcmp(C, '[]'), ieeg.label),:) = [];
        ieeg.label(cellfun(@(C) strcmp(C, '[]'), ieeg.label)) = [];       

        %% Preprocess iEEG data

        % Plot before and after 
        if options.visualize_preproc 
            figure
            view_spec(ieeg.data(plot_ch,:), fs, 1, 300)
            
            figure
            plot(ieeg.data(plot_ch,:))         
        end
        
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
        end
        
        if options.visualize_preproc 
            figure
            view_spec(ieeg.data(plot_ch,:), fs, 1, 300)
            
            figure
            plot(ieeg.data(plot_ch,:))
        end
        
        %% Rerefercencing 
        % Uses Stefans code which is based on fieltrip structure
        
        ecog.ftrip.fsample = fs;
        ecog.ftrip.nChans = length(ieeg.label);
        ecog.ftrip.label = ieeg.label;
        ecog.ftrip.trial{1} = ieeg.data;
        ecog.ftrip.time = ieeg.time;
        
        ecog.ref.contributing_chans = [];
        cfg.includeChannels = get_indElecNames(ecog.ftrip.label);
        ecog = ecog_bipolarize(ecog, cfg);
        
        ieeg.label = ecog.ftrip.label;
        ieeg.data = ecog.ftrip.trial{1};
        
        clearvars ecog

        % Highpass filter
        for ch = 1:size(ieeg.data)
            ieeg.data(ch,:) = filtfilt(sos, g, double(ieeg.data(ch,:))')';
            ieeg.data(ch,:) = HpFilter(ieeg.data(ch,:), 5, 0.5, fs)';
        end

        if options.visualize_preproc 
            
            hold on
            plot(ieeg.data(plot_ch,:))
            
            figure
            view_spec(ieeg.data(plot_ch,:), fs, 1, 300)
         
            plot_label = ieeg.label(plot_ch);

        end
        
        %% Remove artfacts 
        % Subtract mean (might not be necessary after bipolar rereferencing 
        ieeg.data = bsxfun(@minus, ieeg.data, mean(ieeg.data,2));

        % Step function to close gaps
        step1 = ones(round(0.3*fs));
        step1 = step1/sum(step1);

        % Step funciton to extend regions 
        step2 = ones(round(0.05*fs));
        step2 = step2/sum(step2);

        % Window for smoothing edges
        h = hann(round(0.05*fs));

        win = zeros(size(ieeg.data'));

        for ch = 1:size(ieeg.data,1)

            % Connect gaps (region growing)
            win(:,ch) = conv(double(abs(ieeg.data(ch,:)) > 5*iqr(ieeg.data(ch,:))), step1, 'same');
            win(win(:,ch) > 0, ch) = 1; 

            win(:,ch) = conv(win(:,ch), step1, 'same');
            win(win(:,ch) < 0.99999, ch) = 0; 

            % Add step to extend over edges
            win(:,ch) = conv(win(:,ch), step2, 'same');
            win(win(:,ch) > 0, ch) = 1;     

            % Smooth edges
            win(:,ch) = 1 - conv(win(:,ch), h/sum(h), 'same');

        end

        % Apply window to data
        ieeg.data = ieeg.data'.*win;
        
        artifacts = win ~= 1;
        
        clearvars win
        
        if options.visualize_preproc 
            
            plot_idx = find(ismember(ieeg.label, plot_label));
            
            hold on
            plot(ieeg.data(:,plot_idx))
        
            figure
            view_spec(ieeg.data(:,plot_idx), fs, 1, 300)

        end
            
        %% Save the Data

        if exist(sprintf('%s/%s/matlab_data/Neural_prep', options.data_dir, options.patients(pat).name), 'dir') == 0
            mkdir(sprintf('%s/%s/matlab_data/Neural_prep', options.data_dir, options.patients(pat).name))
        end

        save(out_file, 'ieeg', '-v7.3')
        
        save(sprintf('%s/%s/matlab_data/Neural_prep/%s_artifacts.mat', options.data_dir, options.patients(pat).name, ... 
            strrep(data_files(v).name, '.mat', '')), 'artifacts', '-v7.3')
        
        close all
            
    end
 
end

end