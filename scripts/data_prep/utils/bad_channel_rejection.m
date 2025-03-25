%% Plot one recording per subject for visual inspection and identification of bad channels 

function bad_channel_rejection(options)

%% Loop to extract band power from each band 

for pat = 1:length(options.patients)
    
    data_files = dir(sprintf('%s/%s/matlab_data/Neural', options.data_dir, options.patients(pat).name));
    data_files([data_files.isdir]) = [];
    
    % Load the Present if available because it is shortest, otherwise Monkey2
    if sum(cellfun(@(C) contains(C, 'Present_Rep_2'), extractfield(data_files, 'name'))) == 1
        data_file_plot = data_files(cellfun(@(C) contains(C, 'Present_Rep_2'), extractfield(data_files, 'name')));
    elseif sum(cellfun(@(C) contains(C, 'Present_Rep_1'), extractfield(data_files, 'name'))) == 1
        data_file_plot = data_files(cellfun(@(C) contains(C, 'Present_Rep_1'), extractfield(data_files, 'name')));
    else 
        data_file_plot = data_files(cellfun(@(C) contains(C, 'Monkey2_Rep_1'), extractfield(data_files, 'name')));
    end
    
    % Output file
    out_file = sprintf('%s/Bad_channels/bad_channels.mat', strrep(data_file_plot.folder, 'Neural', ''));
    if exist(out_file, 'file') ~= 0, continue, end
    
    fprintf('Patient %s \n', options.patients(pat).name)
    
    %% Load ECoG data
    clearvars ieeg
    
    load(sprintf('%s/%s', data_file_plot.folder, data_file_plot.name), 'ieeg')
    
    if sum(isnan(ieeg.data(:))) ~= 0 && strcmp(options.patients(pat).name, 'NS151')
        
        idx_nan = isnan(ieeg.data(1,:));
        
        ieeg.data(:, idx_nan) = repmat(ieeg.data(:, find(idx_nan, 1, 'first') - 1), 1, sum(isnan(ieeg.data(1,:))));
        
    end

    %% Adjust labels

    % If there are more labels than data files it can be assumed that there are reference channels not included 
    % in the data 
    if length(ieeg.label) > size(ieeg.data,1)  
        
        ieeg.label = ieeg.label(1:size(ieeg.data,1));   
        
        warning('Unnacounted labels without matching data channels!')
        
    % Otherwise, cut the last data channels that are not labeled 
    elseif length(ieeg.label) < size(ieeg.data,1)     
        
        ieeg.data = ieeg.data(1:length(ieeg.label), :);
        
        warning('Unnacounted data channels without matching labels!')
        
    end

    %% Notch filter

    % Often collides with fieltrip function
    pathCell = regexp(path, pathsep, 'split');

    if any(strcmp('/home/max/Documents/MATLAB/Toolboxes/fieldtrip-20180426/external/signal', pathCell))
        rmpath('/home/max/Documents/MATLAB/Toolboxes/fieldtrip-20180426/external/signal')
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

    %% Plot data electrode by electrode
    idx_bad = [];
    
    elec_names = get_indElecNames(ieeg.label);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    for el = 1:length(elec_names)
        
        idx_elec = cellfun(@(C) strcmp(C, elec_names(el).Labels), ...
            cellfun(@(C) C(regexp(C, '\D')), ieeg.label, 'UniformOutput', false));
        
        elec_select = ieeg.label(idx_elec);
        
        if length(elec_select) > 16
            
            idx_block = find(idx_elec);
            
            for chunk = 1:ceil(length(elec_select)/16)
                
                elec_chunk = elec_select((chunk-1)*16+1 : chunk*16);
                
                idx_chunk = idx_block((chunk-1)*16+1 : chunk*16);
                
                for ch = 1:length(elec_chunk)            
                    elec_chunk{ch} = sprintf('%i %s', idx_chunk(ch), elec_chunk{ch});
                end
                
                %% Filter
                [z1,p1,k1] = butter(5, [59.5, 60.5]/(ieeg.fs/2), 'stop');
                [z2,p2,k2] = butter(5, [119.5, 120.5]/(ieeg.fs/2), 'stop');
                [z3,p3,k3] = butter(5, [179.5, 180.5]/(ieeg.fs/2), 'stop');

                Z = [z1;z2;z3];
                P = [p1;p2;p3];
                K = k1*k2*k3;

                [sos,g] = zp2sos(Z,P,K);

                ieeg.data(idx_chunk, :) = filtfilt(sos, g, double(ieeg.data(idx_chunk, :))')';

                %% Highpass filter
                ieeg.data(idx_chunk, :) = HpFilter(ieeg.data(idx_chunk, :), 5, 0.5, ieeg.fs)';
                
                %% Plot
                subplot(211)

                offset = -2.5e-4*[1:length(elec_chunk)]; 
                plot(bsxfun(@plus, ieeg.data(idx_chunk, :)', offset));

                yticks(fliplr(offset))
                yticklabels(flipud(elec_chunk))

                axis tight
                
                subplot(212)

                imagesc(ieeg.data(idx_chunk, :))

                caxis([-5e-5 5e-5])

                yticks(1:length(idx_chunk))
                yticklabels(ieeg.label(idx_chunk))
                
                % Prompt to enter bad channels
                bad_channels = input('Enter any bad channels (e.g 1, 3, [4,5], 3:5): ');
                
                idx_bad = [idx_bad, bad_channels];
            
            end
            
        else

            for ch = 1:length(elec_select)            
                elec_select{ch} = sprintf('%i %s', ch, elec_select{ch});
            end
            
            %% Filter
            [z1,p1,k1] = butter(5, [59.5, 60.5]/(ieeg.fs/2), 'stop');
            [z2,p2,k2] = butter(5, [119.5, 120.5]/(ieeg.fs/2), 'stop');
            [z3,p3,k3] = butter(5, [179.5, 180.5]/(ieeg.fs/2), 'stop');

            Z = [z1;z2;z3];
            P = [p1;p2;p3];
            K = k1*k2*k3;

            [sos,g] = zp2sos(Z,P,K);

            ieeg.data(idx_elec, :) = filtfilt(sos, g, double(ieeg.data(idx_elec, :))')';

            %% Highpass filter
            ieeg.data(idx_elec, :) = HpFilter(ieeg.data(idx_elec, :), 5, 0.5, ieeg.fs)';

            %% Plot

            subplot(211)

            offset = -2.5e-4*[1:sum(idx_elec)]; 
            plot(bsxfun(@plus, ieeg.data(idx_elec, :)', offset));

            yticks(fliplr(offset))
            yticklabels(flipud(elec_select))

            axis tight

            subplot(212)

            imagesc(ieeg.data(idx_elec, :))

            caxis([-5e-5 5e-5])

            yticks(1:length(idx_elec))
            yticklabels(ieeg.label(idx_elec))

            % Prompt to enter bad channels
            bad_channels = input('Enter any bad channels (e.g 1, 3, [4,5], 3:5): ');

            idx_channel = 1:length(ieeg.label);
            idx_elec_select = idx_channel(idx_elec);

            idx_bad = [idx_bad, idx_elec_select(bad_channels)];
            
        end
        
    end
    
    idx_bad = sort(idx_bad)';
    
    close all

    %% Save the bad channel data 
        
    bad_chns_manual = ieeg.label(idx_bad);

    if exist(sprintf('%s/Bad_channels', strrep(data_file_plot.folder, 'Neural', '')), 'dir') == 0
        mkdir(sprintf('%s/Bad_channels', strrep(data_file_plot.folder, 'Neural', '')))
    end
    
    save(out_file, 'bad_chns_manual')

end

end