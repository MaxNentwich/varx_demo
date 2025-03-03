% Load NS files from TDT folders

function ieeg = load_NS(data_dir, patient_name, task, ecog_folder)

    directory = sprintf('%s/%s/Neural/%s/%s', data_dir, patient_name, task, ecog_folder);
    
    % Excel file with channel labels and definitions (which are bad, outside brain, within soz etc.)
    label_file = sprintf('%s/%s/Anatomy/%s_Electrode_Labels_TDT.xlsx', data_dir, patient_name, patient_name);
    
    % File with analog data labels
    analog_chn_file = sprintf('%s/%s/%s_ExtChannels.xlsx', data_dir, patient_name, patient_name);
       
    %% Import TDT    
    EEG1 = TDTbin2mat(directory, 'STORE', 'EEG1');
    EEG2 = TDTbin2mat(directory, 'STORE', 'EEG2');
    Wav5 = TDTbin2mat(directory, 'STORE', 'Wav5');
    
    if isempty(EEG1.streams)
        EEG1 = TDTbin2mat(directory, 'STORE', 'RAWx');
        EEG1.streams.EEG1.data = EEG1.streams.RAWx.data;
        EEG1.streams.EEG1.fs = EEG1.streams.RAWx.fs;
    end
        
    % Save info parameter
    ieeg.params = EEG1.info;
   
    % Load headers
    data_hdr = TDTbin2mat(directory, 'HEADERS', 1);
    
    %% Import labels
    label_import = import_label(label_file);
    
    %% Import all depth and strips/grids, based on xls sheet.
    ieeg.data = EEG1.streams.EEG1.data;         
    ieeg.fs = EEG1.streams.EEG1.fs;  
     
    clearvars EEG1

    % Combine first and second block
    if ~isempty(EEG2.streams)
        
        data2 = EEG2.streams.EEG2.data;
        fs2 = EEG2.streams.EEG2.fs;
        
        clearvars EEG2
        
        if fs2 ~= ieeg.fs
            error('Sampling rates of EEG data blocks are different! \n')
        end

        if size(data2, 2) == size(ieeg.data, 2)
            ieeg.data = [ieeg.data; data2]; 

        % Pad or cut data in second block if necessary
        elseif size(ieeg.data, 2) > size(data2, 2)
            ieeg.data = [ieeg.data; [data2, nan(size(data2,1), size(ieeg.data,2) - size(data2,2))]];

        elseif size(ieeg.data, 2) < size(data2, 2)
            ieeg.data = [ieeg.data; data2(:, 1:size(ieeg.data,2))];

        end
    
    end

    clearvars data1 data2
            
    % Channels labeled as bad, outside the skull, seizure onset zone 
    ieeg.szr_onset_chans = label_import.Label(label_import.SOZ == 1);
    ieeg.bad_chans = [label_import.Label(label_import.BAD == 1); label_import.Label(label_import.Out == 1)];
    ieeg.spike_chans = label_import.Label(label_import.Spikey == 1);
 
    ieeg.label = label_import.Label;
    
    % Triggers 
    if isfield(data_hdr.stores, 'PtC2')
        ieeg.triggers = unique([data_hdr.stores.PtC2.onset, data_hdr.stores.PtC4.onset, data_hdr.stores.PtC6.onset])';
    elseif isfield(data_hdr.stores, 'PC0_')
        ieeg.triggers = unique([data_hdr.stores.PC0_.onset, data_hdr.stores.PtC4.onset, data_hdr.stores.PtC6.onset])';
    elseif isfield(data_hdr.stores, 'PC2_')
        ieeg.triggers = unique([data_hdr.stores.PC2_.onset, data_hdr.stores.PC4_.onset, data_hdr.stores.PC6_.onset])';
    end
    
    % Time vector
    ieeg.time = 0 : 1/ieeg.fs : (length(ieeg.data)-1)/ieeg.fs; 
    
    %% Import Analog Channels      
        
    if exist(analog_chn_file, 'file') && ~isempty(Wav5.streams)
        
        analogChns = importAnalog(sprintf('%s/%s/%s_ExtChannels.xlsx', ...
            data_dir, patient_name, patient_name));

        % Define the following fields if you know them. If not just leave them undefined:
        analog.resp = find(strcmp(analogChns.Label, 'respiration')); 
        analog.micro = find(strcmp(analogChns.Label, 'Mic_R')); 
        analog.aud = find(strcmp(analogChns.Label, 'PC_R')); 
        analog.hr = find(strcmp(analogChns.Label, 'heartrate')); 
        % If you named some specific analog channels above. Say here if you wan't the rest deleted.
        analog.deleteChnsWithoutLabel = 1; % Only keep analog channels you gave a number too (default = keep all)

        ana = Wav5.streams.Wav5.data;
        ieeg.analog.fs = Wav5.streams.Wav5.fs;

        ieeg.analog.trial{1} = ana;

        % Make labels
        for i=1:size(ieeg.analog.trial{1},1)
            analog_label{i} = ['AnalogCh' num2str(i)];
        end

        ana2keep=[];
        if isfield(analog,'resp')
            if ~isempty(analog.resp)
                analog_label{analog.resp} = 'resp';
                ana2keep = [ana2keep; analog.resp];
            end
        end

        if isfield(analog,'micro')
            if ~isempty(analog.micro)
                analog_label{analog.micro} = 'micro';
                ana2keep = [ana2keep; analog.micro];
            end
        end

        if isfield(analog,'hr') 
            if ~isempty(analog.hr)
                analog_label{analog.hr} = 'hr';
                ana2keep = [ana2keep; analog.hr];
            end
        end

        if isfield(analog,'aud')
            if ~isempty(analog.aud)
                analog_label{analog.aud} = 'aud';
                ana2keep = [ana2keep; analog.aud];
            end
        end

        % Only keep analog channels that were listed in paramsfile?
        if isfield(analog, 'deleteChnsWithoutLabel')
            if analog.deleteChnsWithoutLabel == 1
                ieeg.analog.label = analog_label(ana2keep);
                ieeg.analog.trial{1} = ana(ana2keep,:);
            else
                ieeg.analog.label = analog_label;
                ieeg.analog.trial{1} = ana;
            end
        end
    end
    
end