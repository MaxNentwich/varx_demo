function ecog=ecog_bipolarize(ecog_raw,cfg)
% function ecog=ecog_bipolarize(ecog_raw,cfg)
%
% ecog_raw is ecog variable or fieldtrip variable.
%
% If no cfg given, a gui will pop up to give options what channels to bipolarize.
%
% Input:
% cfg.grid_lines - output of derive_grid_lines
%                  eg. grid_lines=derive_grid_lines('RGRd',[8 8]); %MuK
%
% cfg.includeChannels - array of unique stems of electrode labels
%                       eg. if labels are "LG1" "LG2" "LS1" use includeChannels = {'LG' 'LS'};
%                       Tip: You can use indElecs = get_indElecNames(allElecs) to get them from a
%                       label array with electrode names.
%
% Output:
% ecog  - ecog variables in bipolar montage
%
% Dependencies:
% cueEvents is from Penn toolbox eeg_toolbox_v1_3_2 (in IEEG/epiproc/EXTERNAL folder)

%% Check inputs
if nargin<2
    cfg=[];
end

check_cfg(cfg,'ecog_bipolarize.m');
if ~isfield(cfg,'grid_lines'), cfg.grid_lines='no'; end
if ~isfield(cfg,'includeChannels'),   cfg.includeChannels='no'; end

if ~strcmp(cfg.grid_lines,'no') && ~strcmp(cfg.includeChannels,'no')
    error('Can''t define both cfg.grid_lines and cfg.includeChannels');
end

%% Check if data is ecog format and revert data to original reference

if ~isfield(ecog_raw,'ftrip')
    isecog = 0;
    
    fprintf(['No ecog.ftrip field present. ASSUME DATA IS NOT AN ECOG VARIABLE BUT ORIGINAL FIELDTRIP VARIABLE.'...
    'FOR EXAMPLE FROM ft_preprocessing.\n']);
    fprintf('ASSUME DATA IS (AND HAS TO BE IN!!!) ORIGINAL RECORDING REFERENCE.\n');

    orig_input = ecog_raw;
    ecog_raw = [];
    ecog_raw.ftrip=orig_input;
    
    % initiate output variable
    ecog=ecog_raw;
    ecog.edf_labels=ecog_raw.ftrip.label;

else
    isecog = 1;
    
    if isfield(ecog_raw,'ref')
        if isequal(ecog_raw.ref,'bipolar')
            error('ecog_raw variable is ALREADY bipolar.');
        elseif sum(ecog_raw.ref.contributing_chans)
            fprintf('Data are being temporarily reverted to original recording reference.\n');
            orig_reffed_chans=find(ecog_raw.ref.contributing_chans);
            n_trial=length(ecog_raw.ftrip.trial);
            for a=1:n_trial
                ecog_raw.ftrip.trial{a}(orig_reffed_chans,:)=ecog_raw.ftrip.trial{a}(orig_reffed_chans,:)+...
                    repmat(ecog_raw.ref.signal{a},length(orig_reffed_chans),1);
            end
        else
            fprintf('Data were passed to this function with original recording reference.\n');
        end
    else
        fprintf('No ecog.ref field present. ASSUME DATA IS (AND HAS TO BE IN!!!) ORIGINAL RECORDING REFERENCE.\n');
    end
    
    % initiate output variable
    ecog=ecog_raw;
    ecog.edf_labels=ecog_raw.ftrip.label;
    
end
    
n_trial=length(ecog_raw.ftrip.trial);
n_chan=size(ecog_raw.ftrip.trial{1},1);

    
%% Check which method was used to choose channels to bipolarize

% get unique label stems
indlabel = get_indElecNames(ecog_raw.ftrip. label);

% % % No electrodes to be bipolarized are given as input. Open a gui to select.
% % if any(strcmp(cfg.grid_lines,'no')) && any(strcmp(cfg.includeChannels,'no'))
% %     
% %     % make cell array
% %     allLabels=getStructField(indlabel,'Labels');
% %     % choose elecs to use
% %     elec2use = electrode_select(allLabels,'Click Electrodes to include.');
% %     
% % % Input is selected electrodes that should be bipolarized
% % elseif ~strcmp(cfg.includeChannels,'no')    
% %     elec2use = cfg.includeChannels;
% %    
% % % Gridlines is input. Just bipolarize all there is. Treat electrodes in 
% % % cfg.gridlines as grid, the rest as strips/depth    
% % elseif ~strcmp(cfg.grid_lines,'no')
% %     elec2use = 'all';
% % else
    elec2use = 'all';
% % end

%% Get bipolar neighbors
% 
% if strcmp(elec2use,'all') % this is not even correct output type
%     bip_nbors=get_bipolar_nbors(ecog_raw,grid_lines);
% else
    grid_lines={};
    %define dimensions of strips / grids / depths
    for i=1 : length(cfg.includeChannels)
        
        % get number of electrodes
        cueEvents = cfg.includeChannels(i);%filterStruct(indlabel,['strcmp(Labels,'''  elec2use{i} ''')']);
        
        checkthanumbers = 1;
        while checkthanumbers == 1
%             prompt = {'Set dimensions (eg. 8 8 or 1 6). Preset is min and max. # of rows (1 if depth or strip)'; '# of columns'};
            dlg_title = ['Electrode ' cueEvents.Labels];
            defaultans = {num2str(min(cueEvents.ElecNumbers)), num2str(max(cueEvents.ElecNumbers))};
            num_lines=1;
%             answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
%             nrow = str2num(answer{1});
%             ncol = str2num(answer{2});

            % This doesn't make sense -> changed on 9/17/2019
            % nrow = str2num(defaultans{1});
            nrow = 1;
            ncol = str2num(defaultans{2});
            if nrow*ncol > max(cueEvents.ElecNumbers)
               fprintf('HEY!! The number of rows*col can''t be higher than the highest electrode number in the data...try again');
            else
                checkthanumbers=0;
            end                
        end
        
        % Collect across electrodes
        % (theoretically grid_lines is needed as input to get_bipolar_nbors only if there are
        % grids. Theoretically not needed for strips and depth. we do it for all electrodes here, just
        % because it's easier coding wise)
        grid_lines_temp={};
        grid_lines_temp=derive_grid_lines(cueEvents.Labels,[nrow ncol]);
        grid_lines = [grid_lines; grid_lines_temp];
    end
    
    bip_nbors=get_bipolar_nbors(ecog_raw,grid_lines,elec2use);
% end

n_pairs=size(bip_nbors,1);
ecog.ftrip.label=cell(n_pairs,1);


%% loop through pairs
n_trials=length(ecog_raw.ftrip.trial);
for a=1:n_trials
    n_tpts=length(ecog_raw.ftrip.time);
    ecog.ftrip.trial{a}=zeros(n_pairs,n_tpts);
end

used_pairs=zeros(n_pairs,1);
for a=1:n_pairs
    id1=findstr_in_cell(bip_nbors{a,1},ecog_raw.ftrip.label);
    id2=findstr_in_cell(bip_nbors{a,2},ecog_raw.ftrip.label);
    if isempty(id1)
        warning('Could not find electrode %s in ecog_raw var.',bip_nbors{a,1});
    else
        if isempty(id2)
            warning('Could not find electrode %s in ecog_raw var.',bip_nbors{a,1});
        else
            used_pairs(a)=1;
            for b=1:n_trials
                ecog.ftrip.trial{b}(a,:)=ecog_raw.ftrip.trial{b}(id1,:)-ecog_raw.ftrip.trial{b}(id2,:);
            end
            ecog.ftrip.label{a}=[bip_nbors{a,1} '-' bip_nbors{a,2}];
        end
    end
end


%% Remove missed elecs
missed_ids=find(used_pairs==0);
fprintf('%d possible pairs missed.\n',length(missed_ids));
ecog.ftrip.label=ecog.ftrip.label(find(used_pairs));
for a=1:n_trials
    ecog.ftrip.trial{a}=ecog.ftrip.trial{a}(find(used_pairs),:);
end

%% erase old data and affected fields; 
%  If not ecog variable make sure output stays in fieldtrip format
if isecog
    ecog.ref='bipolar';
    ecog.psd=[];
else
    temp=ecog;
    ecog = ecog.ftrip;
end
