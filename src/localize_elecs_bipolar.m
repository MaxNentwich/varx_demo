%% Localize the electrodes in cortical lobes

function [loc, fsaverage_coords, atlas] = localize_elecs_bipolar(electrode_names, atlas_select)

    atlas = {'DK_Lobe', 'AparcAseg_Atlas', 'DK_Atlas', 'D_Atlas', 'D_Full'};

    % Consult Noah's table for electrode location
    load('../data/movie_subs_table.mat')

    loc = cell(length(electrode_names),5,2);
    fsaverage_coords = nan(length(electrode_names),3,2);

    for e = 1:length(electrode_names)

        labels_indiv = strsplit(electrode_names{e}, '-');
        
        for l = 1:length(labels_indiv)
            
            sub_elec = strsplit(labels_indiv{l}, '_');   

            if length(sub_elec) == 2
                sub_name = sub_elec(1);
                ch_name = sub_elec(2);
            elseif length(sub_elec) == 3
                sub_name = {[sub_elec{1}, '_', sub_elec{2}]};
                ch_name = sub_elec(3);
            end    
            if strcmp(sub_name, 'NS144_02') && contains(ch_name, 'RGrid')
                ch_name = strrep(ch_name, 'RGrid', 'LGrid');
            end
            if strcmp(sub_name, 'NS148_02') && contains(ch_name, 'RGrid')
                ch_name = strrep(ch_name, 'RGrid', 'RGridHD');
            end
            if strcmp(sub_name, 'NS134') && contains(ch_name, 'RFl')
                ch_name = strrep(ch_name, 'RFl', 'RFL');
            end
            if strcmp(sub_name, 'NS134') && contains(ch_name, 'RIa')
                ch_name = strrep(ch_name, 'RIa', 'Ria');
            end
            if strcmp(sub_name, 'NS138') && contains(ch_name, 'LIa')
                ch_name = strrep(ch_name, 'LIa', 'Lia');
            end

            idx_elec = ismember(movie_subs_table.SubID, sub_name) & ismember(movie_subs_table.Contact, ch_name);

            if sum(idx_elec) ~= 0           
                loc{e,1,l} = movie_subs_table.DK_Lobe{idx_elec};
                loc{e,2,l} = movie_subs_table.AparcAseg_Atlas{idx_elec};
                loc{e,3,l} = movie_subs_table.DK_Atlas{idx_elec};   
                loc{e,4,l} = movie_subs_table.D_Atlas{idx_elec};
                loc{e,5,l} = movie_subs_table.D_Full{idx_elec};

                fsaverage_coords(e,:,l) = movie_subs_table.FSAverage(idx_elec,:);
            else
                loc{e,1,l} = 'Unknown';  
                loc{e,2,l} = 'Unknown'; 
                loc{e,3,l} = 'Unknown'; 
                loc{e,4,l} = 'Unknown'; 
                loc{e,5,l} = 'Unknown'; 
            end
        
        end

    end
    
    if strcmp(atlas_select, 'dk')
        
        loc_dk = loc(:, ismember(atlas, 'DK_Atlas'));
        loc_aparc = loc(:, ismember(atlas, 'AparcAseg_Atlas'));

        idx_unknown_dk = cellfun(@(C) strcmpi(C, 'unknown'), loc_dk);
        idx_unknown_aparc = cellfun(@(C) strcmpi(C, 'unknown'), loc_aparc);

        loc = loc_dk;
        loc(idx_unknown_dk & ~idx_unknown_aparc) = loc_aparc(idx_unknown_dk & ~idx_unknown_aparc);

        loc(cellfun(@(C) contains(C, 'Amygdala'), loc)) = {'amygdala'};
        loc(cellfun(@(C) contains(C, 'Cerebral-White-Matter'), loc)) = {'whitematter'};
        loc(cellfun(@(C) contains(C, 'Inf-Lat-Vent'), loc)) = {'ventricle'};
        loc(cellfun(@(C) contains(C, 'Hippocampus'), loc)) = {'hippocampus'};
        loc(cellfun(@(C) contains(C, 'choroid-plexus'), loc)) = {'choroidplexus'};
        loc(cellfun(@(C) contains(C, 'Brain-Stem'), loc)) = {'brainstem'};
        loc(cellfun(@(C) contains(C, 'Accumbens-area'), loc)) = {'accumbens'};
        loc(cellfun(@(C) contains(C, 'parahippocampal'), loc)) = {'parahippocampal'};
        loc(cellfun(@(C) contains(C, 'insula'), loc)) = {'insula'};
        loc(cellfun(@(C) contains(C, 'Caudate'), loc)) = {'caudate'};
        loc(cellfun(@(C) contains(C, 'Unknown'), loc)) = {'unknown'};

    elseif strcmp(atlas_select, 'lobes')
        
        %% Find electrodes in MTL
        idx_mtl = sum(squeeze(cellfun(@(C) contains(C, {'Amygdala', 'Hippocampus', 'entorhinal', 'parahippocampal'}), loc(:,2,:))), 2) ~= 0;
        loc_mtl = squeeze(loc(:,2,:));
        loc_mtl(~idx_mtl, :) = repmat({'Unknown', 'Unknown'}, sum(~idx_mtl), 1);
        
        loc = squeeze(loc(:,1,:));
        
        loc(~cellfun(@(C) contains(C, 'Unknown', 'IgnoreCase', true), loc_mtl)) = {'MTL'};
       
    elseif strcmp(atlas_select, 'AparcAseg_Atlas')
        
        loc_dk = squeeze(loc(:,3,:));
        
        loc = squeeze(loc(:,2,:));
        
        loc = cellfun(@(C) strrep(C, 'Right-', ''), loc, 'UniformOutput', false);
        loc = cellfun(@(C) strrep(C, 'Left-', ''), loc, 'UniformOutput', false);
        loc = cellfun(@(C) strrep(C, 'ctx-rh-', ''), loc, 'UniformOutput', false);
        loc = cellfun(@(C) strrep(C, 'ctx-lh-', ''), loc, 'UniformOutput', false);
        
        % Use DK_atlas label if the AparcAseg_Atlas label is unknown
        idx_unknown = ismember(loc, 'Unknown') | ismember(loc, 'unknown');
        
        loc(idx_unknown) = loc_dk(idx_unknown);
        
        % Use DK_atlas label if the AparcAseg_Atlas label is White-Matter and AparcAseg_Atlas is not unknown
        idx_wm = cellfun(@(C) contains(C, 'White-Matter'), loc);
        idx_dk_unknown = ismember(loc_dk, 'Unknown') | ismember(loc_dk, 'unknown');
        idx_wm(idx_dk_unknown) = false;
        
        loc(idx_wm) = loc_dk(idx_wm);
        
    end
    
end