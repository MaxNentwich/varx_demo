function bip_nbors=get_bipolar_nbors(ecog,grid_lines,elec2use)
%function bip_nbors=get_bipolar_nbors(ecog,grid_lines)
%
% Input:
%  ecog       - Mehtalab ecog variable (Mandatory)
%  grid_lines - Use if there is a grid. Output of derive_grid_lines.m. A 2D cell array. 
%               The first column is the grid_stem, the second is a vector of
%               neighbors in a row or column in the grid.
%  elec2use   - Stem of electrodes to use. 
%               eg. if labels are "LG1" "LG2" "LS1" use {'LG' 'LS'} 
%               some helpers to get stems for label array: 
%               indlabel = get_indElecNames(ecog_raw.ftrip. label)
%               elec2use=getStructField(indlabel,'Labels');
%
% Output:
%  bip_nbors - 2D cell array of electrode pairs for doing bipolar
%              referencing. Each row is a different pair.
%
% Examples:
% 1) you want to include all labels and there's no grid:
% bip_nbors=get_bipolar_nbors(ecog);
%
% 2) include all labels, there's a grid
% grid_lines = derive_grid_lines(ecog,'GridLabelStem');
% bip_nbors=get_bipolar_nbors(ecog,grid_lines); 
%
% 3) include selected, no grid 
%    bip_nbors=get_bipolar_nbors(ecog,[],{'LD' 'LS'}); 
%
%
% Author:
% David Groppe
% 2/2014
% adapted SB: 8/2018
%
% Works if electrode coordinates not part of ecog var, but you don't get
% figs.


% make work for cut grids?


%%

if nargin>1
    n_grid_lines=size(grid_lines,1);
    for a=1:n_grid_lines
        grid_stems{a}=grid_lines{a,1};
    end
    grid_stems=unique(grid_stems);
else
    grid_stems=[];
    n_grid_lines=[];
    grid_lines=[];
end

onlyusethese=[];
if nargin == 3
    onlyusethese = elec2use;
end

%% first figure out what electrodes we have
indlabel = get_indElecNames(ecog.ftrip.label);
% make cell array
all_stems = getStructField(indlabel,'Labels');

% Get min and max of each stem
stem_min=zeros(length(indlabel),1);
stem_max=zeros(length(indlabel),1);
for i=1:length(indlabel)
    if ~isempty(indlabel(i).ElecNumbers)
        stem_min(i) = min(indlabel(i).ElecNumbers);
        stem_max(i) = max(indlabel(i).ElecNumbers);
    end
end

% Report to command line
for a=1:length(all_stems)
    fprintf('%s: %d to %d\n',all_stems{a},stem_min(a),stem_max(a));
end

if ~isempty(onlyusethese)
    [sel1,~] = match_str(all_stems,onlyusethese);
    stem_min = stem_min(sel1);
    stem_max = stem_max(sel1);
    all_stems = all_stems(sel1);
end
n_stem = length(all_stems);

disp(['You chose to use only '  all_stems]);

%% Now go through and pair neighbors
bip_nbors=[];
bip_nbors_dist=[]; %distance check
n_pair=0;
skipped_notice=0;
% do strips first
for a=1:n_stem
    for b=stem_min(a):(stem_max(a)-1)
        if ~ismemberi(all_stems{a},grid_stems)
            temp1=[all_stems{a} num2str(b)];
            temp2=[all_stems{a} num2str(b+1)];
            id1=findstr_in_cell(temp1,ecog.ftrip.label);
            id2=findstr_in_cell(temp2,ecog.ftrip.label);
            if isempty(id1) && isempty(id2)
                fprintf('Could not find channels %s and %s. Skipping that pair.\n',temp1,temp2);
                skipped_notice=1;
            elseif isempty(id1)
                fprintf('Could not find channels %s. Skipping its possible pairs.\n',temp1);
                skipped_notice=1;
            elseif isempty(id2)
                fprintf('Could not find channels %s. Skipping its possible pairs.\n',temp2);
                skipped_notice=1;
            else
                n_pair=n_pair+1;
                bip_nbors{n_pair,1}=temp1;
                bip_nbors{n_pair,2}=temp2;
                if isfield(ecog,'presnap_coorRAS')
                    bip_nbors_dist(n_pair)=sqrt(sum( (ecog.presnap_coorRAS(id1,:)-ecog.presnap_coorRAS(id2,:)).^2) );
                else
                    bip_nbors_dist(n_pair)=NaN;
                end
            end
        end
    end
end

%now do grid(s)
if ~isempty(grid_lines)
    for a=1:n_grid_lines
        for b=1:(length(grid_lines{a,2})-1)
            temp1=[grid_lines{a,1} num2str(grid_lines{a,2}(b))];
            temp2=[grid_lines{a,1} num2str(grid_lines{a,2}(b+1))];
            id1=findstr_in_cell(temp1,ecog.ftrip.label);
            id2=findstr_in_cell(temp2,ecog.ftrip.label);
            if isempty(id1) && isempty(id2)
                fprintf('Could not find channels %s and %s. Skipping that pair.\n',temp1,temp2);
                skipped_notice=1;
            elseif isempty(id1)
                fprintf('Could not find channels %s. Skipping its possible pairs.\n',temp1);
                skipped_notice=1;
            elseif isempty(id2)
                fprintf('Could not find channels %s. Skipping its possible pairs.\n',temp2);
                skipped_notice=1;
            else
                n_pair=n_pair+1;
                bip_nbors{n_pair,1}=temp1;
                bip_nbors{n_pair,2}=temp2;
                if isfield(ecog,'presnap_coorRAS')
                    bip_nbors_dist(n_pair)=sqrt(sum( (ecog.presnap_coorRAS(id1,:)-ecog.presnap_coorRAS(id2,:)).^2) );
                else
                    bip_nbors_dist(n_pair)=NaN;
                end
            end
        end
    end
end

if  skipped_notice==1
    fprintf('NOTE! SOME ELECTRODES WERE SKIPPED BECAUSE THE PAIRS WERE NOT FOUND.');
    fprintf('Maybe because non-contiguous numbers are not bipolarized (eg LS1, LS3, LS5).')
    fprintf('Unless it''s a high density grid, bipolarizing them likely does not make sense.');
end

%%
if isfield(ecog,'presnap_coorRAS')
    figure; clf; hold on;
    for a=1:n_pair
        h=plot(a,bip_nbors_dist(a),'*');
        click_text(h,[bip_nbors{a,1} '-' bip_nbors{a,2}]);
    end
    xlabel('Pair');
    ylabel('Euc Dist (mm)');
    title('Click asterisk to get corresponding pair name');
end

%%
if isfield(ecog,'presnap_coorRAS')
    figure; clf;
    for a=1:n_pair
        id1=findstr_in_cell(bip_nbors{a,1},ecog.ftrip.label,1);
        id2=findstr_in_cell(bip_nbors{a,2},ecog.ftrip.label,1);
        h=plot3([ecog.presnap_coorRAS(id1,1) ecog.presnap_coorRAS(id2,1)], ...
            [ecog.presnap_coorRAS(id1,2) ecog.presnap_coorRAS(id2,2)], ...
            [ecog.presnap_coorRAS(id1,3) ecog.presnap_coorRAS(id2,3)],'o-');
        hold on;
        click_text(h,[bip_nbors{a,1} '-' bip_nbors{a,2}]);
        click_text3d(h,[bip_nbors{a,1} '-' bip_nbors{a,2}],3);
    end
    % xlabel('Pair');
    % ylabel('Euc Dist (mm)');
    title('Click bar to get corresponding pair name');
    axis square;
    
    v=axis;
    %ax_len=min(abs(v));
    plot3([0 0],[0 0],[0 v(6)],'k-');
    text(0,0,v(6)*1.1,'S');
    plot3([0 0],[0 v(4)],[0 0],'k-');
    text(0,v(4)*1.1,0,'A');
    if abs(v(2))<abs(v(1))
        plot3([0 v(1)],[0 0],[0 0],'k-');
        text(v(1)*1.1,0,0,'L');
    else
        plot3([0 v(2)],[0 0],[0 0],'k-');
        text(v(2)*1.1,0,0,'R');
    end
end
