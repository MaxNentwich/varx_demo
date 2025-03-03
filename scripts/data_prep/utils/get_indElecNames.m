function indElecs = get_indElecNames(allElecs)
% Get unique electrode names from a electrode label array.
% This is used by the functions that makes bipolar montages.
%
% Input is: a cell array with electrode names (such as after import with
% ft_preprocessing. 
% Output is:
% indElecs.labels: unique labels in allElecs
% indElecs.numbers: numbers of the unique labels that were present in data
% 
% (UniqueCell is from Psychtoolbox; match_str from Fieldtrip)
%
% SB: 07/2018

nelec = length(allElecs);

nonum=[];

for i=1:nelec
    
    % extract number in string
    if contains(allElecs{i}, '_')
        B = regexp(allElecs{i},'\_d*','Match');
    else
        B = regexp(allElecs{i},'\d*','Match');   
    end
    %      B=str2double(regexp(allElecs{12},'[\d.]+','match'))
    
    if ~isempty(B)
        
        nonum{i,1} = allElecs{i}(1:strfind(allElecs{i},B{:})-1);
        
        if contains(allElecs{i}, '_')
            nonum{i,2} = str2double(allElecs{i}(strfind(allElecs{i},B{:})+1:end));
        else
            nonum{i,2} = str2double(B);
        end
    else
        nonum{i,1} = allElecs{i};
    end
    
end

% Get unique labels 
[labels,~,~] = UniqueCell(nonum(:,1)); 
indElecs = struct('Labels',labels);

% Get numbers for each unique label in data
for i=1:length(indElecs)
    [sel1, ~] = match_str(nonum(:,1), indElecs(i).Labels);
    indElecs(i).ElecNumbers = cell2mat(nonum(sel1,2))';
end