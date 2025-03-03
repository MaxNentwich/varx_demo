function T = import_label(workbook_file)

    % Read file into table
    T = readtable(workbook_file);

    % Check if all columns are correct
    if sum(ismember(T.Properties.VariableNames, ...
            {'Label', 'TDTChan', 'Out', 'BAD', 'SOZ', 'Spikey', 'Spec', 'FS_label_brain'})) ~= 8

        warning('Column Labels are not correct!\n')

    end

    %% Set NaNs to zero
    T.Out(isnan(T.Out)) = 0;
    T.BAD(isnan(T.BAD)) = 0;
    T.SOZ(isnan(T.SOZ)) = 0;
    T.Spikey(isnan(T.Spikey)) = 0;

    % Make sure channels are sorted by TDT index
    [~,idx] = sort(T.TDTChan); 
    
    if sum(idx' -( 1:length(idx))) ~= 0 
        T = T(idx,:);
        warning('TDT channel labels have been sorted!')
    end
    
end