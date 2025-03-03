% Lowpass filter function

function data_filt = LpFilter(data, order, pass, fs)

    % Often collides with fieltrip function
    pathCell = regexp(path, pathsep, 'split');

    if any(strcmp('/home/max/Documents/MATLAB/Toolboxes/fieldtrip-20201117/external/signal', pathCell))
        rmpath('/home/max/Documents/MATLAB/Toolboxes/fieldtrip-20201117/external/signal')
    end

    if size(data,2) > size(data,1)
        data = data';
    end

    [z,p,k] = butter(order, pass/(fs/2), 'low');

    [sos,g]  = zp2sos(z,p,k);

    % fvtool(sos);

    data_filt = filtfilt(sos, g, double(data));

end