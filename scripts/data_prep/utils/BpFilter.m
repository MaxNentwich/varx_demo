% Bandpass filter function

function data_filt = BpFilter(data, order, pass, fs)

% Often collides with fieltrip function
pathCell = regexp(path, pathsep, 'split');

if any(strcmp('/home/max/Documents/MATLAB/Toolboxes/fieldtrip-20180426/external/signal', pathCell))
    rmpath('/home/max/Documents/MATLAB/Toolboxes/fieldtrip-20180426/external/signal')
end

% Transpose data
if size(data,2) > size(data,1)
    data = data';
end

% Filter
[z,p,k] = butter(order, pass/(fs/2), 'bandpass');

[sos,g] = zp2sos(z,p,k);

data_filt = filtfilt(sos, g, double(data));
