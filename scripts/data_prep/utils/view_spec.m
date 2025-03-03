function view_spec(data, fs, win, varargin)

    [spec, w, t_spec] = spectrogram(data-mean(data), ...
        round(fs)*win, [], [], fs);
    
    if nargin == 4 && varargin{1} < fs/2
        spec = spec(w < varargin{1}, :);
        w = w(w < varargin{1});
    end
 
    imagesc([min(t_spec) max(t_spec)], [max(w) min(w)], ...
        10*log10(abs(flipud(spec)))); 
    
%     caxis([0 35])
    colorbar
    
    set(gca,'YDir','normal')
%     cb2 = colorbar;
%     colorTitleHandle = get(cb2,'Title');
%     titleString = 'Power/frequency [dB/hz]';
%     set(colorTitleHandle ,'String',titleString);
    xlabel('Time [s]'); ylabel('Frequency [Hz]')
   
end
    