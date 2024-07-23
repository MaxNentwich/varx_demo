%% Compute length of responses

function fwhm = resp_length(resp, fs, filt_pass, filt_order)

    % Low-pass filter
    [z,p,k] = butter(filt_order, filt_pass/(fs/2), 'low');
    [sos,g]  = zp2sos(z,p,k);
    
    resp_filt = filtfilt(sos, g, resp')';
    
    % Initialize array
    fwhm = nan(1, size(resp,1));
    
    for ch = 1:size(resp,1)
    
        % Find peaks and get full-width-half-max (FWHM)
        [~, ~, peak_fwhm, peak_prom] = findpeaks(abs(resp_filt(ch,:)),'Annotate','extents');
        
        % Select the most prominent peak
        [~, idx_peak_max] = max(peak_prom);
        peak_fwhm = peak_fwhm(idx_peak_max);
        
        % Save FWHM
        fwhm(ch) = peak_fwhm;
    
    end

end