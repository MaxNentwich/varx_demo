    
function pupil = process_pupil(eye, options)

    % Find next odd sample 
    s_edge = ceil(options.t_edge_pupil * eye.fs) / 2 * 2 + 1;
    
    pupil_l = eye.left(12,:);
    pupil_r = eye.right(12,:);
    
    bad_l = eye.left(13,:) > 1;
    bad_r = eye.right(13,:) > 1;
    
    % Extend around the bad segments
    bad_l = conv(bad_l, ones(s_edge,1), 'same') > 0;
    bad_r = conv(bad_r, ones(s_edge,1), 'same') > 0;
    
    pupil_samples = 1:length(pupil_l);
    
    if sum(bad_l) / length(bad_l) < 0.5
        pupil_l(bad_l) = interp1(pupil_samples(~bad_l), pupil_l(~bad_l), pupil_samples(bad_l), 'linear', 'extrap');
    end

    if sum(bad_r) / length(bad_r) < 0.5
        pupil_r(bad_r) = interp1(pupil_samples(~bad_r), pupil_r(~bad_r), pupil_samples(bad_r), 'linear', 'extrap');
    end

    % Linear interpolation for edges
    if sum(bad_l) / length(bad_l) < 0.5
            
        clust_l = bwlabel(bad_l);
        
        id_clust = unique(clust_l);
        id_clust = id_clust(id_clust ~= 0);
        
        clust_l(ismember(clust_l, setdiff(id_clust, [clust_l(1), clust_l(end)]))) = 0;
        
        pupil_l(clust_l~=0) = interp1(pupil_samples(clust_l==0), pupil_l(clust_l==0), pupil_samples(clust_l~=0), 'nearest', 'extrap');
        
    end

    % Right
    if sum(bad_r) / length(bad_r) < 0.5
    
        clust_r = bwlabel(bad_r);
    
        id_clust = unique(clust_r);
        id_clust = id_clust(id_clust ~= 0);
        
        clust_r(ismember(clust_r, setdiff(id_clust, [clust_r(1), clust_r(end)]))) = 0;
        
        pupil_r(clust_r~=0) = interp1(pupil_samples(clust_r==0), pupil_r(clust_r==0), pupil_samples(clust_r~=0), 'nearest', 'extrap');
        
    end

    % Average
    if sum(bad_l) / length(bad_l) < 0.5 && sum(bad_r) / length(bad_r) >= 0.5
        pupil = pupil_l;
    elseif sum(bad_r) / length(bad_r) < 0.5 && sum(bad_l) / length(bad_l) >= 0.5
        pupil = pupil_r;
    elseif sum(bad_r) / length(bad_r) < 0.5 && sum(bad_l) / length(bad_l) < 0.5
        pupil = mean([pupil_l; pupil_r]);
    end

    if sum(bad_r) / length(bad_r) >= 0.5 && sum(bad_l) / length(bad_l) >= 0.5
        pupil = nan(1, length(pupil_samples));
    
    else

        % Interpolate NaNs
        pupil(isnan(pupil)) = interp1(pupil_samples(~isnan(pupil)), pupil_l(~isnan(pupil)), pupil_samples(isnan(pupil)), 'nearest', 'extrap');
    
        % Smooth
        pupil = LpFilter(pupil, 5, options.pupil_lp, eye.fs);

    end

end