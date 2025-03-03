
function [saccade_onset, xy, saccade_amplitude, saccade_speed, v, fixation_onset, pos_pre, pos_post, distance_screen] = ...
    detect_saccade_onset(eye, options, visualize)

    % Exclude data outside the screen    
    % Set validity colum to 4 so it gets interpolated later 
    eye.left(13, eye.left(7,:) > 1 | eye.left(7,:) < 0) = 4;
    eye.left(13, eye.left(8,:) > 1 | eye.left(8,:) < 0) = 4;

    eye.right(13, eye.right(7,:) > 1 | eye.right(7,:) < 0) = 4;
    eye.right(13, eye.right(8,:) > 1 | eye.right(8,:) < 0) = 4;

    % Interpolate samples where eye isn't detected (or outside the screen)
    xy_l = interpolateUsingLastGoodValue({eye.left'}, 7:8, 13, [0,1]);
    xy_r = interpolateUsingLastGoodValue({eye.right'}, 7:8, 13, [0,1]);

    % Average left and right eye
    xy = [mean([xy_l(:,1), xy_r(:,1)],2), mean([xy_l(:,2), xy_r(:,2)],2)];

    % x_f = signal.medfilt(xy[:,1], 15)
    xy(eye.left(13,:) > 1 & eye.right(13,:) <= 1) = xy_r(eye.left(13,:) > 1 & eye.right(13,:) <= 1);
    xy(eye.left(13,:) <= 1 & eye.right(13,:) > 1) = xy_l(eye.left(13,:) <= 1 & eye.right(13,:) > 1);

    % Median filter to remove noise
    xy = medfilt1(xy, 20);

    % Velocity (vector sum of x and y direction)
    v = sqrt(([diff(xy(:,1)); 0]*eye.fs).^2 + ([diff(xy(:,2)); 0]*eye.fs).^2);

    %% Detect Saccades (Saccade detection algorithm might need to be refined -> maybe ask Marcin about this some time)
    % Saccades defined as time points with a velocity higher than the standard deviation of the velocity
    saccades = v > options.vel_th * std(v);

    % Dialte and erode -> saccades and overshoot get combined

    % Kernel size for convolution
    k_size = 5;

    saccades = conv(saccades, ones(k_size,1)/k_size, 'same');
    saccades(saccades > 0) = 1;

    saccades = conv(saccades, ones(k_size,1)/k_size, 'same');
    saccades(saccades < 1) = 0;

    % Exclude saccades around excluded samples 
    % Kernel size defining samples before and after blinks and artifacts that should not be included as saccades 
    % ~80 ms

    bad_samples = mean([eye.left(13,:); eye.right(13,:)]);

    % Exlude samples that have one good eye
    bad_samples(eye.left(13,:) > 1 & eye.right(13,:) <= 1) = 0;
    bad_samples(eye.left(13,:) <= 1 & eye.right(13,:) > 1) = 0;

    k_bad = 50;

    bad_samples = conv(bad_samples, ones(k_bad,1)/k_bad, 'same');
    bad_samples(bad_samples > 0) = 1;

    saccades(bad_samples == 1) = 0;

    % Find onset
    saccade_onset = zeros(size(v));
    [~,saccade_idx] = findpeaks(saccades);
    
    % Remove saccades that are too close together 
    saccade_idx(diff(saccade_idx) < options.min_diff*eye.fs) = [];
    
    %% Determine the saccade amplitude and speed
    % this is a relative distance and not converted to pixels or mm
    
    % Cut the window around the saccade -33 to 100 ms (longer after the onset)
    sacc_win = round(eye.fs * options.peri_saccade_window);
    
    % Samples pre and post saccades
    pos_peri_idx = saccade_idx + (-sacc_win(1) : sacc_win(2));
    
    peri_time = (-sacc_win(1) : sacc_win(2))/eye.fs;
    
    % Exlude windows that go outside the edges
    idx_out_of_bounds = sum(pos_peri_idx < 1 | pos_peri_idx > length(xy), 2) ~= 0;
    
    pos_x_peri = nan(size(pos_peri_idx));
    pos_y_peri = nan(size(pos_peri_idx));
    v_peri = nan(size(pos_peri_idx));
    
    x = xy(:,1);
    y = xy(:,2);
    
    pos_x_peri(~idx_out_of_bounds, :) = x(pos_peri_idx(~idx_out_of_bounds, :));
    pos_y_peri(~idx_out_of_bounds, :) = y(pos_peri_idx(~idx_out_of_bounds, :));
    
    v_peri(~idx_out_of_bounds, :) = v(pos_peri_idx(~idx_out_of_bounds, :));
    
    pos_pre = nan(length(saccade_idx), 2);
    pos_post = nan(length(saccade_idx), 2);
    fixation_idx = nan(length(saccade_idx), 1);
    
    for s = 1:length(saccade_idx)
        
        if sum(isnan(v_peri(s,:))) ~= 0
            continue
        end
        
        cutoff = v_peri(s,:) < prctile(v_peri(s,:), 70);
        cutoff = conv(cutoff, ones(5,1), 'same');
        cutoff = cutoff > 1;
        
        cutoff = conv(cutoff, ones(5,1), 'same');
        cutoff = cutoff > 3;
        
        [clusters, n_clusters] = bwlabel(cutoff);
        
        pre_saccade = clusters == 1;
        post_saccade = clusters == n_clusters;
        
        pos_pre(s, 1) = mean(pos_x_peri(s, pre_saccade));
        pos_post(s, 1) = mean(pos_x_peri(s, post_saccade));
        
        pos_pre(s, 2) = mean(pos_y_peri(s, pre_saccade));
        pos_post(s, 2) = mean(pos_y_peri(s, post_saccade));
        
        fixation_idx(s) = pos_peri_idx(s, find(post_saccade, 1, 'first'));
        
%         plot(v_peri(s,:))
%         
%         hold on
%         
%         plot(pre_saccade)
%         plot(post_saccade)
%            
%         hold off
%         pause
        
    end
    
    % Remove empty value at the edges
    idx_nan = isnan(fixation_idx);
    
    % Remove duplicates
    % Correct duplicate frames
    [~, idx_unique] = unique(fixation_idx, 'first');
    idx_dup = 1:length(fixation_idx);
    idx_dup(idx_unique) = [];
    
    idx_dup_fix = false(size(idx_nan));
    idx_dup_fix(idx_dup) = true;
    
    [~, idx_unique] = unique(saccade_idx, 'first');
    idx_dup = 1:length(saccade_idx);
    idx_dup(idx_unique) = [];
    
    idx_dup_sac = false(size(idx_nan));
    idx_dup_sac(idx_dup) = true;
    
    % Find saccades where the fixation is found before the saccade onset
    idx_swap = fixation_idx - saccade_idx < 0;
    
    idx_remove = idx_nan | idx_dup_fix | idx_dup_sac | idx_swap;
    
    saccade_idx(idx_remove) = [];
    pos_pre(idx_remove, :) = [];
    pos_post(idx_remove, :) = [];
    fixation_idx(idx_remove) = [];
    pos_x_peri(idx_remove, :) = []; 
    pos_y_peri(idx_remove, :) = []; 
    pos_peri_idx(idx_remove, :) = []; 
    idx_out_of_bounds(idx_remove) = [];
     
    % Fixation onset
    fixation_onset = zeros(size(v));
    fixation_onset(fixation_idx) = 1;
    
    saccade_amplitude_x = range(pos_x_peri,2);
    saccade_amplitude_y = range(pos_y_peri,2);
    
    %% Convert saccade aplitude to mm  
    % Scale by screen size
    saccade_amplitude_x = saccade_amplitude_x * options.screen_size(1);
    saccade_amplitude_y = saccade_amplitude_y * options.screen_size(2);

    % Convert to mm
    mm_per_pix = mean(options.screen_dimension ./ options.screen_size);
    
    saccade_amplitude_x = saccade_amplitude_x * mm_per_pix;
    saccade_amplitude_y = saccade_amplitude_y * mm_per_pix;

    % Compute amplitude in mm
    saccade_amplitude = sqrt(saccade_amplitude_x.^2 + saccade_amplitude_y.^2);
    
    %% Convert to DVA   
    % Get distance to screen
    eye_pos_left = interpolateUsingLastGoodValue({eye.left'}, 1:3, 13, [0,1]);
    eye_pos_right = interpolateUsingLastGoodValue({eye.right'}, 1:3, 13, [0,1]);
    
    gaze_pos_left = interpolateUsingLastGoodValue({eye.left'}, 9:11, 13, [0,1]);
    gaze_pos_right = interpolateUsingLastGoodValue({eye.right'}, 9:11, 13, [0,1]);
    
    dist_left = mean(sqrt(sum((eye_pos_left - gaze_pos_left).^2, 2)));
    dist_right = mean(sqrt(sum((eye_pos_right - gaze_pos_right).^2, 2)));
    
    distance_screen = mean([dist_left, dist_right]);
    
    % Convert saccade amplitude 
    saccade_amplitude = rad2deg(atan(saccade_amplitude / distance_screen));
    
    %% Determine the saccade speed
    saccade_speed = max(v(pos_peri_idx), [], 2);
    
    %% Now fill the saccade onset vector
    saccade_onset(saccade_idx) = 1;
    
    %% The relationship between saccade amplitude and speed should be linear
    if visualize
        
        % Linear regression
        X = [ones(length(saccade_amplitude),1) saccade_amplitude];
        b = X(~idx_out_of_bounds, :) \ saccade_speed(~idx_out_of_bounds);

        reg_line = X*b;

        figure
        hold on

        plot(saccade_amplitude, saccade_speed, '.')
        plot(saccade_amplitude, reg_line,'--')

        xlabel('Saccade Amplitude')
        ylabel('Saccade Speed')
    
    end
    
end