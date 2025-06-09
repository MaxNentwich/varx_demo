%% Compute all VARX models

%% Define paths and parameters
data_dir = '/media/max/Workspace/Data/movies_legacy_sample/';
et_dir = '/media/max/12T_HDD/ECoGData/Tobii/Patients';
cut_dir = '../data/scene_cuts_frames';
contrast_dir = '../data/ContrastFull';
motion_dir = '../data/Optic_flow';
speech_dir = '../data/speech_files';
sample_data_dir = '/media/max/Workspace/Data/varx_data';

addpath(genpath('../src'))

% models
model_dir = '../results/models_revision_1';
if exist(model_dir, 'dir') == 0
    mkdir(model_dir)
end

if exist(sample_data_dir, 'dir') == 0
    mkdir(sample_data_dir)
end

compute_models = true;

% Model parameters
tl=0.6; % length of the mTRF model filter in seconds
ta=0.1; % length of the VARX model filters in seconds
ta_hfa=0.1;
tb=0.6; % length of the FIR model filter in seconds
gamma = 0.3; 

stim_features = {'fixations', 'fixation_novelty', 'film_cuts', ...
    'motion', 'audio_env', 'acoustic_edges'};

features_sort = {'film_cuts', 'audio_env', 'acoustic_edges', ...
    'fixation_novelty', 'fixations', 'motion'};

videos = {...
    'Despicable_Me_English.mat', ...
    'Despicable_Me_Hungarian.mat', ...
    'Monkey1_Rep_1.mat', ...
    'Monkey2_Rep_1.mat', ...
    'Monkey5_Rep_1.mat', ...
    'The_Present_Rep_1.mat',...
    'The_Present_Rep_2.mat', ...
    'Inscapes.mat',...
    'Resting_fixation.mat', ...
    'Eyes_Closed_Rest.mat', ...
    'Eyes_Closed_Rest_2.mat'};

% List of videos for individual analysis
vids_select = {'Despicable_Me_English', 'Resting_fixation', 'Inscapes', 'Monkey1_Rep_1', 'Eyes_Closed_Rest'};

% Saccade novelty
saccade_features = readtable('../data/saccade_features_and_distance.csv', Delimiter=',');

% Parts of filenames to find frame indices
fr_str = 'fr_';
file_str = '.jpg';

%%
patients = dir(data_dir);
patients = patients(3:end);

%% Compute all Models
for pat = 1:length(patients)

    out_file = sprintf('%s/%s_varx_models.mat', model_dir, patients(pat).name);

    % if compute_models && exist(out_file, 'file') ~= 0
    %     continue
    % end

    data_pat_dir = sprintf('%s/%s', sample_data_dir, patients(pat).name);
    
    if exist(data_pat_dir, 'dir') == 0
        mkdir(data_pat_dir)
    end

    X = cell(length(videos),1);
    Y = cell(length(videos),1);
    HFA = cell(length(videos), 1);
    vid_recs = cell(length(videos),1);

    for vid=1:length(videos)
        
        video = videos{vid};
        [~,vid_name_legacy] = fileparts(video);
    
        vid_file = [data_dir patients(pat).name '/LFP/' video];
        if exist(vid_file, 'file') == 0, continue, end
    
        label_dir = sprintf('%s/%s/matlab_data/Envelope_phase/BHA', et_dir, patients(pat).name);
    
        load([data_dir patients(pat).name '/LFP/' video],'lfp','fs');
        load([data_dir patients(pat).name '/BHA/' video],'envelope');
        load([data_dir patients(pat).name '/Stimuli/' video], 'saccades', 'fixations', 'frame_saccade', 'frame_fixation');
        load(sprintf('%s/%s', label_dir, video), 'labels')

        if strcmp(patients(pat).name, 'NS155') && strcmp(video, 'Resting_fixation.mat')
            idx_extra = find(ismember(labels, 'ROs15-ROs16'));
            lfp(:,idx_extra) = [];
            envelope(:,idx_extra) = [];
            labels(idx_extra) = [];
        end

        time = (0:length(lfp)-1) / fs;

        if strcmp(patients(pat).name, 'NS138') && strcmp(video, 'Eyes_Closed_Rest_2.mat')
            idx_rest = time > 2*fs & time <= 5*fs;
            lfp = lfp(idx_rest, :);
            envelope = envelope(idx_rest, :);   
            time = time(idx_rest);
            time = time - time(1);
        end

        if strcmp(video, 'Eyes_Closed_Rest.mat') || strcmp(video, 'Eyes_Closed_Rest_2.mat')
            load([data_dir patients(pat).name '/Stimuli/' 'Despicable_Me_English.mat'], 'saccades', 'fixations', 'frame_saccade', 'frame_fixation');
            saccades = saccades(1:length(time));
            fixations = fixations(1:length(time));
        end

        %% Remove fixations and saccade outliers
        if length(time) > length(saccades)
            time = time(1:length(saccades));
        end

        %% Load fixation novelty 
        idx_pat = cellfun(@(C) contains(C, patients(pat).name), saccade_features.pre_saccade_file);
        saccade_features_pat = saccade_features(idx_pat, :);

        vid_name_nov = strrep(vid_name_legacy, '_Rep_1', '');
        vid_name_nov = strrep(vid_name_nov, '_Rep_2', '');

        idx_vid = cellfun(@(C) contains(C, sprintf('%s_%s', patients(pat).name, vid_name_nov)), saccade_features_pat.pre_saccade_file);

        if sum(idx_vid) ~= 0

            saccade_features_vid = saccade_features_pat(idx_vid, :);

            frame_pre = cellfun(@(C) str2double(C(regexp(C, fr_str)+length(fr_str):regexp(C, file_str)-1)), ...
                saccade_features_vid.pre_saccade_file);
            frame_post = cellfun(@(C) str2double(C(regexp(C, fr_str)+length(fr_str):regexp(C, file_str)-1)), ...
                saccade_features_vid.post_saccade_file);

            if strcmp(vid_name_nov, 'The_Present')

                idx_split = find(diff(frame_post) < 0);

                if contains(vid_name_legacy, 'Rep_1')

                    saccade_features_vid = saccade_features_vid(1:idx_split, :);

                    frame_pre = frame_pre(1:idx_split);
                    frame_post = frame_post(1:idx_split);

                elseif contains(vid_name_legacy, 'Rep_2')

                    saccade_features_vid = saccade_features_vid(idx_split+1:end, :);

                    frame_pre = frame_pre(idx_split+1:end);
                    frame_post = frame_post(idx_split+1:end);

                end

            end
    
            if contains(video, 'Monkey')
                frame_saccade = round(frame_saccade/2);
                frame_fixation = round(frame_fixation/2);
            end

            saccade_features_vid = saccade_features_vid(...
            ismember(frame_pre, frame_saccade) & ismember(frame_post, frame_fixation), :);

            novelty = saccade_features_vid.distance;
    
            % Keep only saccades that match
            frame_pre = cellfun(@(C) str2double(C(regexp(C, fr_str)+length(fr_str):regexp(C, file_str)-1)), ...
                saccade_features_vid.pre_saccade_file);
            frame_post = cellfun(@(C) str2double(C(regexp(C, fr_str)+length(fr_str):regexp(C, file_str)-1)), ...
                saccade_features_vid.post_saccade_file);
    
            idx_saccade = ismember(frame_saccade, frame_pre);
            idx_fixation = ismember(frame_fixation, frame_post);
    
            sample_saccade = find(saccades);
            sample_fixation = find(fixations);
    
            saccades(sample_saccade(~idx_saccade)) = 0;
            fixations(sample_fixation(~idx_fixation)) = 0;

            fixation_novelty = fixations;
    
            fixation_novelty(fixation_novelty == 1) = novelty;

        else

            % For recordings without novelty (inscapes, rest) fill in
            % uncorrelated features
            fixation_novelty = fixations;
            fixation_novelty(fixation_novelty == 1) = saccade_features(1:sum(fixations==1), :).distance;

        end

        %% Add saccade vector to stimuli

        x = fixations;
        y = lfp; 
        fs_neural = fs;
        hfa = envelope;
    
        if length(x) < length(y)
            y = y(1:length(x), :);
            lfp = lfp(1:length(x), :);
            hfa = hfa(1:length(x), :);
        end

        if length(hfa) < length(x)
            hfa = [hfa; repmat(hfa(end,:), length(x)-length(hfa), 1)];
        end

        if size(x,2) > size(x,1), x = x'; end

        %% Add fixation onset and novelty
        x = [x, fixation_novelty];

        %% Load film cuts
        if strcmp(vid_name_legacy, 'Inscapes') || strcmp(vid_name_legacy, 'Resting_fixation') ...
                || strcmp(vid_name_legacy, 'Eyes_Closed_Rest')  || strcmp(vid_name_legacy, 'Eyes_Closed_Rest_2') 
            vid_name_scenes = 'Despicable_Me_English';
            vid_name_eye = 'Despicable_Me_English.mat';
        else
            vid_name_scenes = strrep(vid_name_legacy, '_Rep_1', '');
            vid_name_scenes = strrep(vid_name_scenes, '_Rep_2', '');
            vid_name_eye = video;
        end
   
        n_samples = length(time);
    
        load(sprintf('%s/%s/matlab_data/Eyetracking/%s', et_dir, patients(pat).name, vid_name_eye), 'eye')      
        scene_table = readtable(sprintf('%s/%s_scenes.xlsx', cut_dir, vid_name_scenes));
        
        if strfind(vid_name_scenes, 'Monkey') == 1
            scene_frame = scene_table.Frame * 2;
            scene_frame = scene_frame(scene_frame < length(eye.frame_time));
        else
            scene_frame = scene_table.Frame;
        end

        scenes_time = eye.frame_time(scene_frame);
        scenes_samp = round(interp1(time-time(1), 1:n_samples, scenes_time));

        scenes_samp(isnan(scenes_samp)) = [];

        scenes_vec = zeros(1, n_samples);
        scenes_vec(scenes_samp) = 1;
    
        if size(scenes_vec,1) ~= size(x,1)
            scenes_vec = scenes_vec';
        end

        x = [x, scenes_vec];
    
        %% Motion (optical flow)
        motion_vid = sprintf('%s/%s.mat', motion_dir, vid_name_scenes);

        if exist(motion_vid, 'file') == 0
            motion_vid = sprintf('%s/Despicable_Me_English.mat', motion_dir);
        end

        load(motion_vid, 'optic_flow', 'fr');
        optic_flow = optic_flow(1:length(eye.frame_time));
        frame_time = eye.frame_time;

        idx_dup = diff(eye.frame_time) == 0;
        frame_time(idx_dup) = [];
        optic_flow(idx_dup) = [];

        motion = interp1(frame_time, optic_flow, time, 'linear', 'extrap')';

        x = [x, motion];

        clearvars eye

        %% Load audio
        if strcmp(vid_name_legacy, 'Inscapes')
            speech_vid = sprintf('%s/%s.mat', speech_dir, vid_name_legacy);
        else
            speech_vid = sprintf('%s/%s.mat', speech_dir, vid_name_scenes);
        end

        if exist(speech_vid, 'file') == 0
            speech_vid = sprintf('%s/Despicable_Me_English.mat', speech_dir);
        end
    
        load(speech_vid, 'audio', 'fs')
        time_audio = (0:size(audio, 1)-1) / fs;

        audio_env = abs(hilbert(audio));

        audio_env = audio_env(time_audio < time(end));
        time_audio = time_audio(time_audio < time(end));

        audio_ds = resample(audio_env, fs_neural, fs);

        if length(audio_ds) < length(x)
            audio_ds = [audio_ds, zeros(1, length(x) - length(audio_ds))];
        elseif length(audio_ds) > length(x)
            audio_ds = audio_ds(1:length(audio_ds)-length(x));
        end
        
        if size(audio_ds,1) ~= size(x,1)
            audio_ds = audio_ds';
        end

        x = [x, audio_ds];

        %% Find acoustic edges
        [~, idx_edge] = findpeaks([0; diff(audio_ds)], 'MinPeakProminence', 0.025);
    
        edge_vec = zeros(size(audio_ds));
        edge_vec(idx_edge) = 1;

        x = [x, edge_vec];

        %%
        % subtract mean to avoid dealing with it in the regression
        x = x-mean(x); y=y-mean(y); hfa = hfa-mean(hfa);
    
        % normalize scale for output so that penalty affects them equally
        y = y./std(y);
        x = x./std(x);
        hfa = hfa./std(hfa);
    
        X{vid}=x; Y{vid}=y; HFA{vid} = hfa;

        vid_recs{vid} = vid_name_legacy;
            
        % Save data
        fs = fs_neural;

        if strcmp(vid_name_legacy, 'Resting_fixation')
            lfp = lfp(1:fs*300, :);
            hfa = hfa(1:fs*300, :);
            fixations = fixations(1:fs*300);
            fixation_novelty = fixation_novelty(1:fs*300);
            save(sprintf('%s/%s.mat', data_pat_dir, vid_name_legacy), ...
                'lfp', 'hfa', 'fixations', 'fixation_novelty', 'fs');
        elseif strcmp(vid_name_legacy, 'Eyes_Closed_Rest')
            
            if size(lfp,1) > fs*300
                nd = round((size(lfp,1) - fs*300) / 2);
                lfp = lfp(nd:nd-end,:);
                hfa = hfa(nd:nd-end,:);
                fixations = fixations(nd:nd-end,:);
                fixation_novelty = fixation_novelty(nd:nd-end,:);
            end

            save(sprintf('%s/%s.mat', data_pat_dir, vid_name_legacy), ...
                'lfp', 'hfa', 'fixations', 'fixation_novelty', 'fs');

        elseif strcmp(vid_name_legacy, 'Inscapes')
            audio_env = audio_ds';
            save(sprintf('%s/%s.mat', data_pat_dir, vid_name_legacy), ...
                'lfp', 'hfa', 'audio_env', 'edge_vec', 'fixations', ...
                'fixation_novelty', 'fs');
        else 
            audio_env = audio_ds';
            film_cuts = scenes_vec';
            save(sprintf('%s/%s.mat', data_pat_dir, vid_name_legacy), ...
                'lfp', 'hfa', 'film_cuts', 'motion', 'audio_env', ...
                'edge_vec', 'fixations', 'fixation_novelty', 'fs');

            if strcmp('Despicable_Me_English', vid_name_legacy)
                lfp = lfp(1:fs*300, :);
                hfa = hfa(1:fs*300, :);
                audio_env = audio_env(1:fs*300);
                edge_vec = edge_vec(1:fs*300);
                film_cuts = film_cuts(1:fs*300);
                motion = motion(1:fs*300);
                fixations = fixations(1:fs*300);
                fixation_novelty = fixation_novelty(1:fs*300);
                save(sprintf('%s/Despicable_Me_English_5min.mat', data_pat_dir), ...
                    'lfp', 'hfa', 'film_cuts', 'motion',  ...
                    'audio_env', 'edge_vec', 'fixations',  ...
                    'fixation_novelty', 'fs');
            end
        end

    end
    
    if compute_models 

        % Remove empty entries
        idx_empty = cellfun(@(C) isempty(C), X);
        X(idx_empty) = [];
        Y(idx_empty) = [];
        HFA(idx_empty) = [];
        vid_recs(idx_empty) = [];
     
        %% Compute VARX model
        na=ta*fs_neural;
        na_hfa = ta_hfa * fs_neural;
        nb=tb*fs_neural; 
        L=tl*fs_neural;

        idx_comb = ~ismember(vid_recs, {'Eyes_Closed_Rest', 'Eyes_Closed_Rest_2'});
        
        % base = basis(n_base,nb,'normal');
        m_varx = varx(Y(idx_comb), na, X(idx_comb), nb, gamma);
        m_varx_hfa = varx(HFA(idx_comb), na_hfa, X(idx_comb), nb, gamma);

        % derive the corresponding impulse responses
        H = varx_trf(m_varx.B, m_varx.A, L);
        H_hfa = varx_trf(m_varx_hfa.B, m_varx_hfa.A, L);

        % now lets do just an AR model
        m_var = varx(Y(idx_comb),na,[],[],gamma);
        m_var_hfa = varx(HFA(idx_comb),na_hfa,[],[],gamma);
        
        %% Remove one feature at a time (permute it)
        m_varx_features = cell(length(stim_features),1);
        m_varx_features_hfa = cell(length(stim_features),1);

        for feat = 1:length(stim_features)

            X_shift = cell(1,length(X));

            for ix = 1:length(X)
                idx_shift = round(length(X{ix})/2);
                X_shift{ix} = X{ix};
                X_shift{ix}(:,feat) = [X_shift{ix}(idx_shift+1:end, feat); X_shift{ix}(1:idx_shift, feat)];
            end

            m_varx_features{feat} = varx(Y(idx_comb),na,X_shift(idx_comb),nb,gamma);
            m_varx_features_hfa{feat} = varx(HFA(idx_comb),na_hfa,X_shift(idx_comb),nb,gamma);

        end

        %% Add one feature at a time
        m_varx_additive = cell(length(stim_features)+1, 1);
        m_varx_additive_hfa = cell(length(stim_features)+1, 1);

        for a = 1:length(stim_features)+1

            X_shift = cell(1,length(X));

            col_shift = ~ismember(stim_features, features_sort(1:a-1));

            for ix = 1:length(X)
                idx_shift = round(length(X{ix})/2);
                X_shift{ix} = X{ix};
                X_shift{ix}(:,col_shift) = [X_shift{ix}(idx_shift+1:end, col_shift); X_shift{ix}(1:idx_shift, col_shift)];
            end

            m_varx_additive{a} = varx(Y(idx_comb),na,X_shift(idx_comb),nb,gamma);
            m_varx_additive_hfa{a} = varx(HFA(idx_comb),na_hfa,X_shift(idx_comb),nb,gamma);

        end

        %% Example of correlated feature
        example_features = {{'audio_env'}, {'audio_env', 'acoustic_edges'}, {'audio_env', 'acoustic_edges', 'fixations'}};
        m_varx_example = cell(length(example_features), 1);

        for ex = 1:length(example_features)

            X_shift = cell(1,length(X));

            col_shift = ~ismember(stim_features, example_features{ex});

            for ix = 1:length(X)
                idx_shift = round(length(X{ix})/2);
                X_shift{ix} = X{ix};
                X_shift{ix}(:,col_shift) = [X_shift{ix}(idx_shift+1:end, col_shift); X_shift{ix}(1:idx_shift, col_shift)];
            end

            m_varx_example{ex} = varx(Y(idx_comb),na,X_shift(idx_comb),nb,gamma);

        end

        %% Example of removing some channels 

        % Remove half of the channels that don't respond to the auditory
        % envelope
        idx_non_sig = m_varx.B_pval(:, strcmp(stim_features, 'audio_env')) > 0.05;
        ch_id = 1:length(idx_non_sig);

        ch_id_non_sig = ch_id(idx_non_sig);
        ch_remove = ch_id_non_sig(1:round(length(idx_non_sig)/2));

        ch_keep = setdiff(ch_id, ch_remove);

        Y_cut = cell(1,length(Y));

        for iy = 1:length(X)
            Y_cut{iy} = Y{iy}(:, ch_keep);
        end

        m_varx_ch_cut = varx(Y_cut(idx_comb), na, X(idx_comb), nb, gamma);

        %% Compute models for each movie/condition separately
        m_varx_mov = cell(length(vids_select),1);
        m_varx_mov_hfa = cell(length(vids_select),1);
        vids_indiv = vids_select;

        for iv = 1:length(vids_select)

            v = find(ismember(vid_recs, vids_select{iv}));

            if isempty(v), continue, end

            m_varx_mov{iv} = varx(Y{v},na,X{v},nb,gamma);
            m_varx_mov_hfa{iv} = varx(HFA{v},na_hfa,X{v},nb,gamma);

        end

        % If Despicable Me English was recorded also compute model for a 5 min version
        if ismember('Despicable_Me_English', vid_recs)

            idx_vid = ismember(vid_recs, 'Despicable_Me_English');

            Y_dme5 = Y{idx_vid}(1:fs*300, :);
            HFA_vid = HFA{idx_vid}(1:fs*300, :);
            X_dme5 = X{idx_vid}(1:fs*300, :);

            m_varx_mov{end+1} = varx(Y_dme5,na,X_dme5,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_vid,na_hfa,X_dme5,nb,gamma);
            vids_indiv{end+1} = 'Despicable_Me_English_5min';

            % Unaligned features
            idx_shift = round(length(X_dme5)/2);
            X_dme5_shift = [X_dme5(idx_shift+1:end, :); X_dme5(1:idx_shift, :)];

            m_varx_mov{end+1} = varx(Y_dme5,na,X_dme5_shift,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_vid,na_hfa,X_dme5_shift,nb,gamma);
            vids_indiv{end+1} = 'Despicable_Me_English_5min_shift';

            % Second half of Despicable Me
            Y_dme5_end = Y{idx_vid}(end-fs*300:end-1, :);
            HFA_vid_end = HFA{idx_vid}(end-fs*300:end-1, :);
            X_dme5_end = X{idx_vid}(end-fs*300:end-1, :);

            m_varx_mov{end+1} = varx(Y_dme5_end,na,X_dme5_end,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_vid_end,na_hfa,X_dme5_end,nb,gamma);
            vids_indiv{end+1} = 'Despicable_Me_English_last_5min';

            % Unaligned features
            idx_shift = round(length(X_dme5_end)/2);
            X_dme5_end_shift = [X_dme5_end(idx_shift+1:end, :); X_dme5_end(1:idx_shift, :)];

            m_varx_mov{end+1} = varx(Y_dme5_end,na,X_dme5_end_shift,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_vid_end,na_hfa,X_dme5_end_shift,nb,gamma);
            vids_indiv{end+1} = 'Despicable_Me_English_last_5min_shift';

            % Edit model with resting state data to contain the same features
            if ismember('Resting_fixation', vid_recs)

                idx_rest = ismember(vid_recs, 'Resting_fixation');
                idx_rest_indiv = ismember(vids_indiv, 'Resting_fixation');

                Y_rest = Y{idx_rest}(1:fs*300, :);
                HFA_rest = HFA{idx_rest}(1:fs*300, :);
                X_rest = X{idx_rest};
                X_rest = X_rest(1:fs*300, :);

                m_varx_mov{idx_rest_indiv} = varx(Y_rest,na,X_rest,nb,gamma);
                m_varx_mov_hfa{idx_rest_indiv} = varx(HFA_rest,na_hfa,X_rest,nb,gamma);

                % Unaligned features (also fixations, rest is already
                % unaligned)
                idx_shift = round(length(X_rest)/2);
                X_rest_shift = [X_rest(idx_shift+1:end, :); X_rest(1:idx_shift, :)];

                m_varx_mov{end+1} = varx(Y_rest, na, X_rest_shift, nb, gamma);
                m_varx_mov_hfa{end+1} = varx(HFA_rest, na_hfa, X_rest_shift, nb, gamma);
                vids_indiv{end+1} = 'Resting_fixation_shift';

            end

        end

        if ismember('Inscapes', vid_recs)

            idx_vid = ismember(vid_recs, 'Inscapes');

            % First half
            Y_insc = Y{idx_vid}(1:fs*300, :);
            HFA_insc = HFA{idx_vid}(1:fs*300, :);
            X_insc = X{idx_vid}(1:fs*300, :);

            m_varx_mov{end+1} = varx(Y_insc,na,X_insc,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_insc,na_hfa,X_insc,nb,gamma);
            vids_indiv{end+1} = 'Inscapes_5min';

            % Unaligned features
            idx_shift = round(length(X_insc)/2);
            X_insc_shift = [X_insc(idx_shift+1:end, :); X_insc(1:idx_shift, :)];

            m_varx_mov{end+1} = varx(Y_insc,na,X_insc_shift,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_insc,na_hfa,X_insc_shift,nb,gamma);
            vids_indiv{end+1} = 'Inscapes_5min_shift';

            % Second half
            Y_insc_end = Y{idx_vid}(end-fs*300:end-1, :);
            HFA_insc_end = HFA{idx_vid}(end-fs*300:end-1, :);
            X_insc_end = X{idx_vid}(end-fs*300:end-1, :);

            m_varx_mov{end+1} = varx(Y_insc_end,na,X_insc_end,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_insc_end,na_hfa,X_insc_end,nb,gamma);
            vids_indiv{end+1} = 'Inscapes_last_5min';

        end

        idx_monkey = cellfun(@(C) contains(C, 'Monkey'), vid_recs);

        if sum(idx_monkey) ~= 0

            idx_vid = find(idx_monkey, 1, 'first');

            Y_monkey = Y{idx_vid}(1:fs*300, :);
            HFA_monkey = HFA{idx_vid}(1:fs*300, :);
            X_monkey = X{idx_vid}(1:fs*300, :);

            m_varx_mov{end+1} = varx(Y_monkey,na,X_monkey,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_monkey,na_hfa,X_monkey,nb,gamma);
            vids_indiv{end+1} = 'Monkey_5min';

            % Unaligned features
            idx_shift = round(length(X_monkey)/2);
            X_monkey_shift = [X_monkey(idx_shift+1:end, :); X_monkey(1:idx_shift, :)];

            m_varx_mov{end+1} = varx(Y_monkey,na,X_monkey_shift,nb,gamma);
            m_varx_mov_hfa{end+1} = varx(HFA_monkey,na_hfa,X_monkey_shift,nb,gamma);
            vids_indiv{end+1} = 'Monkey_5min_shift';


        end

        if ismember('Eyes_Closed_Rest', vid_recs) || ismember('Eyes_Closed_Rest_2', vid_recs)

            idx_rest = cellfun(@(C) contains(C, 'Eyes_Closed_Rest'), vid_recs);
            idx_rest_indiv = ismember(vids_indiv, 'Eyes_Closed_Rest');

            for ir = find(idx_rest)'

                if size(X{ir},1) < fs*300
                    nd = round((size(X{ir},1) - fs*150) / 2);
                else
                    nd = round((size(X{ir},1) - fs*300) / 2);
                end

                Y{ir} = Y{ir}(nd:end-nd, :);
                HFA{ir} = HFA{ir}(nd:end-nd, :);
                X{ir} = X{ir}(nd:end-nd, :);

            end

            m_varx_mov{idx_rest_indiv} = varx(Y(idx_rest),na,X(idx_rest),nb,gamma);
            m_varx_mov_hfa{idx_rest_indiv} = varx(HFA(idx_rest),na_hfa,X(idx_rest),nb,gamma);

        end

        vid_recs = cellfun(@(C) strrep(C, '.mat', ''), vids_indiv, 'UniformOutput', false);

        idx_empty = cellfun(@(C) isempty(C), m_varx_mov);
        m_varx_mov(idx_empty) = [];
        m_varx_mov_hfa(idx_empty) = [];
        vid_recs(idx_empty) = [];


        %% Compute MIMO FIR model (mTRF)
        Rxx=0; Rxy=0;
        for i=1:length(X)
            if idx_comb(i)
                [Rxx_,Rxy_] = myxcorr(X{i},Y{i},L);
                Rxx = Rxx+Rxx_; Rxy = Rxy+Rxy_;
            end
        end
        mTRF = inv(Rxx)*Rxy;
        xdim=size(X{1},2);
        ydim=size(Y{1},2);
        mTRF = permute(reshape(mTRF,L,xdim,ydim),[1 3 2]);

        % mTRF for HFA
        Rxx=0; Rxy=0;
        for i=1:length(X)
            if idx_comb(i)
                [Rxx_,Rxy_] = myxcorr(X{i},HFA{i},L);
                Rxx = Rxx+Rxx_; Rxy = Rxy+Rxy_;
            end
        end
        mTRF_hfa = inv(Rxx)*Rxy;
        xdim=size(X{1},2);
        ydim=size(HFA{1},2);
        mTRF_hfa = permute(reshape(mTRF_hfa,L,xdim,ydim),[1 3 2]);

        %% Save data
        save(out_file, 'm_varx', 'm_varx_hfa', 'm_var', 'm_var_hfa', ...
            'H', 'H_hfa', 'mTRF', 'mTRF_hfa', 'm_varx_mov', 'm_varx_mov_hfa', ...
            'm_varx_features', 'm_varx_features_hfa', ...
            'm_varx_additive', 'm_varx_additive_hfa', 'm_varx_example', ...
            'm_varx_ch_cut', 'ch_keep', ...
            'example_features', 'stim_features', 'features_sort', ...
            'vid_recs', 'L', 'fs_neural', 'labels')

    end

end 