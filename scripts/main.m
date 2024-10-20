%% Compute all VARX models

%% Define paths and parameters
data_dir = '/media/max/Workspace/Data/movies_legacy_sample/';
et_dir = '/media/DATA/ECoGData/Tobii/Patients';
cut_dir = '../data/scene_cuts_frames';
speech_dir = '../data/speech_files';
sample_data_dir = '/media/max/Workspace/Data/varx_data';

addpath(genpath('../src'))

% models
model_dir = '../results/models';
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

stim_features = {'fixations'};
str = {'fixations'};

videos = {...
    'Despicable_Me_English.mat', ...
    'Despicable_Me_Hungarian.mat', ...
    'Monkey1_Rep_1.mat', ...
    'Monkey2_Rep_1.mat', ...
    'Monkey5_Rep_1.mat', ...
    'The_Present_Rep_1.mat',...
    'The_Present_Rep_2.mat', ...
    'Inscapes.mat',...
    'Resting_fixation.mat'};

% List of videos for individual analysis
vids_select = {'Despicable_Me_English', 'Resting_fixation'};

patients = dir(data_dir);
patients = patients(3:end);

%% Compute all Models
for pat = 1:length(patients)

    out_file = sprintf('%s/%s_varx_models.mat', model_dir, patients(pat).name);

    if compute_models && exist(out_file, 'file') ~= 0
        continue
    end

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
    
        label_dir = sprintf('/media/DATA/ECoGData/Tobii/Patients/%s/matlab_data/Envelope_phase/BHA', patients(pat).name);
    
        load([data_dir patients(pat).name '/LFP/' video],'lfp','fs');
        load([data_dir patients(pat).name '/BHA/' video],'envelope');
        load([data_dir patients(pat).name '/Stimuli/' video], stim_features{:});
        load(sprintf('%s/%s', label_dir, video), 'labels')

        if strcmp(patients(pat).name, 'NS155') && strcmp(video, 'Resting_fixation.mat')
            idx_extra = find(ismember(labels, 'ROs15-ROs16'));
            lfp(:,idx_extra) = [];
            envelope(:,idx_extra) = [];
            labels(idx_extra) = [];
        end

        x = [eval(stim_features{1})];
        y = lfp; 
        fs_neural = fs;
        hfa = envelope;
    
        if length(x) < length(y)
            y = y(1:length(x), :);
            lfp = lfp(1:length(x), :);
            hfa = hfa(1:length(x), :);
        end

        time = (0:length(x)-1) / fs_neural;
    
        if size(x,2) > size(x,1), x = x'; end

        %% Load film cuts
        if strcmp(vid_name_legacy, 'Inscapes') || strcmp(vid_name_legacy, 'Resting_fixation') 
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
    
        clearvars eye
    
        %% Load audio
        speech_vid = sprintf('%s/%s.mat', speech_dir, vid_name_legacy);
    
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
            save(sprintf('%s/%s.mat', data_pat_dir, vid_name_legacy), ...
                'lfp', 'hfa', 'fixations', 'fs');
        elseif strcmp(vid_name_legacy, 'Inscapes')
            audio_env = audio_ds';
            save(sprintf('%s/%s.mat', data_pat_dir, vid_name_legacy), ...
                'lfp', 'hfa', 'audio_env', 'fixations', 'fs');
        else 
            audio_env = audio_ds';
            film_cuts = scenes_vec';
            save(sprintf('%s/%s.mat', data_pat_dir, vid_name_legacy), ...
                'lfp', 'hfa', 'film_cuts', 'audio_env', 'fixations', 'fs');

            if strcmp('Despicable_Me_English', vid_name_legacy)
                lfp = lfp(1:fs*300, :);
                hfa = hfa(1:fs*300, :);
                audio_env = audio_env(1:fs*300);
                film_cuts = film_cuts(1:fs*300);
                fixations = fixations(1:fs*300);
                save(sprintf('%s/Despicable_Me_English_5min.mat', data_pat_dir), ...
                    'lfp', 'hfa', 'film_cuts', 'audio_env', 'fixations', 'fs');
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
        
        % base = basis(n_base,nb,'normal');
        m_varx = varx(Y,na,X,nb,gamma);
        m_varx_hfa = varx(HFA,na_hfa,X,nb,gamma);

        % derive the corresponding impulse responses
        H = varx_trf(m_varx.B, m_varx.A, L);
        H_hfa = varx_trf(m_varx_hfa.B, m_varx_hfa.A, L);
        
        % now lets do just an AR model
        m_var = varx(Y,na,[],[],gamma);
        m_var_hfa = varx(HFA,na_hfa,[],[],gamma);

        % VARX models with individual features
        feature_combos = {'Fixations and Cuts', 'Fixations and Sound', 'Cuts and Sound', ...
                    'Fixations', 'Cuts', 'Sound', 'None'};

        col_shift = {3,2,1,2:3,[1,3],1:2,1:3};

        m_varx_features = cell(length(feature_combos),1);
        m_varx_hfa_features = cell(length(feature_combos),1);

        for feat = 1:length(feature_combos)

            X_shift = cell(1,length(X));
            for ix = 1:length(X)
                idx_shift = round(length(X{ix})/2);
                X_shift{ix} = X{ix};
                X_shift{ix}(:,col_shift{feat}) = [X_shift{ix}(idx_shift+1:end, col_shift{feat}); X_shift{ix}(1:idx_shift, col_shift{feat})];
            end

            m_varx_features{feat} = varx(Y,na,X_shift,nb,gamma);
            m_varx_hfa_features{feat} = varx(HFA,na_hfa,X_shift,nb,gamma);

        end
     
        % Compute models for each movie/condition separately
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

            % Edit model with resting state data to contain the same features
            if ismember('Resting_fixation', vid_recs)

                idx_rest = ismember(vid_recs, 'Resting_fixation');
                idx_rest_indiv = ismember(vids_indiv, 'Resting_fixation');

                Y_rest = Y{idx_rest}(1:fs*300, :);
                HFA_rest = HFA{idx_rest}(1:fs*300, :);
                X_rest = X{idx_rest};
                X_rest = X_rest(1:fs*300, :);
                X_rest(:,2:3) = X_dme5(:,2:3);

                m_varx_mov{idx_rest_indiv} = varx(Y_rest,na,X_rest,nb,gamma);
                m_varx_mov_hfa{idx_rest_indiv} = varx(HFA_rest,na_hfa,X_rest,nb,gamma);
            
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


        end

        vid_recs = cellfun(@(C) strrep(C, '.mat', ''), vids_indiv, 'UniformOutput', false);
        
        idx_empty = cellfun(@(C) isempty(C), m_varx_mov);
        m_varx_mov(idx_empty) = [];
        m_varx_mov_hfa(idx_empty) = [];
        vid_recs(idx_empty) = [];

        %% Compute MIMO FIR model (mTRF)
        Rxx=0; Rxy=0;
        for i=1:length(X)
            [Rxx_,Rxy_] = myxcorr(X{i},Y{i},L);
            Rxx = Rxx+Rxx_; Rxy = Rxy+Rxy_;
        end
        mTRF = inv(Rxx)*Rxy;
        xdim=size(X{1},2);
        ydim=size(Y{1},2);
        mTRF = permute(reshape(mTRF,L,xdim,ydim),[1 3 2]);

        % mTRF for HFA
        Rxx=0; Rxy=0;
        for i=1:length(X)
            [Rxx_,Rxy_] = myxcorr(X{i},HFA{i},L);
            Rxx = Rxx+Rxx_; Rxy = Rxy+Rxy_;
        end
        mTRF_hfa = inv(Rxx)*Rxy;
        xdim=size(X{1},2);
        ydim=size(HFA{1},2);
        mTRF_hfa = permute(reshape(mTRF_hfa,L,xdim,ydim),[1 3 2]);
    
        %% Save data
        save(out_file, 'm_varx', 'm_var', 'H', 'mTRF', 'm_varx_mov', ...
            'm_varx_hfa', 'm_var_hfa', 'H_hfa', 'mTRF_hfa', 'm_varx_mov_hfa', ...
            'm_varx_features', 'm_varx_hfa_features', 'vid_recs', 'feature_combos', ...
            'L', 'fs_neural', 'labels')
        
    end

end 