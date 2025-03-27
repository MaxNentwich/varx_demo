addpath('../src/varx')
addpath('../src/Violinplot-Matlab-master/'); % https://github.com/bastibe/Violinplot-Matlab/tree/master

subjects = dir('/media/max/Workspace/Data/varx_data/NS*');
condition = {'Despicable_Me_English_5min.mat', 'Resting_fixation.mat',};
res_dir = '../results/innovation';

if exist(res_dir, 'dir') == 0
    mkdir(res_dir)
end

% We will be plotting one subject at a time, and at the end add a summary.
% picking nice example subject to go last, so we have it in summayr figure.
last = 1;
subject_order = 1:length(subjects); subject_order(last)=[]; subject_order=[subject_order last];

font_size = 11;

% Poster settings
poster_size = false;

if poster_size
    font_size = 26;
end

y2_diff = cell(1, length(subject_order));
s2_diff = cell(1, length(subject_order));
Change_responsive = cell(1, length(subject_order));
Change_not_responsive = cell(1, length(subject_order));

for subj = subject_order

    % first load two conditions, and find electrodes where power 
    % does not change by more than 20dB, i.e. rule out contact problem or
    % movement artefacts
    y2 = cell(1,length(condition));
    for i=1:length(condition)
        file = [subjects(subj).folder '/' subjects(subj).name '/' condition{i}];
        if exist(file, 'file') == 0, continue, end
        load(file,'hfa');
        y2{i} = var(hfa);
    end

    if sum(cellfun(@(C) isempty(C), y2)) ~= 0, continue, end
    
    goodchannels = abs(db(y2{1}./y2{2})/2)<20; 

    % now analyze the data in the good channels
    model = cell(1,length(condition));

    for i=1:length(condition)
    
        warning off % i know some of these variables are missing in rest condition
        load([subjects(subj).folder '/' subjects(subj).name '/' condition{i}],'hfa','fs','fixations','film_cuts','audio_env');
        warning on

        % The order of loading is important, first movie then rest, so that
        % film_cuts and audio_env are defined also for rest. We want to
        % regress out the same signal in both conditions to have same
        % number of free parameters in both conditions.
        if size(fixations,2) > size(fixations,1), fixations = fixations'; end
        if size(film_cuts,2) > size(film_cuts,1), film_cuts = film_cuts'; end
        if size(audio_env,2) > size(audio_env,1), audio_env = audio_env'; end

        x = [fixations film_cuts audio_env]; 
        y = hfa(:,goodchannels); clear hfa;
        
        % subtract mean to avoid dealing with it in the regression
        x = x-mean(x); y=y-mean(y);

        % power of the signal
        y2{i} = var(y);

        % normalize scale for output so that penalty affects them equally
        y = y./sqrt(y2{1});
        x = x./std(x);

        % varx model
        na = 4; nb=fs/2; lambda=0;
        model{i} = varx(y,na,x,nb,lambda);

    end

    noise_diff = db(model{1}.s2./y2{1}') - db(model{2}.s2./y2{2}');

    % remember for each subject change in power for signal and noise
    y2_diff{subj} = db(y2{1})-db(y2{2});
    s2_diff{subj} = db(model{1}.s2)-db(model{2}.s2);

    % report change in power for responsive and non-responsive electrodes
    responsive=sum(model{1}.B_pval<0.01/size(y,2)/size(x,2),2)>0; % bonferony corrected
    Change_responsive    {subj} = noise_diff( responsive);
    Change_not_responsive{subj} = noise_diff(~responsive);
    if sum(responsive)>0
        pval = ranksum(noise_diff(responsive),noise_diff(~responsive));
    else pval=1; end
    fprintf([subjects(subj).name ...
        ': responsive '   num2str(median(noise_diff( responsive)),2), ...
        ' non responsive ' num2str(median(noise_diff(~responsive)),2), ...
        ', p=' num2str(pval,2) ...
        ' N=' num2str(sum(responsive)) ' responsive electrodes.\n'])


end

% Plot last patient
figure(1); tiledlayout(3,3)

h1(1)=nexttile(1,[1 2]);

plot(model{1}.B_Rvalue, 'LineWidth', 1.5); axis tight
ylabel('R')
title('B effect during Movie');
fontsize(gca, font_size, 'points')
grid on 

% display relative change in signal power for all electrodes
h1(2)=nexttile(4,[1 2]); 
plot(noise_diff, 'LineWidth', 1.5); axis tight; title('e^2/y^2 (Movie - Rest): BHA')
ylabel('\DeltadB')
ax = axis; hold on; plot(ax(1:2),[0 0],'k'); hold off

fontsize(gca, font_size, 'points')
grid on 

%% summary of change in power of relative noise in responsive and non responsive electrodes
h1(4)=nexttile(3,[3 1]);
N=length(Change_responsive);
clear tmp1; for i=1:N, tmp1(i)=median(Change_responsive{i}); end; tmp1=tmp1';
clear tmp2; for i=1:N, tmp2(i)=median(Change_not_responsive{i}); end; tmp2=tmp2'; 
pval = ranksum(tmp1,tmp2);
fprintf(['Change in power ratio responsive vs non responsive channesl, p=' num2str(pval,2) ', N=' num2str(N) ' subjects.'])

violinplot([tmp1;tmp2],[repmat("responsive",size(tmp1,1),1);repmat("not responsive",size(tmp2,1),1)]); 
ax=axis; plot(ax(1:2),[0 0],'k'); hold off 
title('Median Channel')
ylabel('e^2/y^2 change (dB)')

fontsize(gca, font_size, 'points')
grid on 

%% simulation with gain adaptation
gamma = 0.001; e_std=0.5;
y_stim = varx_simulate(model{subj}.B,model{subj}.A,x,e_std,gamma);  
model_stim = varx(y_stim,na,x,nb);
y_rest = varx_simulate(model{subj}.B,model{subj}.A,zeros(size(x)),e_std,gamma); 
model_rest = varx(y_rest,na,x,nb);

% Compute change in noise level
noise_diff = db(model_stim.s2'./var(y_stim)) - db(model_rest.s2'./var(y_stim));

%
h1(3)=nexttile(7,[1 2]); 
plot(noise_diff, 'LineWidth', 1.5); axis tight; title('e^2/y^2 (Movie - Rest): Simulated')
xlabel('Channels'); ylabel('\DeltadB')
ax = axis; hold on; plot(ax(1:2),[0 0],'k'); hold off

fontsize(gca, font_size, 'points')
grid on 

if poster_size
    set(gcf, 'Units', 'inches', 'Position', [1,1,15,10])
else
    set(gcf, 'Units', 'inches', 'Position', [5.5,4,6,4])
end

exportgraphics(gcf, '../results/figures/noise-power-relative-change_responsive-HFA.png', ...
    'Resolution', 600)
