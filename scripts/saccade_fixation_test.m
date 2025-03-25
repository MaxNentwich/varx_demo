
% Load data
load('example_data.mat', 'X', 'Y', 'fs', 'stim_features')

% Figure directory
fig_dir = '/home/max/Desktop/saccade_fixation_test';
if exist(fig_dir, 'dir') == 0, mkdir(fig_dir), end

% Model parameters
ta=0.1; % length of the VARX model filters in seconds
tb=0.6; % length of the FIR model filter in seconds
gamma = 0; 

na=ta*fs;
nb=tb*fs; 

% Simulation parameters
fix_delay = 0.2;
fix_jitter = 0.07;

time = (0:length(X{2})-1) / fs;

% Saccades
idx_saccades = ismember(stim_features, 'saccades');
idx_fixations = ismember(stim_features, 'fixations');

m_varx_saccades = varx(Y{2}, na, X{2}(:,idx_saccades), nb, gamma);

figure
imagesc([0, tb], [1, size(Y{2},2)], m_varx_saccades.B')
xlabel('Time [s]')
ylabel('Channel')
title('Saccades')

exportgraphics(gcf, sprintf('%s/B_saccade_only.png', fig_dir))

% Saccades and fixations
idx_sac_fix = ismember(stim_features, {'saccades', 'fixations'});

m_varx_sac_fix = varx(Y{2}, na, X{2}(:,idx_sac_fix), nb, gamma);

figure
imagesc([0, tb], [1, size(Y{2},2)], m_varx_sac_fix.B(:,:,1)')
xlabel('Time [s]')
ylabel('Channel')
title('Saccades + Fixations')

exportgraphics(gcf, sprintf('%s/B_saccade_and_fixation.png', fig_dir))

t_saccades = time(X{2}(:,idx_saccades) > 1);
t_fixations = time(X{2}(:,idx_fixations) > 1);

figure
hist(1000*(t_fixations - t_saccades), 100)
xlim([0, 250])
xlabel('Delay Saccades - Fixations [ms]')
ylabel('#')
title('Real Fixations')

exportgraphics(gcf, sprintf('%s/hist_saccade_and_fixation.png', fig_dir))

% Simulate fixations with a certain delay and jitter
saccade_vec = X{2}(:,idx_saccades);

t_saccade = time(saccade_vec > 1);

rng(3)
t_fixation = t_saccade + fix_delay + randn(length(t_saccade),1)' * fix_jitter;

fixation_vec = zeros(size(saccade_vec));

sample_fixation = round(interp1(time, 1:length(time), t_fixation));
fixation_vec(sample_fixation) = 1;
fixation_vec = zscore(fixation_vec);

x = [saccade_vec, fixation_vec];

m_varx_fix_sim = varx(Y{2}, na, x, nb, gamma);

figure
imagesc([0, tb], [1, size(Y{2},2)], m_varx_fix_sim.B(:,:,1)')
hold on
y_lim = get(gca, 'YLim');
plot(fix_delay*ones(2,1), y_lim, 'r')
xlabel('Time [s]')
ylabel('Channel')
title(sprintf('Saccades + Simulated Fixations (delay=%ims, jitter=%ims)', 1000*fix_delay, 1000*fix_jitter))

exportgraphics(gcf, sprintf('%s/B_saccade_and_fixation_sim_delay_%i_jitter_%i.png', fig_dir, 1000*fix_delay, 1000*fix_jitter))

figure
hist(1000*(time(sample_fixation) - t_saccade), 100)
xlim([0, 250])
xlabel('Delay Saccades - Fixations [ms]')
ylabel('#')
title(sprintf('Simulated Fixations (delay=%ims, jitter=%ims)', 1000*fix_delay, 1000*fix_jitter))

exportgraphics(gcf, sprintf('%s/hist_saccade_and_fixation_sim_delay_%i_jitter_%i.png', fig_dir, 1000*fix_delay, 1000*fix_jitter))
