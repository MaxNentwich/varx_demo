
%% Example of input filters 

% Select a patient
pat = 'NS127_02';

% Select signal type
signal_type = 'LFP';

% Name of inputs
input_names = {'Fixation Onset', 'Film Cuts', 'Audio Envelope'}; 

% Font size
fig_font = 9;

% Figure directory
fig_dir = '../results/figures';

%% Load data
if strcmp(signal_type, 'LFP')
    load(sprintf('../results/models/%s_varx_models.mat', pat), 'm_varx', 'labels', 'fs_neural')
elseif strcmp(signal_type, 'HFA')
    load(sprintf('../results/models/%s_varx_models.mat', pat), 'm_varx_hfa', 'labels', 'fs_neural')
    m_varx = m_varx_hfa;
end

% Find example channels (biggest effect size)
[~, ch_max] = max(m_varx.B_Rvalue);

% Time 
time = (1:size(m_varx.B,1)) / fs_neural;

n_feature = size(m_varx.B,3);

figure('Units','inches', 'Position', [10,9,10,2])
fontsize(gcf, fig_font, 'points')

tiledlayout(1, n_feature);

for f = 1:n_feature

    ax = nexttile;
    
    plot(time, m_varx.B(:,ch_max(f),f), 'LineWidth', 1)
    
    xlabel('Delay to input [s]')
    ylabel('Weight')
    title(input_names{f})
    grid on
    grid minor

end

exportgraphics(gcf, sprintf('%s/fig3_addon_filters_B_%s.png', fig_dir, signal_type), 'Resolution', 600)
