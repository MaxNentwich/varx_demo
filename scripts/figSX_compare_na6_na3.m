
%% Compare models with different na

% Data and figure directory
data_dir = '/media/max/Workspace/Code/varx_demo/results';
fig_dir = '../results/figures';

% Select patient
pat = 'NS127_02';

% Set na to plot
na = {'1', '3', '6', '12', '30'};

% Font 
fig_font = 16;

% Signal Type 
signal_type = 'LFP';

%% Load data
mvarx_na = cell(1,length(na));

for i = 1:length(na)

    if strcmp(signal_type, 'LFP')

        load(sprintf('%s/models_na%s/%s_varx_models.mat', data_dir, na{i}, pat), 'm_varx')
        mvarx_na{i} = m_varx;

    elseif strcmp(signal_type, 'HFA')

        load(sprintf('%s/models_na%s/%s_varx_models.mat', data_dir, na{i}, pat), 'm_varx_hfa')
        mvarx_na{i} = m_varx_hfa;

    end

end

%% R-values
varx_R_na = zeros(size(mvarx_na{1}.A_Rvalue,1), size(mvarx_na{1}.A_Rvalue,2),length(na));

for i = 1:length(na)
    R = mvarx_na{i}.A_Rvalue;
    R(eye(size(R)) == 1) = 0;
    varx_R_na(:,:,i) = R;
end

%% Plot
plot_max = max(varx_R_na(:));
plot_range = [0, plot_max];

figure('Position', [3000,950,1350,300]);

tiledlayout(1, length(na));

for i = 1:length(na)

    ax = nexttile;
    
    imagesc(real(varx_R_na(:,:,i)))
    
    clim(0.2*plot_range)
    axis square
    colormap(ax, slanCM('Reds'))
    xlabel('Channels')
    ylabel('Channels')
    title(sprintf('n_a = %s', na{i}))
    
end

cb = colorbar(); 
ylabel(cb,'R' ,'Rotation',90)

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/figSX_na_R_%s_%s.png', ...
    fig_dir, pat, signal_type), 'Resolution', 600)

%% p-values
varx_p_na = zeros(size(mvarx_na{1}.A_pval,1), size(mvarx_na{1}.A_pval,2),length(na));

for i = 1:length(na)
    varx_p_na(:,:,i) = log10(mvarx_na{i}.A_pval);
end

plot_range = [min(varx_p_na(:)), max(varx_p_na(:))];

figure('Position', [3000,500,1350,300]);

t = tiledlayout(1,length(na));

for i = 1:length(na)

    ax = nexttile;
    
    imagesc(varx_p_na(:,:,i))
    
    clim(ax, plot_range)
    axis square
    colormap(ax, slanCM('summer'))
    xlabel('Channels')
    ylabel('Channels')
    title(sprintf('n_a = %s', na{i}))

end

cb = colorbar(); 
ylabel(cb,'log p-value','Rotation',90)

fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/figSX_na_p_%s_%s.png', ...
    fig_dir, pat, signal_type), 'Resolution', 600)

%% Summary plots
na_x = cellfun(@(C) str2double(C), na);

figure('Position', [4350,950,425,300])
plot(na_x, squeeze(mean(mean(real(varx_R_na),1),2)), '*-')
xticks(na_x)
xlabel('n_a')
ylabel('Mean R')
grid on
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/figSX_na_R_summary_%s_%s.png', ...
    fig_dir, pat, signal_type), 'Resolution', 600)

figure('Position', [4350,500,425,300])
plot(na_x, squeeze(sum(sum(varx_p_na < log10(0.001),1),2)) / (size(mvarx_na{1}.A_pval,1) * size(mvarx_na{1}.A_pval,2)), '*-')
xticks(na_x)
xlabel('n_a')
ylabel('Ratio sig. Connections')
grid on
fontsize(gcf, fig_font, 'points')

exportgraphics(gcf, sprintf('%s/figSX_na_p_summary_%s_%s.png', ...
    fig_dir, pat, signal_type), 'Resolution', 600)

