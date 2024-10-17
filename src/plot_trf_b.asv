%% Compare B filter from VARX model and H (including internal dynamics
function plot_trf_b(m_varx, H, fs, idx_visual, idx_stim, x_label_str, title_str, poster_size, fig_font, out_file)
    
    [nb, ~, ~] = size(m_varx.B);
    time_B = (1:nb)/fs;                     
    time_H = (1:size(H,1))/fs;

    plot_range = max(max(abs(H(:,idx_visual,idx_stim))));
    plot_clim = [-plot_range, plot_range];

    if poster_size
        fig = figure('Units', 'inches', 'Position', [1,1,8,5]);
    else
        fig = figure('Position', [150,300,850,425]);
    end

    if poster_size
        t = tiledlayout(1,2, 'TileSpacing', 'compact');
    else
        t = tiledlayout(1,14, 'TileSpacing', 'compact');
    end

    if ~poster_size

        ax1 = nexttile(1, [1,2]);
        imagesc(log10(m_varx.B_pval(idx_visual, idx_stim)))
        
        xticks([])
        set(ax1,'YAxisLocation','right')
        ylabel('Channel')
        cb = colorbar('Location', 'westoutside'); 
        
        ylabel(cb,'log p-value','Rotation',90)
        colormap(ax1, slanCM('summer'))

    end

    if poster_size
        ax2 = nexttile;
    else
        ax2 = nexttile(3, [1,6]);
    end

    imagesc([time_B(1), time_B(end)], [1, sum(idx_visual)], m_varx.B(:,idx_visual,idx_stim)')
    
    clim(ax2, plot_clim)
    xlim([time_H(1), time_H(end)])
    title('B')
    xlabel(x_label_str)
    yticks([])
    clim(ax2, plot_clim)
    colormap(ax2, slanCM('bwr'))
    grid on

    if poster_size
        ax3 = nexttile;
    else
        ax3 = nexttile(9, [1,6]);
    end

    imagesc([time_H(1), time_H(end)], [1, sum(idx_visual)], H(:,idx_visual,idx_stim)')
    
    clim(ax3, plot_clim)
    xlim([time_H(1), time_H(end)])
    title('H')
    xlabel(x_label_str)
    cb = colorbar();
    yticks([])
    ylabel(cb,'Weight','Rotation',90)
    clim(ax2, plot_clim)
    colormap(ax3, slanCM('bwr'))
    grid on
    
    title(t, title_str)
    fontsize(fig, fig_font, 'points')
    
    exportgraphics(fig, out_file, 'Resolution', 600)

    %% Difference of power in significant channels
    B_plot = m_varx.B(:,idx_visual,idx_stim)';
    H_plot = H(:,idx_visual,idx_stim)';

    B_avg_pow = mean(B_plot.^2, 2);
    H_avg_pow = mean(H_plot.^2, 2);

    figure('Position', [500,275,300,420])
    hold on 

    for i = 1:length(B_avg_pow)
        
        plot([1,2], [B_avg_pow(i), H_avg_pow(i)], 'k')
        plot(1, B_avg_pow(i), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 7)
        plot(2, H_avg_pow(i), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 7)

    end

    plot([1,2], [mean(B_avg_pow), mean(H_avg_pow)], 'k', 'LineWidth', 4)
    plot(1, mean(B_avg_pow), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    plot(2, mean(H_avg_pow), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k')

    xlim([0.75, 2.25])
    xticks([1,2])
    xticklabels({'B', 'H'})

    ylabel('Average Power [a.u.]')
    title('Power')

    grid on

    ax = ancestor(gca, 'axes');
    ax.YAxis.Exponent = 0;
    ytickformat('%0.4f')

    fontsize(16, 'points')
    
    %% Length of responses
    fwhm_H = resp_length(H_plot, fs, 10, 5) / fs * 1e3;
    fwhm_B = resp_length(B_plot, fs, 10, 5) / fs * 1e3;
    
    figure('Position', [500,275,300,420])
    hold on 

    for i = 1:length(fwhm_B)
        
        plot([1,2], [fwhm_B(i), fwhm_H(i)], 'k')
        plot(1, fwhm_B(i), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 7)
        plot(2, fwhm_H(i), 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 7)

    end

    plot([1,2], [mean(fwhm_B), mean(fwhm_H)], 'k', 'LineWidth', 4)
    plot(1, mean(fwhm_B), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
    plot(2, mean(fwhm_H), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k')

    xlim([0.75, 2.25])
    xticks([1,2])
    xticklabels({'B', 'H'})

    ylabel('Length of responses [ms]')
    title('Length')
    
    grid on

    fontsize(16, 'points')
    
    exportgraphics(gcf, strrep(out_file, '.png', '_fwhm.png'), 'Resolution', 300)

end