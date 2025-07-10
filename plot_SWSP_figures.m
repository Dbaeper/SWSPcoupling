function plot_SWSP_figures(savepath, fname_write, window_size, nBins, EEG, interval, rep_bins, rep_bins_sp_percent)
% plot_SWSP_figures - Generates and saves visualizations for SW-SP coupling.
%
% Description:
%   This function produces various plots related to slow wave-spindle coupling, including:
%   - Time-domain averaged SW-SP windows.
%   - Event correlation histograms.
%   - Other relevant statistical visualizations.
%
% Usage:
%   plot_SWSP_figures(All_subjects_table, savepath);
%
% Parameters:
%   All_subjects_table - Struct containing detection results for visualization.
%   savepath           - Path where figures will be saved.
%
% Outputs:
%   - Saves plots as `.tiff` and `.png` in the specified directory.
%
%
% Author: Daniel Baena  
% Email: dbaenape@uottawa.ca - dbaeper@gmail.com  
% Affiliation: University of Ottawa  
% -------------------------------------------------------------------------
% This script is part of the SW-SP Coupling Toolbox
%
% If you use this software or its methods in your research, please cite:
%
% Baena, D., Ray, L.B., & Fogel, S.M. (2025).
% A novel adaptive time‑window method for detecting slow wave–spindle coupling:
% Comparison of temporal co‑occurrence and phase–amplitude coupling approaches.
% Journal of Neuroscience Methods, 422, 110526.
% https://doi.org/10.1016/j.jneumeth.2025.110526
% -------------------------------------------------------------------------

    timeLim = window_size;
    sampFreq = EEG.srate;
    time = -timeLim:(1/sampFreq):timeLim;
    
    % --- 1. Plot SW-SP Coupling Window ---
    figure('Visible', 'off');
    data_mean = mean(interval);
    plot(time, data_mean, '-r', 'LineWidth', 2);
    title(['SW-SP coupling window for ', strrep(fname_write, '_', '\_')], 'Interpreter', 'none');
    xlabel('Seconds');
    ylabel('Micro Volts');
    xx = xlim;
    line([xx(1) xx(2)], [0 0], 'LineStyle', '--');
    saveas(gcf, fullfile(savepath, [fname_write(1:end-4), '_SW-SP_coupling_window.tiff']), 'tif');
    % close;
    
    % --- 2. Event Correlation Histogram (Hz) ---
    xAxis = linspace(-window_size, window_size, nBins);
    figure('Visible', 'off');
    bar(xAxis, rep_bins);
    axis tight;
    title(['Event correlation histogram for ', strrep(fname_write, '_', '\_')], 'Interpreter', 'none');
    xlabel('Time (sec)');
    ylabel('Hz');
    saveas(gcf, fullfile(savepath, [fname_write(1:end-4), '_Event_correlation_histogram_sw.tiff']), 'tif');
    % close;

    % --- 3. Event Correlation Histogram (Spindle Percent) ---
    figure('Visible', 'off');
    bar(xAxis, rep_bins_sp_percent);
    axis tight;
    title(['Event correlation histogram for ', strrep(fname_write, '_', '\_')], 'Interpreter', 'none');
    xlabel('Time (sec)');
    ylabel('Sp percent');
    saveas(gcf, fullfile(savepath, [fname_write(1:end-4), '_Event_correlation_histogram_sp.tiff']), 'tif');
    % close;
    
    disp(['Figures saved for ', fname_write]);

end
