function finalize_SWSP_results(savepath, filenames, All_subjects_table, ...
    All_subjects_csp, All_subjects_uncsp, all_rep_bins, all_rep_bins_SP, ...
    All_subjects_csw, All_subjects_uncsw, window_size, nBins)
% finalize_SWSP_results - Aggregates and summarizes SW-SP coupling results.
%
% Description:
%   This function compiles all processed subjects' data into a summary table.
%   It computes overall statistics, stores detection results, and exports structured
%   reports in Excel and .mat formats.
%
% Usage:
%   finalize_SWSP_results(All_subjects_table, savepath);
%
% Parameters:
%   All_subjects_table - Cell array containing event statistics for all processed subjects.
%   savepath           - Path where summary files will be saved.
%
% Outputs:
%   - Saves an Excel file (`All_subjects_summary.xlsx`) containing compiled results.
%   - Saves `.mat` files with relevant detection metrics.
%
% Notes:
%   - The output table includes key parameters like SW count, spindle count, coupling ratios, etc.
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


    %% Define Header for Summary Table
    header = {'Subject', 'Sw_neg count', 'SP_count', 'Uncoupled SPS', 'Coupled SPs', ...
        'Coupled SPs total', 'Coupled SPs downstate start to zero', 'Coupled SPs upstate zero to end', ...
        'Coupled SPs downstate -0.5s to zero', 'Coupled SPs upstate zero to 0.5s', ...
        'Coupled SPs amplitude', 'Coupled SPs duration', 'Coupled SPs density', ...
        'Uncoupled SPs amplitude', 'Uncoupled SPs duration', 'Uncoupled SPs density', ...
        'Coupled SWs', 'Coupled SWs duration', 'Coupled SWs density', ...
        'Uncoupled SWs', 'Uncoupled SWs duration', 'Uncoupled SWs density'};

    %% Combine All Data into a Table
    complete_sp_resume = [All_subjects_table, squeeze(struct2cell(All_subjects_csp))', ...
        squeeze(struct2cell(All_subjects_uncsp))', squeeze(struct2cell(All_subjects_csw))', ...
        squeeze(struct2cell(All_subjects_uncsw))'];

    %% Save the Summary Table
    summary_filename_xlsx = fullfile(savepath, 'All_subjects_summary.xlsx');
    summary_filename_mat  = fullfile(savepath, 'All_subjects_rep_bins.mat');
    writetable(cell2table([header; complete_sp_resume]), summary_filename_xlsx);
    save(summary_filename_mat, 'all_rep_bins');

    disp(['Final summary saved: ', summary_filename_xlsx]);

    %% Save Event Correlation Data
    rep_bins_filename_xlsx = fullfile(savepath, 'All_subjects_rep_bins.xlsx');
    rep_bins_sp_filename_xlsx = fullfile(savepath, 'All_subjects_rep_bins_sppercent.xlsx');
    save(fullfile(savepath, 'All_subjects_rep_bins.mat'), 'all_rep_bins');
    save(fullfile(savepath, 'All_subjects_rep_bins_sppercent.mat'), 'all_rep_bins_SP');
    writetable(array2table(squeeze(all_rep_bins)), rep_bins_filename_xlsx);
    writetable(array2table(squeeze(all_rep_bins_SP)), rep_bins_sp_filename_xlsx);

    %% Create Group-Level Event Correlation Histogram
    xAxis = linspace(-window_size, window_size, nBins);
    figure;
    bar(xAxis, squeeze(mean(all_rep_bins_SP)));
    axis tight;
    title({'Event correlation histogram for all subjects', ''});
    xlabel({'', 'Time (s)'});
    ylabel({'SW-SP co-occurrence (%)', ''});
    ylim auto;
    saveas(gcf, fullfile(savepath, 'SW-SP_coupling_all_subjects.tiff'), 'tif');
    close;

    disp('Final group summary and figures saved.');

end
