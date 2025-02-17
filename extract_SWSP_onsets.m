function extract_SWSP_onsets(savepath, fname_write, coupled_spindles, coupled_SW, Uncoupled_SW)
% extract_SWSP_onsets - Extracts slow wave-spindle event onsets for further analyses.
%
% Description:
%   This function extracts key event latencies (onsets), durations, and other parameters
%   for both coupled and uncoupled spindles. The extracted data is saved in a .mat file.
%
% Usage:
%   extract_SWSP_onsets(EEG, coupled_spindles, savepath);
%
% Parameters:
%   EEG              - EEGLAB dataset structure.
%   coupled_spindles - Struct containing detected spindles and coupling information.
%   savepath         - Path where extracted onset files will be saved.
%
% Outputs:
%   - Saves a .mat file containing spindle onset latencies, durations, and coupling status.
%
% Notes:
%   - Extracted latencies are referenced to the original EEG sampling rate.
%
% Author: Daniel Baena  
% Email: dbaenape@uottawa.ca  
% Affiliation: University of Ottawa  
% Date: 2025-02-06


    %% Process Coupled Spindles
    if ~isempty(coupled_spindles)
        CSP_table = struct2table(coupled_spindles);
        CSP_table = sortrows(CSP_table, 'latency', 'ascend'); % Sort by latency
        save(fullfile(savepath, [fname_write(1:end-4), '_SPOnsets.mat']), 'CSP_table');
        writetable(CSP_table, fullfile(savepath, [fname_write(1:end-4), '_SPOnsets.xlsx']));
    else
        warning('No coupled spindles found for %s', fname_write);
    end

    %% Process Uncoupled Spindles
    UNCSP = coupled_spindles([coupled_spindles.coup] == 0);
    if ~isempty(UNCSP)
        UNCSP_table = struct2table(UNCSP);
        UNCSP_table = sortrows(UNCSP_table, 'latency', 'ascend');
        save(fullfile(savepath, [fname_write(1:end-4), '_UNCSPOnsets.mat']), 'UNCSP_table');
        writetable(UNCSP_table, fullfile(savepath, [fname_write(1:end-4), '_UNCSPOnsets.xlsx']));
    end

    %% Process Coupled Slow Waves
    if ~isempty(coupled_SW)
        CSW_table = struct2table(coupled_SW);
        CSW_table = sortrows(CSW_table, 'latency', 'ascend');
        save(fullfile(savepath, [fname_write(1:end-4), '_CSWOnsets.mat']), 'CSW_table');
        writetable(CSW_table, fullfile(savepath, [fname_write(1:end-4), '_CSWOnsets.xlsx']));
    end

    %% Process Uncoupled Slow Waves
    if ~isempty(Uncoupled_SW)
        UNCSW_table = struct2table(Uncoupled_SW);
        UNCSW_table = sortrows(UNCSW_table, 'latency', 'ascend');
        save(fullfile(savepath, [fname_write(1:end-4), '_UNCSWOnsets.mat']), 'UNCSW_table');
        writetable(UNCSW_table, fullfile(savepath, [fname_write(1:end-4), '_UNCSWOnsets.xlsx']));
    end

    disp(['Onset data saved for ', fname_write]);

end
