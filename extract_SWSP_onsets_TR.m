function extract_SWSP_onsets_TR(savepath, fname_write, coupled_spindles, coupled_SW, Uncoupled_SW, EEG, computeTR, TR, dur2TR)
% extract_SWSP_onsets_TR - Extracts SW-SP event onsets and converts latencies to TR.
%
% Description:
%   This function extracts slow wave-spindle (SW-SP) event latencies and durations.
%   If fMRI alignment is required (computeTR=1), events are aligned with a reference
%   (V1 event), converted to TR units, and saved for further analyses.
%
% Usage:
%   extract_SWSP_onsets_TR(savepath, fname_write, coupled_spindles, coupled_SW, Uncoupled_SW, EEG, computeTR, TR, dur2TR);
%
% Parameters:
%   savepath         - Path where extracted onset files will be saved.
%   fname_write      - Name of the dataset being processed.
%   coupled_spindles - Struct containing detected coupled spindles.
%   coupled_SW       - Struct containing coupled slow waves.
%   Uncoupled_SW     - Struct containing uncoupled slow waves.
%   EEG              - EEGLAB dataset structure.
%   computeTR        - Boolean flag (1 = align with V1 and convert to TR, 0 = keep raw).
%   TR               - TR value for fMRI alignment (in seconds).
%   dur2TR           - Boolean flag (1 = convert durations to TR, 0 = keep in seconds).
%
% Outputs:
%   - Saves .mat and .xlsx files with extracted event data.
%   - If computeTR=1, events are converted to TR units before saving.
%
% Notes:
%   - If no V1 event is found, data is saved without alignment.
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

    if computeTR == 1
        %% Find the V1 event (reference)
        try
        V1idx = find(ismember({EEG.event.type}, {'V  1'}), 1, 'first'); % Find first 'V 1'
        catch
        V1idx = find(ismember({EEG.event.type}, {'R149'}), 1, 'first'); % Find first 'V 1'
        end
        if isempty(V1idx)
            warning('No V1 event found. Data will not be aligned to TR.');
            V1lat = 0; % No reference, keep data as is
        else
            V1lat = EEG.event(V1idx).latency;
        end
    end

    %% Process Coupled Spindles
        CSP = coupled_spindles([coupled_spindles.coup] == 1);
    if ~isempty(CSP)
        CSP_table = struct2table(CSP);
        
        if computeTR == 1
            CSP_table.latency = (CSP_table.latency - V1lat) / EEG.srate / TR;
            if dur2TR
                CSP_table.duration = (CSP_table.duration / EEG.srate) / TR;
                CSP_table.lag = (CSP_table.lag / EEG.srate) / TR;
            else
                CSP_table.duration = (CSP_table.duration / EEG.srate);
                CSP_table.lag = (CSP_table.lag / EEG.srate);
            end
        end
        
        CSP_table = sortrows(CSP_table, 'latency', 'ascend');
        save(fullfile(savepath, [fname_write(1:end-4), '_CSPOnsets.mat']), 'CSP_table');
        writetable(CSP_table, fullfile(savepath, [fname_write(1:end-4), '_CSPOnsets.xlsx']));
    else
        warning('No coupled spindles found for %s', fname_write);
    end

    %% Process Uncoupled Spindles
    UNCSP = coupled_spindles([coupled_spindles.coup] == 0);
    if ~isempty(UNCSP)
        UNCSP_table = struct2table(UNCSP);
        
        if computeTR == 1
            UNCSP_table.latency = (UNCSP_table.latency - V1lat) / EEG.srate / TR;
            if dur2TR
                UNCSP_table.duration = (UNCSP_table.duration / EEG.srate) / TR;
            else
                UNCSP_table.duration = (UNCSP_table.duration / EEG.srate);
            end
        end
        
        UNCSP_table = sortrows(UNCSP_table, 'latency', 'ascend');
        save(fullfile(savepath, [fname_write(1:end-4), '_UNCSPOnsets.mat']), 'UNCSP_table');
        writetable(UNCSP_table, fullfile(savepath, [fname_write(1:end-4), '_UNCSPOnsets.xlsx']));
    end

    %% Process Coupled Slow Waves
    if ~isempty(coupled_SW)
        CSW_table = struct2table(coupled_SW);
        
        if computeTR == 1
            CSW_table.latency = (CSW_table.latency - V1lat) / EEG.srate / TR;
            if dur2TR
                CSW_table.duration = (CSW_table.duration / EEG.srate) / TR;
            else
                CSW_table.duration = (CSW_table.duration / EEG.srate);
            end
        end
        
        CSW_table = sortrows(CSW_table, 'latency', 'ascend');
        save(fullfile(savepath, [fname_write(1:end-4), '_CSWOnsets.mat']), 'CSW_table');
        writetable(CSW_table, fullfile(savepath, [fname_write(1:end-4), '_CSWOnsets.xlsx']));
    end

    %% Process Uncoupled Slow Waves
    if ~isempty(Uncoupled_SW)
        UNCSW_table = struct2table(Uncoupled_SW);
        
        if computeTR == 1
            UNCSW_table.latency = (UNCSW_table.latency - V1lat) / EEG.srate / TR;
            if dur2TR
                UNCSW_table.duration = (UNCSW_table.duration / EEG.srate) / TR;
            else
                UNCSW_table.duration = (UNCSW_table.duration / EEG.srate);
            end
        end
        
        UNCSW_table = sortrows(UNCSW_table, 'latency', 'ascend');
        save(fullfile(savepath, [fname_write(1:end-4), '_UNCSWOnsets.mat']), 'UNCSW_table');
        writetable(UNCSW_table, fullfile(savepath, [fname_write(1:end-4), '_UNCSWOnsets.xlsx']));
    end

    disp(['Onset data saved for ', fname_write]);

end
