function All_subjects_table = SWSPcoupling_adapt(EEG, loadpath, filenames, stage_tag, savepath)
% SWSPcoupling_adapt - Performs slow wave-spindle coupling analysis using an adaptive time window.
%
% Description:
%   Similar to `SWSPcoupling_fix`, but instead of using a **fixed time window**, this
%   function dynamically adjusts the time window **based on each detected slow wave duration**.
%
% Usage:
%   [All_subjects_table] = SWSPcoupling_adapt(EEG, loadpath, filenames, stage_tag, savepath);
%
% Parameters:
%   EEG        - EEGLAB dataset structure.
%   loadpath   - Path to the EEG dataset files.
%   filenames  - Cell array of EEG dataset filenames to process.
%   stage_tag  - Structure containing GUI-selected parameters (see `SWSPcoupling_fix.m`).
%   savepath   - Path where results will be saved.
%
% Outputs:
%   All_subjects_table - Cell array containing detected event statistics.
%
% Notes:
%   - Results are saved in multiple formats: .mat and .xlsx for further processing.
%   - Uses **adaptive windowing**, meaning each slow wave event gets a time window
%     proportional to its duration.
%
% Author: Daniel Baena  
% Email: dbaenape@uottawa.ca  
% Affiliation: University of Ottawa  
% Date: 2025-02-06

filenames = filenames';
All_subjects_table = {zeros(size(filenames, 1), 3)};
for s = 1:size(filenames, 1)
    try
        %% USER-DEFINED PARAMETERS
        ChOI = stage_tag.ChOI; % e.g., {'Fz','Cz','Pz','C3','C4'}
        eventName = stage_tag.eventName; % e.g., 'SWneg'
        stages = stage_tag.stages;       % Sleep stages of interest
        window_size = stage_tag.window_size;  % In seconds (used if fixed window mode is selected)
        nBins = stage_tag.nBins;         % Number of bins for histogram
        
        %% FILTER SETTINGS (for phase extraction)
        LPorder = 32; % lowpass filter order
        LPfreq = 4.3; % lowpass frequency cutoff
        HPorder = 64; % highpass filter order
        HPfreq = 0.46; % highpass frequency cutoff
        filtAttn = 80; % stopband attenuation (dB)
        
        %% Load the EEGlab dataset
        EEG = pop_loadset([loadpath, filenames{s}]);
        EEG.srate=round(EEG.srate);
        %% Find channels of interest
        ChName = {EEG.chanlocs.labels};
        nChOI = false(size(ChName));
        for n = 1:length(ChOI)
            chIdx(n) = find(strcmp(ChName, ChOI{n}));
            nChOI = logical(nChOI + strcmp(ChName, ChOI{n}));
        end
        if ~any(nChOI)
            error('Missing or incorrect channel label input')
        end
        data = double(EEG.data(chIdx, :)); % Use double precision for filtering.
        clear n nChOI
        
        %% Filter the EEG channels of interest
        disp('Filtering the data...')
        lp = design(fdesign.lowpass('N,Fst,Ast', LPorder, LPfreq, filtAttn, EEG.srate), 'cheby2');
        hp = design(fdesign.highpass('N,Fst,Ast', HPorder, HPfreq, filtAttn, EEG.srate), 'cheby2');
        datafilt = zeros(size(data)); % Pre-allocate
        parfor n = 1:size(data, 1)
            filtCh = filtfilt(lp.sosMatrix, lp.ScaleValues, data(n, :));
            filtCh = filtfilt(hp.sosMatrix, hp.ScaleValues, filtCh);
            datafilt(n, :) = filtCh;
        end
        clear hp lp filtCh n data
        
        %% Select SW and spindle events from the dataset
        SW_idx = find(ismember({EEG.event.type}, eventName));
        Sw_wind = EEG.event(SW_idx);
        Spindles_idx = find(ismember({EEG.event.type}, 'Spindle'));
        Spindles = EEG.event(Spindles_idx);
        
        %% Keep only events in channels of interest
        for chid = 1:length(chIdx)
            if chid == 1
                Chann_idxs = find(strcmp({Sw_wind.channel}, ChOI(chid)) == 1)';
                Chann_idxs_sp = find(strcmp({Spindles.channel}, ChOI(chid)) == 1)';
            else
                Chann_select = find(strcmp({Sw_wind.channel}, ChOI(chid)) == 1)';
                Chann_idxs = [Chann_idxs; Chann_select];
                clear Chann_select
                Chann_select_sp = find(strcmp({Spindles.channel}, ChOI(chid)) == 1)';
                Chann_idxs_sp = [Chann_idxs_sp; Chann_select_sp];
                clear Chann_select_sp
            end
        end
        Sw_wind = Sw_wind(Chann_idxs);
        Spindles = Spindles(Chann_idxs_sp);
        
        %% Remove duplicate spindles (if detected in different channels)
        Spindles_latency = [Spindles.latency]';
        [~, unique_idx] = unique(Spindles_latency, 'stable');
        Spindles = Spindles(unique_idx);
        
        %% Keep only events in the specified sleep stages
        SW_idx = find(~ismember({Sw_wind.SleepStage}, stages));
        Sw_wind(SW_idx) = [];
        SP_idx = find(~ismember({Spindles.SleepStage}, stages));
        Spindles(SP_idx) = [];

        %% Compute total sleep duration (in minutes) for selected sleep stages
        sleep_events = EEG.event(ismember({EEG.event.type}, stages)); % Select events matching sleep stages
        sleep_duration_total_datapoint = sum([sleep_events.duration]); % Sum their durations (in seconds)
        sleep_duration_total_sec = sleep_duration_total_datapoint / EEG.srate; % Convert to minutes
        sleep_duration_total_min = sleep_duration_total_sec / 60;
        %% Start coupled spindles identification
        run = 1;
        global_run = 1;
        coup_histogram = zeros(1, [window_size * 2 * EEG.srate + 1]);
        interval = zeros(1, window_size * 2 * EEG.srate + 1);
        
        Current_event_window_orig = Sw_wind;
        Current_spindle_orig = Spindles;
        for Stage_id = 1:length(stages)
            event_window_idx = find(ismember({Current_event_window_orig.SleepStage}, stages{Stage_id}));
            Current_event_window = Current_event_window_orig(event_window_idx);
            Sp_idx = find(ismember({Current_spindle_orig.SleepStage}, stages{Stage_id}));
            Current_spindle = Current_spindle_orig(Sp_idx);
            
            %% Build time windows around each SW event using fixed or adaptive mode
            range_windows = zeros(size(Current_event_window, 2), 5); % Columns: [start, end, center, channel index, duration]
            if ~isempty(Current_spindle) && ~isempty(Current_event_window)
                for x = 1:size(Current_event_window, 2)
                    center = Current_event_window(x).latency + Current_event_window(x).peak;
                    chanIdx = find(strcmp(Current_event_window(x).channel, ChName));
                         % Adaptive mode: use event's own duration
                        halfDuration = round(Current_event_window(x).duration / 2);
                        range_windows(x, 1) = center - halfDuration;
                        range_windows(x, 2) = center + halfDuration;
                        range_windows(x, 3) = center;
                        range_windows(x, 4) = chanIdx;
                        range_windows(x, 5) = Current_event_window(x).duration;

                end

                for sps = 1:size(Current_spindle, 2)
                    % Define spindle peak as the middle of its duration, only if it does not exist
                    if ~isfield(Current_spindle(sps), 'peak') || isempty(Current_spindle(sps).peak)
                        Current_spindle(sps).peak = round(Current_spindle(sps).latency + (Current_spindle(sps).duration / 2));
                    end
                    for win = 1:size(range_windows, 1)
                        if range_windows(win, 1) < Current_spindle(sps).peak && Current_spindle(sps).peak < range_windows(win, 2)
                            if range_windows(win, 2) > size(datafilt, 2)
                                break
                            end
                            channel = find(range_windows(win, 4) == chIdx);
                            if channel == 0
                                continue
                            end
                            sizetoadd = size(datafilt(channel, range_windows(win, 1):range_windows(win, 2)), 2);
                            midpoint = round(size(interval, 2) / 2);
                            leftside = round(midpoint - (sizetoadd / 2));
                            rightside = round(midpoint + (sizetoadd / 2));
                            
                            try
                                interval(global_run, leftside:rightside) = datafilt(channel, range_windows(win, 1):range_windows(win, 2));
                            catch
                                rightside = rightside - 1;
                                interval(global_run, leftside:rightside) = datafilt(channel, range_windows(win, 1):range_windows(win, 2));
                            end
                            
                            [~, peak_int_pos] = find(range_windows(win, 1):range_windows(win, 2) == Current_spindle(sps).peak);
                            sw_phase = angle(hilbert(datafilt(channel, range_windows(win, 1):range_windows(win, 2))));
                            coup_histogram(1, leftside + peak_int_pos - 1) = coup_histogram(1, leftside + peak_int_pos - 1) + 1;
                            
                            Current_spindle(sps).coup = 1; % Mark as coupled
                            Current_spindle(sps).phase_angle = sw_phase(peak_int_pos);
                            Current_spindle(sps).lag = Current_spindle(sps).peak - range_windows(win, 3);
                            Current_spindle(sps).sw_wind_start = range_windows(win, 1);
                            Current_spindle(sps).sw_wind_end = range_windows(win, 2);
                            Current_spindle(sps).sw_neg_peak = range_windows(win, 3);
                            
                            Current_event_window(win).coup = 1;
                            break
                        else
                            Current_spindle(sps).coup = 0;
                            Current_spindle(sps).phase_angle = [];
                            Current_spindle(sps).lag = nan;
                            Current_spindle(sps).sw_wind_start = [];
                            Current_spindle(sps).sw_wind_end = [];
                            Current_spindle(sps).sw_neg_peak = [];
                        end
                    end
                    global_run = global_run + 1;
                end
                
                if run == 1
                    complete_detect = Current_spindle;
                    try
                        cp_SW_idx = find(~cellfun(@isempty, {Current_event_window.coup}));
                        coupled_SW = Current_event_window(cp_SW_idx);
                        uncp_SW_idx = find(cellfun(@isempty, {Current_event_window.coup}));
                        Uncoupled_SW = Current_event_window(uncp_SW_idx);
                    catch
                        coupled_SW = [];
                        Uncoupled_SW = Current_event_window;
                    end
                    run = 2;
                else
                    complete_detect = [complete_detect, Current_spindle];
                    try
                        cp_SW_idx = find(~cellfun(@isempty, {Current_event_window.coup}));
                        coupled_SW_r = Current_event_window(cp_SW_idx);
                        coupled_SW = rmfield(coupled_SW, 'coup');
                        coupled_SW_r = rmfield(coupled_SW_r, 'coup');
                        coupled_SW = [coupled_SW, coupled_SW_r];
                        uncp_SW_idx = find(cellfun(@isempty, {Current_event_window.coup}));
                        Uncoupled_SW_r = Current_event_window(uncp_SW_idx);
                        Uncoupled_SW = rmfield(Uncoupled_SW, 'coup');
                        Uncoupled_SW_r = rmfield(Uncoupled_SW_r, 'coup');
                        Uncoupled_SW = [Uncoupled_SW, Uncoupled_SW_r];
                    catch
                        Uncoupled_SW_r = Current_event_window;
                        if isfield(Uncoupled_SW, 'coup')
                            Uncoupled_SW = rmfield(Uncoupled_SW, 'coup');
                        end
                        if isfield(Uncoupled_SW_r, 'coup')
                            Uncoupled_SW_r = rmfield(Uncoupled_SW_r, 'coup');
                        end
                        Uncoupled_SW = [Uncoupled_SW, Uncoupled_SW_r];
                    end
                end
            end
        end
        
        %% Save subject detection results
        rep_binsChann = sum(reshape(coup_histogram(1:end-1), size(coup_histogram(1:end-1), 2)/nBins, nBins));
        binSize = (window_size * window_size) / nBins;
        rep_bins = rep_binsChann / size(Sw_wind, 2) / binSize;
        rep_bins_sp_percent = rep_binsChann / size(find(ismember({complete_detect.type}, 'Spindle')), 2) / binSize;
        all_rep_bins(s, :, :) = rep_bins;
        all_rep_bins_SP(s, :, :) = rep_bins_sp_percent;
        
        fname_write = filenames{s};
        writetable(struct2table(complete_detect), [savepath, filesep, fname_write(1, 1:end-4), '_coupled_spindles.xlsx']);
        coupled_spindles = complete_detect;
        save([savepath, filesep, fname_write(1, 1:end-4), '_coupled_spindles.mat'], 'coupled_spindles')
        %% Update EEG.event structure and rename spindle events based on coupling status and eventName
        for idx = 1:length(coupled_spindles)
            spindle_latency = coupled_spindles(idx).latency;
            event_idx = find([EEG.event.latency] == spindle_latency, 1);
            if ~isempty(event_idx)
                EEG.event(event_idx).coup = coupled_spindles(idx).coup;
                EEG.event(event_idx).phase_angle = coupled_spindles(idx).phase_angle;
                EEG.event(event_idx).lag = coupled_spindles(idx).lag;
                EEG.event(event_idx).sw_wind_start = coupled_spindles(idx).sw_wind_start;
                EEG.event(event_idx).sw_wind_end = coupled_spindles(idx).sw_wind_end;
                EEG.event(event_idx).sw_neg_peak = coupled_spindles(idx).sw_neg_peak;

                % Rename event.type based on coupling status and eventName directly
                if coupled_spindles(idx).coup == 1
                    EEG.event(event_idx).type = ['spindle_C_', eventName];
                else
                    EEG.event(event_idx).type = ['spindle_U_', eventName];
                end
            end
        end

        % Save the updated EEGLAB dataset (.set file)
        pop_saveset(EEG, 'filename', [fname_write(1:end-4), '_coupled.set'], 'filepath', savepath);

        %% Populate All_subjects table with subject info
        All_subjects_table{s, 1} = fname_write(1, 1:end-4);
        All_subjects_table{s, 2} = size(Sw_wind, 2);
        All_subjects_table{s, 3} = size(find(ismember({complete_detect.type}, 'Spindle')), 2);
        All_subjects_table{s, 4} = size(find([complete_detect.coup] == 0), 2);
        All_subjects_table{s, 5} = size(find([complete_detect.coup] == 1), 2);
        All_subjects_table{s, 6} = size(find([complete_detect.coup] == 1), 2) / size(find(ismember({complete_detect.type}, 'Spindle')), 2);
        All_subjects_table{s, 7} = sum(coup_histogram(1:length(coup_histogram)/2)) / size(find(ismember({complete_detect.type}, 'Spindle')), 2);
        All_subjects_table{s, 8} = sum(coup_histogram((length(coup_histogram)/2 + 1):end)) / size(find(ismember({complete_detect.type}, 'Spindle')), 2);
        All_subjects_table{s, 9} = sum(coup_histogram(round(length(coup_histogram) / 2) - round(EEG.srate/2): round(length(coup_histogram) / 2))) / size(find(ismember({complete_detect.type}, 'Spindle')), 2);
        All_subjects_table{s, 10} = sum(coup_histogram(round(length(coup_histogram) / 2): round(length(coup_histogram) / 2) + round(EEG.srate/2))) / size(find(ismember({complete_detect.type}, 'Spindle')), 2);
        
        CSP_idx = find([complete_detect.coup] == 1);
        UNCSP_idx = find([complete_detect.coup] ~= 1);
        CSP = complete_detect(CSP_idx);
        UNCSP = complete_detect(UNCSP_idx);
        
        try
            All_subjects_csp(s).peakAmplitude = mean([CSP.peakAmplitude]);
        catch
            All_subjects_csp(s).peakAmplitude = mean([CSP.amplitude]);
        end
        All_subjects_csp(s).duration = mean([CSP.duration]) / EEG.srate;
        All_subjects_csp(s).density = size(find([complete_detect.coup] == 1), 2) / sleep_duration_total_min;
        try
            All_subjects_uncsp(s).peakAmplitude = mean([UNCSP.peakAmplitude]);
        catch
            All_subjects_uncsp(s).peakAmplitude = mean([UNCSP.amplitude]);
        end
        All_subjects_uncsp(s).duration = mean([UNCSP.duration]) / EEG.srate;
        All_subjects_uncsp(s).density = size(find([complete_detect.coup] == 0), 2) / sleep_duration_total_min;
        
        if ~isempty(coupled_SW)
            All_subjects_csw(s).number = size(coupled_SW, 2);
            All_subjects_csw(s).duration = mean([coupled_SW.duration]) / EEG.srate;
            All_subjects_csw(s).density = size(coupled_SW, 2) / sleep_duration_total_min;
        else
            All_subjects_csw(s).number = 0;
            All_subjects_csw(s).duration = 0;
            All_subjects_csw(s).density = 0;
        end
        
        All_subjects_uncsw(s).number = size(Uncoupled_SW, 2);
        All_subjects_uncsw(s).duration = mean([Uncoupled_SW.duration]) / EEG.srate;
        All_subjects_uncsw(s).density = size(Uncoupled_SW, 2) / sleep_duration_total_min;
        
        %% Generate plots (if selected)
        if stage_tag.plots == 1
        % Call function to create and save figures
        plot_SWSP_figures(savepath, fname_write, window_size, nBins, EEG, interval, rep_bins, rep_bins_sp_percent);
        end    
        if stage_tag.onsets == 1

%             extract_SWSP_onsets(savepath, fname_write, coupled_spindles, coupled_SW, Uncoupled_SW);
            extract_SWSP_onsets_TR(savepath, fname_write, coupled_spindles, coupled_SW, Uncoupled_SW, EEG, 1, 2.16, 1)
        end
        clearvars -except stage_tag savepath filenames s loadpath All_subjects_table All_subjects_csp All_subjects_uncsp all_rep_bins_SP all_rep_bins All_subjects_csw All_subjects_uncsw window_size nBins;
    catch
        continue
    end
end

finalize_SWSP_results(savepath, filenames, All_subjects_table, ...
    All_subjects_csp, All_subjects_uncsp, all_rep_bins, all_rep_bins_SP, ...
    All_subjects_csw, All_subjects_uncsw, window_size, nBins);
end