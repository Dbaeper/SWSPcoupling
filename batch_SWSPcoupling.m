% batch_SWSPcoupling - Batch processing script for SW-SP coupling analysis.
%
% Description:
%   This script runs **slow wave-spindle (SW-SP) coupling detection** on multiple EEG datasets.
%   Users select files, specify analysis parameters via a GUI, and the script processes
%   each dataset using either:
%     - A **fixed time window** (`SWSPcoupling_fix.m`)
%     - An **adaptive time window** (`SWSPcoupling_adapt.m`)
%
% Usage:
%   Run this script in MATLAB. The user will:
%   1. Select EEG datasets (.set files)
%   2. Choose a save folder
%   3. Define parameters in a GUI
%   4. Process files in batch mode
%
% Parameters (Selected via GUI):
%   - Sleep Stages (`stages`): Sleep stages to analyze (e.g., NREM2, NREM3).
%   - Slow Wave Label (`eventName`): Label for SW negative peak in `EEG.event`.
%   - Window Length (`window_size`): Time (seconds) around SW peak to analyze.
%   - Histogram Bins (`nBins`): Number of bins for event histograms.
%   - Channels of Interest (`ChOI`): EEG channels for analysis.
%   - Generate Plots (`plots`): If enabled, plots will be saved.
%   - Save Event Onsets (`onsets`): If enabled, event onset files will be stored.
%   - Fixed Window (`fixed_window`): If checked, a fixed time window is used; if unchecked, an adaptive window is applied.
%
% Outputs:
%   - Processed EEG datasets with results saved in `EEG.etc.SWSP`
%   - Individual `.mat` files with spindle detection results
%   - Summary `.xlsx` report containing all subjects' data
%   - Figures (`.tiff` and `.png`) showing SW-SP event distributions
%
% Notes:
%   - The EEG datasets must contain event markers (`EEG.event`) with labeled sleep stages.
%   - The analysis parameters are set once in the GUI and applied to all selected datasets.
%
% Dependencies:
%   - EEGLAB (Tested on EEGLAB 2021+)
%   - `SWSPcoupling_fix.m` (fixed window analysis)
%   - `SWSPcoupling_adapt.m` (adaptive window analysis)
%   - `extract_SWSP_onsets.m` (onset extraction)
%   - `plot_SWSP_figures.m` (figure generation)
%   - `finalize_SWSP_results.m` (summary compilation)
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

eeglab nogui % Initialize EEGlab

disp('Please select file(s) to process.');
[filenames, loadpath] = uigetfile( ...
    {'*.set', 'EEGlab Dataset (*.set)'}, 'multiselect', 'on');
if isequal(filenames, 0)
    disp('File selection canceled. Exiting...');
    return;
end
if ischar(filenames) % Only one file was selected
    filenames = cellstr(filenames); % Convert to cell array for consistency
end

timestamp = datestr(now, 'yyyy-mm-dd_HHMM'); % e.g., 2025-07-09_2315
savepath = [loadpath 'SWSPcouplingBatch_' timestamp filesep];
if ~isdir(savepath);mkdir(savepath);end

disp(['Data will be saved at ' savepath]);

% Load EEG dataset before GUI
EEG = pop_loadset([loadpath, filenames{1}]);

% Ensure EEG structure has necessary fields
if ~isfield(EEG, 'event') || isempty(EEG.event)
    warning('The EEG dataset does not contain any events. Sleep Stage and event selection may not work.');
end

% GUI geometry
g = [3, 2];
geometry = {1, [2, 2, 1], [2, 2, 1], g, g, 1, [2, 2, 1], 1, 1};
geomvert = [1, 1, 1, 1, 1, 1, 1, 1, 1];

% Callback for selecting available channels
cb_chan = [
    'EEG = get(gcbf, ''userdata''); ' ...
    'if isfield(EEG, ''chanlocs''); ' ...
    '   pop_chansel(EEG.chanlocs, ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''ChOI'')); ' ...
    'else; ' ...
    '   warndlg(''No channels available.'', ''Error''); ' ...
    'end'
];

% Callback for selecting unique SleepStage values
cb_stage = [
    'EEG = get(gcbf, ''userdata''); ' ...
    'if isfield(EEG, ''event'') && isfield(EEG.event, ''SleepStage''); ' ...
    '   sleepStages = unique({EEG.event.SleepStage}); ' ...
    '   [selection, ok] = listdlg(''PromptString'', ''Select Sleep Stages:'', ''ListString'', sleepStages, ''SelectionMode'', ''multiple''); ' ...
    '   if ok && ~isempty(selection); ' ...
    '       selectedStagesStr = strjoin(sleepStages(selection), '' ''); ' ...
    '       set(findobj(''parent'', gcbf, ''tag'', ''stages''), ''string'', selectedStagesStr); ' ...
    '   end; ' ...
    'else; ' ...
    '   warndlg(''No SleepStage field found in EEG.event.'', ''Error''); ' ...
    'end'
];


% Callback for selecting unique event types
cb_swtag = [
    'EEG = get(gcbf, ''userdata''); ' ...
    'if isfield(EEG, ''event'') && isfield(EEG.event, ''type''); ' ...
    '   eventTypes = unique({EEG.event.type}); ' ...
    '   [selection, ok] = listdlg(''PromptString'', ''Select Event Type:'', ''ListString'', eventTypes); ' ...
    '   if ok && ~isempty(selection); ' ...
    '       set(findobj(''parent'', gcbf, ''tag'', ''eventName''), ''string'', eventTypes{selection}); ' ...
    '   end; ' ...
    'else; ' ...
    '   warndlg(''No event types found in EEG.event.'', ''Error''); ' ...
    'end'
];


% Build the GUI with the new fixed-window checkbox at the end.
uilist = { ...
    {'style', 'text', 'string', 'Parameters for SW-SP coupling detection'}, ...
    {'style', 'text', 'string', 'Label Sleep Stages to include'}, ...
    {'style', 'edit', 'string', 'NREM2 NREM3', 'tag', 'stages'}, ...
    {'style', 'pushbutton', 'string', '...', 'callback', cb_stage}, ...
    {'style', 'text', 'string', 'Label for Slow Wave negative peak'}, ...
    {'style', 'edit', 'string', 'SWneg', 'tag', 'eventName'}, ...
    {'style', 'pushbutton', 'string', '...', 'callback', cb_swtag}, ...
    {'style', 'text', 'string', 'Window length in seconds to build around the SW negative peak'}, ...
    {'style', 'edit', 'string', '2', 'tag', 'window_size'}, ...
    {'style', 'text', 'string', 'Number of time bins for histograms'}, ...
    {'style', 'edit', 'string', '20', 'tag', 'nBins'}, ...
    {'style', 'text', 'string', 'Channel labels or indices'}, ...
    {'style', 'edit', 'string', 'Fz Cz Pz', 'tag', 'ChOI'}, ...
    {'style', 'pushbutton', 'string', '...', 'callback', cb_chan}, ...
    {'style', 'checkbox', 'value', 1, 'string', 'Generate individual plots', 'tag', 'plots'}, ...
    {'style', 'checkbox', 'value', 1, 'string', 'Generate events onset mat files', 'tag', 'onsets'}, ...
    {'style', 'checkbox', 'value', 0, 'string', 'Use fixed time window (ignore event duration)', 'tag', 'fixed_window'} ...
};

% Launch the GUI with the full EEG dataset in userdata
result = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, ...
    'title', 'SW-SP Coupling analyses', 'helpcom', 'pophelp(''SWSPcoupling'')', ...
    'userdata', EEG);
if isempty(result)
    disp('Parameter selection canceled. Exiting...');
    return;
end

% Process user selections
if ~isempty(result)
    Inputarg.stages       = strsplit(result{1});
    Inputarg.eventName    = result{2};
    Inputarg.window_size  = str2double(result{3});
    Inputarg.nBins        = str2double(result{4});
    Inputarg.ChOI         = strsplit(result{5});
    Inputarg.plots        = result{6};
    Inputarg.onsets       = result{7};
    Inputarg.fixed_window = result{8};
else
    return;
end

% Process each file
for s = 1:size(filenames, 1)
    EEG = pop_loadset([loadpath, filenames{s}]);
    
    if ~isempty(result)
        Inputarg.stages       = strsplit(result{1});
        Inputarg.eventName    = result{2};
        Inputarg.window_size  = str2double(result{3});
        Inputarg.nBins        = str2double(result{4});
        Inputarg.ChOI         = strsplit(result{5});
        Inputarg.plots        = result{6};
        Inputarg.onsets       = result{7};
        Inputarg.fixed_window = result{8};
    else
        return;
    end
    
    % Run the correct processing function
    if Inputarg.fixed_window == 1
        [All_subjects_table] = SWSPcoupling_fix(EEG, loadpath, filenames, Inputarg, savepath);
    else
        [All_subjects_table] = SWSPcoupling_adapt(EEG, loadpath, filenames, Inputarg, savepath);
    end
end
