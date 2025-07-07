function [EEG, com] = pop_SWSPcoupling(EEG)
% pop_SWSPcoupling - EEGLAB GUI wrapper for SW-SP coupling analysis.
%
% Description:
%   This function integrates the SW-SP coupling pipeline into EEGLAB’s GUI.
%   Users can select datasets, define analysis parameters, and launch processing
%   through EEGLAB's interface.
%
% Usage:
%   EEG = pop_SWSPcoupling(EEG);
%
% Parameters:
%   EEG - EEGLAB dataset structure.
%
% Outputs:
%   EEG - Updated EEGLAB dataset with additional `EEG.etc.SWSP` field storing results.
%
% Notes:
%   - This function serves as the main entry point for the EEGLAB plugin.
%   - The GUI parameters are passed to either `SWSPcoupling_fix` or `SWSPcoupling_adapt`.
%
% Author: Daniel Baena  
% Email: dbaenape@uottawa.ca  
% Affiliation: University of Ottawa  
% Date: 2025-02-06


com = ''; % Command string for EEGLAB history tracking

% Ensure EEG is provided
if nargin < 1 || isempty(EEG)
    error('EEG dataset is required. Load a dataset before running this function.');
end

% Ensure EEG structure has necessary fields
if ~isfield(EEG, 'event') || isempty(EEG.event)
    warning('The EEG dataset does not contain any events. SleepStage and event selection may not work.');
end

% Define GUI callbacks for interactive selection
cb_chan = [
    'EEG = get(gcbf, ''userdata''); ' ...
    'if isfield(EEG, ''chanlocs''); ' ...
    '   pop_chansel(EEG.chanlocs, ''field'', ''labels'', ''handle'', findobj(''parent'', gcbf, ''tag'', ''ChOI'')); ' ...
    'else; ' ...
    '   warndlg(''No channels available.'', ''Error''); ' ...
    'end'
];

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

% GUI layout
g = [3, 2];
geometry = {1, [2, 2, 1], [2, 2, 1], g, g, 1, [2, 2, 1], 1, 1};
geomvert = [1, 1, 1, 1, 1, 1, 1, 1, 1];

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
    {'style', 'checkbox', 'value', 1, 'string', 'Individual plots', 'tag', 'plots'}, ...
    {'style', 'checkbox', 'value', 1, 'string', 'Generate events onset mat files', 'tag', 'onsets'}, ...
    {'style', 'checkbox', 'value', 0, 'string', 'Use fixed time window (ignore event duration)', 'tag', 'fixed_window'} ...
};

% Launch GUI
result = inputgui('geometry', geometry, 'geomvert', geomvert, 'uilist', uilist, ...
    'title', 'SW-SP Coupling Analysis', 'userdata', EEG);

% Process user selections
if isempty(result)
    return; % User canceled
end

Inputarg.stages       = strsplit(result{1});
Inputarg.eventName    = result{2};
Inputarg.window_size  = str2double(result{3});
Inputarg.nBins        = str2double(result{4});
Inputarg.ChOI         = strsplit(result{5});
Inputarg.plots        = result{6};
Inputarg.onsets       = result{7};
Inputarg.fixed_window = result{8};

% Generate the custom save folder name
methodStr = "Fixed";
if Inputarg.fixed_window == 0
    methodStr = "Adapt";
end
stagesStr = strjoin(Inputarg.stages, "");
channelsStr = strjoin(Inputarg.ChOI, "");
savepath = fullfile(EEG.filepath, sprintf('%s_%s_%s_%s', methodStr, stagesStr, channelsStr, Inputarg.eventName));

% Create the save directory if it doesn't exist
if ~exist(savepath, 'dir')
    mkdir(savepath);
end

% Run the correct processing function
if Inputarg.fixed_window == 1
    [All_subjects_table] = SWSPcoupling_fix(EEG, EEG.filepath, {EEG.filename}, Inputarg, savepath);
else
    [All_subjects_table] = SWSPcoupling_adapt(EEG, EEG.filepath, {EEG.filename}, Inputarg, savepath);
end

% Save results inside EEG structure
EEG.etc.SWSP = struct();
EEG.etc.SWSP.results = 'Processed with SWSP Coupling';
EEG.etc.SWSP.summary_table = All_subjects_table;
EEG.etc.SWSP.savepath = savepath;  % Store save path inside EEG struct

% % Automatically open the results folder after processing
if ispc
    system(['explorer ', strrep(savepath, '/', '\')]);
elseif ismac
    system(['open ', savepath]);
elseif isunix
    system(['xdg-open ', savepath]);
end

% Generate EEGLAB command history
com = sprintf('EEG = pop_SWSPcoupling(EEG);');

% eeglab redraw;

end
