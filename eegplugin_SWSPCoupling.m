function vers = eegplugin_SWSPCoupling(fig, trystrs, catchstrs)
% EEGPLUGIN_SWSPCOUPLING - Adds SW-SP Coupling Analysis to EEGLAB Tools menu.
%
% Usage:
%   This function is automatically executed when EEGLAB starts,
%   adding a new menu item under "Tools".
%
% Inputs:
%   fig      - EEGLAB main GUI figure.
%   trystrs  - EEGLAB error handling strings.
%   catchstrs - EEGLAB error handling strings.
%
% Outputs:
%   vers - Plugin version.
%
% Author: Daniel Baena
% Version: 1.0

    vers = '1.0'; % Define plugin version
    if nargin < 3
        error('eegplugin_SWSPCoupling requires 3 input arguments.');
    end

    % Find the "Tools" menu in EEGLAB
    menu = findobj(fig, 'tag', 'tools'); 
    if isempty(menu)
        error('Could not find EEGLAB "Tools" menu.');
    end

    % Add "SW-SP Coupling" submenu under "Tools"
    uimenu(menu, 'Label', 'SW-SP Coupling Analysis', ...
        'Callback', 'EEG = pop_SWSPcoupling(EEG); eeglab redraw;', ...
        'Separator', 'on'); % Adds a separator above the new menu item

end
