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

    vers = '1.0'; % Initial release
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
