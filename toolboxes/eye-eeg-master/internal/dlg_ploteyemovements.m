% dlg_ploteyemovements - pops user dialogue, called by pop_ploteyemovements.m
%                see >> help pop_ploteyemovements
%
% Copyright (C) 2009-2020 Olaf Dimigen & Ulrich Reinacher, HU Berlin
% olaf.dimigen@hu-berlin.de 

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, 51 Franklin Street, Boston, MA 02110-1301, USA

function [sacstring, fixstring, metric] = dlg_ploteyemovements(callingFcn,EEG)

geometry = { 1 [2 1.3 0.4] [2 1.3 0.4] [2 1.3 0.4]};

%% callbacks
cbevent1 = ['if ~isfield(EEG.event, ''type'')' ...
    '   errordlg2(''No type field'');' ...
    'else' ...
    '   if isnumeric(EEG.event(1).type),' ...
    '        [tmps,tmpstr] = pop_chansel(unique([ EEG.event.type ]));' ...
    '   else,' ...
    '        [tmps,tmpstr] = pop_chansel(unique({ EEG.event.type }));' ...
    '   end;' ...
    '   if ~isempty(tmps)' ...
    '       set(findobj(''parent'', gcbf, ''tag'', ''sevent''), ''string'', tmpstr);' ...
    '   end;' ...
    'end;' ...
    'clear tmps tmpv tmpstr tmpfieldnames;' ];

cbevent2 = ['if ~isfield(EEG.event, ''type'')' ...
    '   errordlg2(''No type field'');' ...
    'else' ...
    '   if isnumeric(EEG.event(1).type),' ...
    '        [tmps,tmpstr] = pop_chansel(unique([ EEG.event.type ]));' ...
    '   else,' ...
    '        [tmps,tmpstr] = pop_chansel(unique({ EEG.event.type }));' ...
    '   end;' ...
    '   if ~isempty(tmps)' ...
    '       set(findobj(''parent'', gcbf, ''tag'', ''fevent''), ''string'', tmpstr);' ...
    '   end;' ...
    'end;' ...
    'clear tmps tmpv tmpstr tmpfieldnames;' ];

%% main menu
uilist = {...
    {'style', 'text', 'string', 'How are your saccade & fixation events called?', 'fontweight','bold'},...
    ...
    {'style', 'text', 'string', 'Saccade event:'},...
    {'style', 'edit', 'string', 'saccade', 'tag', 'sevent' },...
    {'style', 'pushbutton', 'string', '...', 'callback', cbevent1}, ...
    ...
    {'style', 'text', 'string', 'Fixation event:'},...
    {'style', 'edit', 'string', 'fixation', 'tag', 'fevent' }, ...
    {'style', 'pushbutton', 'string', '...', 'callback', cbevent2}, ...
    ...
    {'style', 'text', 'string', 'Unit for saccade amplitudes:'},...
    {'style', 'edit', 'string', 'degree', 'tag', 'metric'}, ...
    {},...
    };


%% make GUI
[results tmp tmp outstruct] = inputgui( 'geometry',geometry, ...
    'uilist',uilist,'helpcom', ['pophelp(''' callingFcn ''');'],...
    'title', ['Plot properties of saccades and/or fixations -- ', callingFcn]);


%% process user input (cancel)
if isempty(results)
    return
end

%% collect dialogue entries
sacstring = outstruct.sevent;
fixstring = outstruct.fevent;
metric    = outstruct.metric;

end