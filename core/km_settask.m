function task = km_settask(task,method)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% check input
if nargin < 2,      method = 'free';    end
if ~iscell(task);   task = {task};      end

% loop over tasks
for t = 1:length(task)
    
    % remove any prepending 'km_'
    task{t} = regexprep(task{t},'km_','');
    
    switch lower(task{t})
        
        % preprocessing
        case {'r','read'}
            task{t} = 'read';
        case {'b','chunk','cutbreaks'}
            task{t} = 'cutbreaks';
        case {'l','log','readlog'}
            task{t} = 'readlog';
        case {'al','alignlog'}
            task{t} = 'alignlog';
        case {'h','hemi','hemifieldcross'}
            task{t} = 'hemifieldcross';
        case {'a','axes'}
            task{t} = 'axes';
        case {'c','calc','calcmarker'}
            task{t} = 'calcmarker';
        case {'i','interp','interpolate'}
            task{t} = 'interp';
        case {'art','artifact','artfctdef'}
            task{t} = 'artfctdef';
        case {'m','mask','maskart','maskartifact'}
            task{t} = 'maskart';
        case {'f','filt','filter'}
            task{t} = 'filter';
        case {'g','grp','grip'}
            task{t} = 'grip';
        case {'d','deriv','derivative','v','vel','velocity','acc','acceleration','velacc',}
            task{t} = 'deriv';
        
        % processing
        case {'md','movdef','movedef','movementdef','defmov','defmove','defmovement'}
            task{t} = 'movdef';
        case {'mc','movcheck','movecheck','movementcheck','checkmov','checkmove','checkmovement'}
            task{t} = 'movcheck';
        case {'movplot'}
            task{t} = 'movplot';
        case {'t','trldef','trialdef'}
            task{t} = 'trldef';
        case {'ct','cuttrl','cuttrials','dividetrials'}
            task{t} = 'cuttrials';
        case {'tc','trlcheck','trialcheck','checktrl','checktrial'}
            task{t} = 'trlcheck';
        case {'ms','movsel','movselect','movementsel','movementselect'};
            task{t} = 'movselect';
        case {'at','adjusttime','timeadjust'}
            task{t} = 'adjusttime';
        case {'mp','movparam','movementparam'}
            task{t} = 'movparam';
        case {'mpci','movparamci','movementparamci'}
            task{t} = 'movparamCI';
        case {'tr','trlrej','trlreject','trialrej','trialreject'}
            task{t} = 'trlrej';
        case {'err','eperr','eposerr','epointerr','endposerr','endpointerr','endpointerror'}
            task{t} = 'eposerr';
        case {'tra','traject','trajectory','trajectana','trajectoryana','trajectoryanalysis'}
            task{t} = 'traject';
            
        % reporting and plotting
        case 'paramreport'
            task{t} = 'paramreport';
        case 'paramplot'
            task{t} = 'paramplot';
        
        % not recognized
        otherwise
            if strcmpi(method,'strict')
                error('KM %s: Task ''%s'' is not supported.',mfilename,task{t});
            elseif strcmpi(method,'loose')
                warning('KM:NotSupported','%s: Task ''%s'' is not supported.',mfilename,task{t});
            end
    end
    
end

% check output
if length(task)==1,   task = task{1};  end
