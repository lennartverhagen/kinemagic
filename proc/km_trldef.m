function cfg = km_trldef(cfg)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end


%% Define trials according to parameter
%----------------------------------------
% initialize new configuration structure
tcfg    = cfg.trldef;
switch lower(cfg.trldef.param)
    
    % evaluate experimentally specific function
    case 'exp'
        tcfg.exp	= cfg.exp;
        tcfg.subj	= cfg.subj;
        tcfg.sess	= cfg.sess;
        funname = ['km_trldef_' cfg.exp];
        trl = feval(funname,tcfg);
    
    % base trial definitions on the logfile
    case 'log'
        % get log
        if isfield(cfg,'log')
            tcfg.log = cfg.log;
        elseif isfield(cfg,'readlog')
            if isfield(cfg.readlog,'log')
                tcfg.log = cfg.readlog.log;
            else
                cfg = km_readlog(cfg);
                tcfg.log = cfg.log;
            end
        else
            error('no log provided');
        end
        
        [trl vars] = km_trldef_log(tcfg);
        
    % evaluate appropriate function
    otherwise
        error('parameter ''%s'' not recognized',cfg.trldef.param);

end

% catenate runs if requested
if istrue(cfg.trldef.catruns)
    trl = vertcat(trl{:});
end

% store trial definition
cfg.trl	= trl;
if exist('vars','var')
    cfg.vars = vars;
end
