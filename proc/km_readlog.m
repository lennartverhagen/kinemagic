function cfg = km_readlog(cfg)
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


%% Read log according to experiment
%----------------------------------------
% evaluate appropriate function
tcfg        = cfg.readlog;
tcfg.subj	= cfg.subj;
tcfg.sess	= cfg.sess;
tcfg.dir	= cfg.dir;
log = feval(cfg.readlog.expfun,tcfg);

% keep only logs of current session
if isfield(log,'sessname')
    sess = cfg.sess;
    % remove first low dash '_' from session name
    if ~isempty(sess) && strcmpi(sess(1),'_')
        sess = sess(2:end);
    end
    idx = strcmpi(log.sessname,sess);
    log.data        = log.data(idx);
    if isfield(log,'runnr')
        log.runnr	= log.runnr(idx);
    end
    % remove redundant fields
    log	= rmfield(log,'sessname');
    if isfield(log,'runname')
        log = rmfield(log,'runname');
    end
end

% check logdata and vars
ndatacol = cellfun(@(x) size(x,2),log.data);
if length(unique([ndatacol length(log.vars)])) ~= 1
    error('the number of log data columns[%s] does not match the number of log variables [%d]',num2str(ndatacol),length(log.vars));
end

% store logfile
cfg.log	= log;




