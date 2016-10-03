function log = km_readlogfun_kinematicTMSEEG(cfg)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------


%% Configuration
%----------------------------------------
% set configuration
cfg.expr.startexp	= 'Time of experiment start:';
cfg.expr.stopexp	= 'eof';
cfg.expr.startrun	= 'Session';
cfg.expr.stoprun	= '\r\n\r\n';
cfg.expr.runname    = 'Session\s*\S*:\s(\S+)\s';
cfg.expr.trigger    = 'session start\s*\S*:\s(\d+)\s';
cfg.expr.vars       = '===\s+(.*)$';
cfg.expr.timevar    = '^t_';
cfg.timscale        = 1e-3;

% check configuration
cfg = km_setcfg(cfg,{'subj','sess','dirlog'});


%% Read logfile
%----------------------------------------
% find logfile
fname = dir(fullfile(cfg.dir.log,['log_' cfg.subj '.txt']));

% read the logfile, based on start and stop expressions
tcfg                = [];
tcfg.logfile        = fullfile(cfg.dir.log,fname.name);
tcfg.expr           = cfg.expr;
[hdr,logdata] = km_readlog_text(tcfg);

% loop over runs, find triggers and run names
runname = cell(1,length(hdr));
trigger = nan(1,length(hdr));
for r = 1:length(hdr)
    tok = regexp(hdr{r},cfg.expr.runname,'tokens','once');
    runname{r} = tok{1};
    tok = regexp(hdr{r},cfg.expr.trigger,'tokens','once');
    trigger(r) = str2double(tok{1});
end

% capture variable names
tok = regexp(hdr{1},cfg.expr.vars,'tokens','once');
vars = regexp(tok{1},'\S*','match');


%% Adapt logfile to fit standards
%----------------------------------------
% correct timing based on the triggers and scale times from ms to s
t_idx = regexp(vars,cfg.expr.timevar,'once');
t_idx = cellfun(@(x) ~isempty(x),t_idx);
for r = 1:length(logdata)
    logdata{r}(:,t_idx) = (logdata{r}(:,t_idx) - trigger(r))*cfg.timscale;
end

% change name of variable 'log' to 'site'
vars{ismember(vars,'loc')} = 'site';

% change the coding of the variable 'site'
% old - 1:V6A, 2:AIP, 3:vertex
% new - 1:vertex, 2:V6A, 3:AIP 
idx = find(strcmpi(vars,'site'));
for r = 1:length(logdata)
    oldsite = logdata{r}(:,idx);
    newsite = zeros(size(oldsite,1),1);
    newsite(oldsite==1) = 2;
    newsite(oldsite==2) = 3;
    newsite(oldsite==3) = 1;
    logdata{r}(:,idx) = newsite;
end

% insert run number
idx = find(strcmpi(vars,'sess'));
vars = [vars(1:idx) {'run'} vars(idx+1:end)];
runnr = nan(1,length(logdata));
for r = 1:length(logdata)
    ntrl = size(logdata{r},1);
    runnr(r) = str2double(runname{r}(end));
    tmp = repmat(runnr(r),ntrl,1);
    logdata{r} = [logdata{r}(:,1:idx) tmp logdata{r}(:,idx+1:end)];
end

% insert early late TMS
idx_tTMS = strcmpi(vars,'TTMS');
idx = find(strcmpi(vars,'site'));
vars = [vars(1:idx) {'TMSel'} vars(idx+1:end)];
for r = 1:length(logdata)
    tTMS = logdata{r}(:,idx_tTMS);
    tmp = zeros(size(tTMS,1),1);
    tmp(tTMS>50 & tTMS<250) = 1;
    tmp(tTMS>250 & tTMS<450) = 2;
    logdata{r} = [logdata{r}(:,1:idx) tmp logdata{r}(:,idx+1:end)];
end
    
% insert verthori
idx = find(strcmpi(vars,'slant'));
vars = [vars(1:idx) {'verthori'} vars(idx+1:end)];
for r = 1:length(logdata)
    slant = logdata{r}(:,idx);
    tmp = zeros(size(slant,1),1);
    tmp(slant<45) = 1;
    tmp(slant>45) = 2;
    logdata{r} = [logdata{r}(:,1:idx) tmp logdata{r}(:,idx+1:end)];
end

% get session names
sessname = cellfun(@(x) x(1:end-1),runname,'UniformOutput',false);

% store information in log structure
log.runname     = runname;
log.sessname    = sessname;
log.runnr       = runnr;
log.vars        = vars;
log.data        = logdata;


