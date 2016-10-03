function log = km_readlogfun_GraspingSemantics(cfg)
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
cfg.expr.startexp	= 'Time of start Experiment:';
cfg.expr.stopexp	= 'eof';
cfg.expr.startrun	= 'Run';
cfg.expr.stoprun	= '\r\n\r\n';
cfg.expr.runname    = 'Run\s(\d+)\s';
cfg.expr.trigger    = 'run start\s*\S*:\s(\d+)\s';
cfg.expr.vars       = '===\s+(.*)$';
cfg.expr.timevar    = '^t_';
cfg.timscale        = 1e-3;

% check configuration
cfg = km_setcfg(cfg,{'subj','dirlog'});


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

% change the coding of the variable 'word' for subjects [1 2 3 99]
% old - 1:GROOT, 2:KLEIN, 3:GEEL, 4:BLAUW
% new - 1:KLEIN, 2:GROOT, 3:GEEL, 4:BLAUW
% you need the subject number to know in which category the subject belongs
subjnr = str2double(regexp(cfg.subj,'\d*','match','once'));    
if ismember(subjnr,[1 2 3 99])
    idx = find(strcmpi(vars,'word'));
    for r = 1:length(logdata)
        oldword = logdata{r}(:,idx);
        newword = oldword;
        newword(oldword==1) = 2;
        newword(oldword==2) = 1;
        logdata{r}(:,idx) = newword;
    end
end

% insert first and second half coding
idx = find(strcmpi(vars,'run'));
vars = [vars(1:idx) {'half'} vars(idx+1:end)];
for r = 1:length(logdata)
    run = logdata{r}(:,idx);
    tmp = zeros(size(run,1),1);
    % first half
    tmp(run<3) = 1;
    % second half
    tmp(run>2) = 2;
    logdata{r} = [logdata{r}(:,1:idx) tmp logdata{r}(:,idx+1:end)];
end

% insert size congruency
idx = find(strcmpi(vars,'word'));
i_size = strcmpi(vars,'size');
i_word = strcmpi(vars,'word');
vars = [vars(1:idx) {'sizecongr'} vars(idx+1:end)];
for r = 1:length(logdata)
    saiz = logdata{r}(:,i_size);
    word = logdata{r}(:,i_word);
    tmp = zeros(size(saiz,1),1);
    % small & KLEIN
    tmp(saiz==1 & word==1) = 1;
    % large & GROOT
    tmp(saiz==2 & word==2) = 1;
    % small & GROOT
    tmp(saiz==1 & word==2) = 2;
    % large & KLEIN
    tmp(saiz==2 & word==1) = 2;
    logdata{r} = [logdata{r}(:,1:idx) tmp logdata{r}(:,idx+1:end)];
end

% insert color congruency
idx = find(strcmpi(vars,'sizecongr'));
i_color = strcmpi(vars,'color');
i_word = strcmpi(vars,'word');
vars = [vars(1:idx) {'colorcongr'} vars(idx+1:end)];
for r = 1:length(logdata)
    color = logdata{r}(:,i_color);
    word = logdata{r}(:,i_word);
    tmp = zeros(size(color,1),1);
    % yellow & GEEL
    tmp(color==1 & word==3) = 1;
    % blue & BLAUW
    tmp(color==2 & word==4) = 1;
    % yellow & BLAUW
    tmp(color==1 & word==4) = 2;
    % blue & GEEL
    tmp(color==2 & word==3) = 2;
    logdata{r} = [logdata{r}(:,1:idx) tmp logdata{r}(:,idx+1:end)];
end

% insert bino-mono
idx = find(strcmpi(vars,'vision'));
vars = [vars(1:idx) {'binomono'} vars(idx+1:end)];
for r = 1:length(logdata)
    vision = logdata{r}(:,idx);
    tmp = zeros(size(vision,1),1);
    % bino
    tmp(vision==1) = 1;
    % mono
    tmp(vision>1)  = 2;
    logdata{r} = [logdata{r}(:,1:idx) tmp logdata{r}(:,idx+1:end)];
end

% store information in log structure
log.runname     = runname;
log.runnr       = cellfun(@str2num,runname);
log.vars        = vars;
log.data        = logdata;


