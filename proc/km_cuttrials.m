function [cfg,data] = km_cuttrials(cfg,data)
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

% FIXME: legacy reasons for app > apt
if isfield(data,'grip')
    if isfield(data.grip,'appaxis')
        data.grip.aptaxis = data.grip.appaxis;
        data.grip = rmfield(data.grip,'appaxis');
    end
    if isfield(data.grip,'app')
        data.grip.apt = data.grip.app;
        data.grip = rmfield(data.grip,'app');
    end
    if isfield(data.grip,'appvel')
        data.grip.aptvel = data.grip.appvel;
        data.grip = rmfield(data.grip,'appvel');
    end
    if isfield(data.grip,'appacc')
        data.grip.aptacc = data.grip.appacc;
        data.grip = rmfield(data.grip,'appacc');
    end
end

% get trial definitions and possibly variable names
if isfield(cfg,'trl')
    runtrl = cfg.trl;
    if isfield(cfg,'vars')
        vars = cfg.vars;
    end
elseif isfield(cfg,'trldef')
    if isfield(cfg.trldef,'trl')
        runtrl = cfg.trldef.trl;
        if isfield(cfg.trldef,'vars')
            vars = cfg.trldef.vars;
        end
    else
        cfg = km_trldef(cfg);
        runtrl = cfg.trl;
        if isfield(cfg,'vars')
            vars = cfg.vars;
        end
    end
else
    error('no trial definitions (trl) provided');
end
if ~iscell(runtrl), runtrl = {runtrl}; end

%% Distribute data over trials
%----------------------------------------
% initiate new data structure
newdata = initnewdata(data);

% initiate new configuration structure
[newcfg cfg] = initnewcfg(cfg);
    
% loop over runs
newtrl = [];
for r = 1:length(data.time)
    % get time
    tim = data.time{r};
    
    % get trial times
    runnr = data.run(r);
    trl = runtrl{runnr};
    
    % restrict to trials within time range
    idx = trl(:,1) >= min(tim) & trl(:,2) <= max(tim);
    ttrl = trl(idx,:);
    trl = time2idx(ttrl(:,1:2),tim);
    trl(:,3) = time2idx(ttrl(:,3),tim);
    
    % distribute data fields over trials
    tcfg	= [];
    tcfg.trl= trl;
    tcfg.r	= r;
    newdata = datarun2trl(tcfg,data,newdata);
    
    % distribute runnr and sampling frequency over trials
    runnrs = repmat(runnr,1,size(trl,1));
    fsamples = repmat(data.fsample(r),1,size(trl,1));
    if isfield(newdata,'run')
        newdata.run = [newdata.run runnrs];
        newdata.fsample = [newdata.fsample fsamples];
    else
        newdata.run = runnrs;
        newdata.fsample = fsamples;
    end
    
    % distribute configuration fields over trials
    tcfg        = [];
    tcfg.trl    = trl;
    tcfg.r      = r;
    tcfg.oldtime= data.time;
    tcfg.newtime= newdata.time;
    newcfg = cfgrun2trl(tcfg,cfg,newcfg);
    
    % construct new trl
    newtrl = [newtrl; ttrl];
    
end

% make sure that all trials have the same number of movements, fill with
% NaNs if needed
if isfield(cfg,'movement')
    cfg.movement = km_setequalnmov(cfg.movement);
end
if isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
    cfg.movdef.movement = km_setequalnmov(cfg.movdef.movement);
end
if isfield(data,'movement')
    data.movement = km_setequalnmov(data.movement);
end

% update trial definitions and timing if possible
% get indeces of time variables in logdata
if ~isfalse(cfg.cuttrials.timevar)
    if ~exist('vars','var')
        error('no variable names (vars) found, so the time can not be adjusted');
    elseif iscell(cfg.cuttrials.timevar)
        t_idx = ismember(vars,cfg.cuttrials.timevar);
    elseif ischar(cfg.cuttrials.timevar)
        t_idx = regexp(vars,cfg.cuttrials.timevar,'once');
        t_idx = cellfun(@(x) ~isempty(x),t_idx);
    else
        t_idx = cfg.cuttrials.timevar;
    end
    newtrl(:,t_idx) = newtrl(:,t_idx) - repmat(newtrl(:,3),1,sum(t_idx));
end

% replace old data and configuration structures with the new ones
newcfg.trl	= newtrl;
cfg         = newcfg;
data        = newdata;
%--------------------------------------------------------------------------


%% function initnewdata
%----------------------------------------
function newdata = initnewdata(data)

newdata = [];

% data fields from first level of data structure
field = {'label','axis','derivaxis','dimord'};
for f = 1:length(field)
    if isfield(data,field{f}),  newdata.(field{f}) = data.(field{f});   end
end

% data fields from grip level of data structure
if isfield(data,'grip')
    field = {'label','aptaxis','oriaxis','dimord'};
    for f = 1:length(field)
        if isfield(data.grip,field{f}),	newdata.grip.(field{f}) = data.grip.(field{f});	end
    end
end


%% function initnewcfg
%----------------------------------------
function [newcfg oldcfg] = initnewcfg(oldcfg)

newcfg = oldcfg;
if isfield(oldcfg,'artifact')
    newcfg.artifact = {};
elseif isfield(oldcfg,'artfctdef') && isfield(oldcfg.artfctdef,'artifact')
    newcfg.artifact = {};
    oldcfg.artfctdef = rmfield(oldcfg.artfctdef,'artifact');
end
if isfield(oldcfg,'movement')
    newcfg.movement = {};
elseif isfield(oldcfg,'movdef') && isfield(oldcfg.movdef,'movement')
    newcfg.movement = {};
    oldcfg.movdef = rmfield(oldcfg.movdef,'movement');
end


%% function datarun2trl
%----------------------------------------
function newdata = datarun2trl(cfg,olddata,newdata)

if nargin < 3,  newdata = [];   end

% trial number of new data
if isfield(newdata,'time')
    tinit = length(newdata.time);
else
    newdata.time = {};
    tinit = 0;
end

% run number of old data
rold = cfg.r;

% trial indeces
trl = cfg.trl;

% return quickly if nothing to do
if isempty(trl)
    return
end

% first devide the time array over trials and apply the offset
for t = 1:size(trl,1)
    % get trial time
    idx_beg	= trl(t,1);
    idx_end = trl(t,2);
    idx_off = trl(t,3);
    tim     = olddata.time{rold}(idx_beg:idx_end);
    offset	= olddata.time{rold}(idx_off);
    % store in new data
    tnew	= tinit + t;
    newdata.time{tnew} = tim - offset;
end

% data from first level of data structure
tcfg        = [];
tcfg.trl	= trl;
tcfg.rold	= rold;
tcfg.tinit	= tinit;
tcfg.field  = {'pos','vel','acc','ori','orivel','oriacc'};
newdata = dat2trl(tcfg,olddata,newdata);

% grip level of data structure
if isfield(olddata,'grip')
    if ~isfield(newdata,'grip'),    newdata.grip = struct;  end
    tcfg.field = {'apt','aptvel','aptacc','ori','orivel','oriacc'};
    newdata.grip = dat2trl(tcfg,olddata.grip,newdata.grip);
end

% distribute indeces fields over trials
tcfg        = [];
tcfg.trl	= trl;
tcfg.rold	= rold;
tcfg.tinit	= tinit;
tcfg.field  = {'event','artifact','movement'};
tcfg.oldtime= olddata.time;
tcfg.newtime= newdata.time;
newdata = idx2trl(tcfg,olddata,newdata);


%% function dat2trl
%----------------------------------------
function newdata = dat2trl(cfg,olddata,newdata)
for f = 1:length(cfg.field)
    if isfield(olddata,cfg.field{f})
        d = size(olddata.(cfg.field{f}){cfg.rold});
        for t = 1:size(cfg.trl,1)
            % get trial data
            idx_beg	= cfg.trl(t,1);
            idx_end	= cfg.trl(t,2);
            dat     = olddata.(cfg.field{f}){cfg.rold}(:,:,idx_beg:idx_end);
            % initialize new data field to remain size constancy
            tnew	= cfg.tinit + t;
            nsamp   = 1 + idx_end-idx_beg;
            newdata.(cfg.field{f}){tnew} = nan(d(1),d(2),nsamp);   
            % store in new data
            newdata.(cfg.field{f}){tnew}(:) = dat;
        end
    end
end


%% function idx2trl
%----------------------------------------
function newdata = idx2trl(cfg,olddata,newdata)
trl = cfg.oldtime{cfg.rold}(cfg.trl(:,1:3));
for f = 1:length(cfg.field)
    if isfield(olddata,cfg.field{f})
        oldidx = olddata.(cfg.field{f}){cfg.rold};
        for t = 1:size(trl,1)
            % get trial data
            if strcmpi(cfg.field{f},'event')
                i = oldidx(:,2) >= trl(t,1) & oldidx(:,2) <= trl(t,2);
            elseif strcmpi(cfg.field{f},'artifact')
                i = any(oldidx >= trl(t,1) & oldidx <= trl(t,2),2);
            else
                i = oldidx(:,1) >= trl(t,1) & oldidx(:,2) <= trl(t,2);
            end
            % store in new data
            tnew = cfg.tinit + t;
            newdata.(cfg.field{f}){tnew} = oldidx(i,:) - trl(t,3);
        end
    end
end

%% function cfgrun2trl
%----------------------------------------
function newcfg = cfgrun2trl(cfg,oldcfg,newcfg)

if nargin < 3,  newcfg = [];   end
tcfg        = [];
tcfg.field  = {'artifact','movement'};

% return quickly if nothing to do
if ~any(ismember(tcfg.field,fieldnames(newcfg))) ||...
        isempty(cfg.trl)
    return
end

% trial number of new cfg
for f = 1:length(tcfg.field)
    if isfield(newcfg,tcfg.field{f})
        tinit = length(newcfg.(tcfg.field{f}));
        break
    end
end

% distribute indeces fields over trials
tcfg.trl	= cfg.trl;
tcfg.rold	= cfg.r;
tcfg.tinit	= tinit;
tcfg.oldtime= cfg.oldtime;
tcfg.newtime= cfg.newtime;
newcfg = idx2trl(tcfg,oldcfg,newcfg);

