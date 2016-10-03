function [cfg,data] = km_adjusttime(cfg,data)
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


%% Non-continuous data
%----------------------------------------
% concatenate trl matrix
trl = vertcat(runtrl{:});
% get time window of data trials
timwin = nan(length(data.time),2);
timwin(:,1) = cellfun(@min,data.time);
timwin(:,2) = cellfun(@max,data.time);

% check if the time windows of the data trials do not overlap
if any(timwin(1:end-1,2) > timwin(2:end,1))
    error('time windows of the data.time{:} are overlapping. This does not allow a correct adjustment.')
end

% loop over the data time windows
idx_trl = zeros(size(timwin,1),1);
for r = 1:size(timwin,1)
    idx = timwin(r,1) <= trl(:,3) & trl(:,3) <= timwin(r,2);
    if any(idx)
        idx_trl(r) = find(idx,1,'first');
    end
end

% reject data runs that do not match the trl matrix
idx_dat = idx_trl>0;
data.time       = data.time(idx_dat);
data.pos        = data.pos(idx_dat);
data.fsample    = data.fsample(idx_dat);
data.run        = data.run(idx_dat);
data.artifact   = data.artifact(idx_dat);
data.vel        = data.vel(idx_dat);
data.movement   = data.movement(idx_dat);
timwin          = timwin(idx_dat,:);

% reject data run fields from the cfg structure
cfg.movement    = cfg.movement(idx_dat);
cfg.artifact    = cfg.artifact(idx_dat);

% reject rows from trl matrix that do not match the data
idx_trl = idx_trl(idx_dat);
trl = trl(idx_trl,:);

% override the on- and off-sets of the trl matrix
% FIXME: this is not strictly necessary, but it is the easiest way
% I should be using nearest, but then I would have to update the indeces
% arrays of data.event, data.movement and data.artifact
trl(:,1:2) = timwin;

% loop over runs and adjust the time
for r = 1:length(data.time)
    offset = trl(r,3);
    data.time{r} = data.time{r} - offset;
    data.movement{r} = data.movement{r} - offset;
    cfg.movement{r} = cfg.movement{r} - offset;
    cfg.artifact{r} = cfg.artifact{r} - offset;
end


% update trial definitions and timing if possible
% get indeces of time variables in logdata
if ~isfalse(cfg.adjusttime.timevar)
    if ~exist('vars','var')
        error('no variable names (vars) found, so the time can not be adjusted');
    elseif iscell(cfg.adjusttime.timevar)
        t_idx = ismember(vars,cfg.adjusttime.timevar);
    elseif ischar(cfg.adjusttime.timevar)
        t_idx = regexp(vars,cfg.adjusttime.timevar,'once');
        t_idx = cellfun(@(x) ~isempty(x),t_idx);
    else
        t_idx = cfg.adjusttime.timevar;
    end
    trl(:,t_idx) = trl(:,t_idx) - repmat(trl(:,3),1,sum(t_idx));
end

% store the trl matrix in the configuration
cfg.trl = trl;
%--------------------------------------------------------------------------