function [cfg,data] = km_alignlog(cfg,data)
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
cfg = km_setcfg(cfg,task,'strict');

% return if requested
if isfalse(cfg.(task))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end


%% Match data runs with logfile runs
%----------------------------------------
datarun = unique(data.run);
logrun = unique(cfg.log.runnr);
dl = setdiff(datarun,logrun);
ld = setdiff(logrun,datarun);
if length(dl) == 1
    error('KM:MissingLogRun','the data contains a run [%s] that is not present in the logfile',num2str(dl));
elseif length(dl) > 1
    error('KM:MissingLogRun','the data contains runs [%s] that are not present in the logfile',num2str(dl));
end
if length(ld) == 1
	warning('KM:MissingDataRun','the logfile contains a run [%s] that is not present in the data',num2str(ld));
elseif length(ld) > 1
	warning('KM:MissingDataRun','the logfile contains runs [%s] that are not present in the data',num2str(ld));
end
% remove runs without data from the logfile
idx = ~ismember(cfg.log.runnr,ld);
cfg.log.runnr = cfg.log.runnr(idx);
cfg.log.data = cfg.log.data(idx);

    
%% Align logfile with data
%----------------------------------------
switch lower(cfg.alignlog.param)
    case 'no'
        return
        
    case 'event'
        
        % loop over runs
        runnr = -Inf;
        nevt = nan(1,length(data.time));
        offset = nan(1,length(data.time));
        for r = 1:length(data.time)
            
            % get events
            idx = data.event{r}(:,1) == cfg.alignlog.evtcode;
            evt = data.event{r}(idx,2);
            nevt(r) = sum(idx);
            if nevt(r) == 0 || isnan(nevt(r))
                warning('KM:NoEventsInRun','no events found in run %d of %s',r,cfg.subj);
            end
            
            % get trigger times as defined by the logfile
            if data.run(r) > runnr
                runnr = data.run(r);
                idxvar = strcmpi(cfg.log.vars,cfg.alignlog.logvar);
                idxrun = cfg.log.runnr == runnr;
                trg = cfg.log.data{idxrun}(:,idxvar);
            end
            
            % get time offset by aligning the trigger with the trial times
            [offset(r) trg] = km_aligntrgevt(cfg.alignlog,trg,evt);
            
        end
        
    otherwise
        error('parameter ''%s'' not recognized',cfg.alignlog.param);
end


%% weight the offsets found on different chunks of the same run
%----------------------------------------
% loop over runs
nevt(isnan(offset)) = 0;
for r = unique(data.run)
    
    % get run specific offsets
    lrunidx = cfg.log.runnr == r;
    drunidx = data.run == r;
    roffset = offset(drunidx);
    rnevt   = nevt(drunidx);
    
    % determine modus of offset
    tmp = [];
    for i = 1:length(roffset)
        tmp = [tmp repmat(round(roffset(i)/cfg.alignlog.err),1,rnevt(i))];
    end
    m = mode(tmp)*cfg.alignlog.err;
    
    % replace outliers with NaNs
    idx = roffset < m-2*cfg.alignlog.err | roffset > m+2*cfg.alignlog.err;
    roffset(idx) = NaN;
    rnevt(idx) = 0;
    
    % get weighted offset
    s = nansum(roffset .* rnevt.^1.3);
    n = nansum(rnevt.^1.3);
    tmp = s/n;
    
    % the offset must be a finite number
    if ~isfinite(tmp)
        cfg.log.data(lrunidx) = [];
        cfg.log.runnr(lrunidx) = [];
        offset(drunidx) = NaN;
        nevt(drunidx) = 0;
    end
    
    % update offset
    offset(drunidx) = tmp;
    
end


%% adjust the timing
%----------------------------------------
switch lower(cfg.alignlog.adjust)
    
    % adjust the timing of the log
    case 'log'
        
        % find time variables in the logdata
        if iscell(cfg.alignlog.timevar)
            t_idx = ismember(cfg.log.vars,cfg.alignlog.timevar);
        elseif ischar(cfg.alignlog.timevar)
            t_idx = regexp(cfg.log.vars,cfg.alignlog.timevar,'once');
            t_idx = cellfun(@(x) ~isempty(x),t_idx);
        else
            t_idx = cfg.alignlog.timevar;
        end
        
        % loop over logruns
        for r = 1:length(cfg.log.data)
            
            % get run offset
            drunidx = data.run == cfg.log.runnr(r);
            tmp = offset(find(drunidx,1,'first'));
            
            % apply offset to the logdata
            cfg.log.data{r}(:,t_idx) = cfg.log.data{r}(:,t_idx) - tmp;
            
        end
        
        % adjust the timing of the data
    case 'data'
        
        % loop over dataruns
        for r = 1:length(cfg.log.data)
        
            % get run offset
            tmp = offset(r);
            
            % make offset concurrent with sampling rate
            tmp = round(tmp*data.fsample(r))/data.fsample(r);
        
            % apply offset to time fields in cfg and data structures
            [cfg data] = cfgdataoffset(cfg,data,tmp,r);
        end
        
    otherwise
        error('parameter ''%s'' not recognized',cfg.alignlog.adjust);
end
%--------------------------------------------------------------------------


%% function cfgdataoffset
%----------------------------------------
function [cfg data] = cfgdataoffset(cfg,data,offset,r)

if isempty(offset), offset = nan;   end

% apply offset to time fields in the data structure
data.time{r} = data.time{r} + offset;
if isfield(data,'event')
    data.event{r}(:,2) = data.event{r}(:,2) + offset;
end
if isfield(data,'artifact')
    data.artifact{r} = data.artifact{r} + offset;
end
if isfield(data,'movement')
    data.movement{r} = data.movement{r} + offset;
end

% apply offset to time fields in the configuration structure
if isfield(cfg,'artifact')
    cfg.artifact{r} = cfg.artifact{r} + offset;
end
if isfield(cfg,'artfctdef') && isfield(cfg.artfctdef,'artifact')
    cfg.artfctdef.artifact{r} = cfg.artfctdef.artifact{r} + offset;
end
if isfield(cfg,'movement')
    cfg.movement{r} = cfg.movement{r} + offset;
end
if isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
    cfg.movdef.movement{r} = cfg.movdef.movement{r} + offset;
end
