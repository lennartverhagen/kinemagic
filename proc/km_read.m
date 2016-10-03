function [cfg,data] = km_read(cfg)
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


%% Read data according to dataype
%----------------------------------------
% evaluate appropriate function
switch lower(cfg.read.type)
    case {'p','pol','polhemus'}
        [cfg,data] = km_read_polhemus(cfg);
    case {'m','mb','minibird'}
        [cfg,data] = km_read_minibird(cfg);
    case {'j','joy','joystick'}
        [cfg,data] = km_read_joystick(cfg);
    otherwise
        [cfg,data] = feval(['km_read_' cfg.read.type],cfg);
end

% evaluate appropriate experimentally specific function
if ~isfalse(cfg.read.expfun)
    [cfg,data] = feval(cfg.read.expfun,cfg,data);
end

% update dataset
cfg.dataset = '_RAW';


%% check crucial fields
%----------------------------------------
% set dimensions right
for r = 1:length(data.time)
    % time point
    [m,n] = size(data.time{r});
    if m > n,   data.time{r} = data.time{r}';   end
    
    % sample number
    if isfield(data,'samp')
        [m,n] = size(data.samp{r});
        if m > n,   data.samp{r} = data.samp{r}';   end
    end
    
    % event stamp
    if isfield(data,'event')
        [m,n] = size(data.event{r});
        if (m==1 || n==1) && m > n,	data.event{r} = data.event{r}';   end
    end
end

% cut data points where the time decreases (time-travel!)
for r = 1:length(data.time)
    % check for data quality
    idx_timetravel = diff(data.time{r})<=0;
    if sum(idx_timetravel) > 1000
        warning('KM:Read:TimingError','Too many samples have a negative time. It is very likely that there is a timing error in your data file. Please check it thoroughly.');
    elseif ~any(idx_timetravel)
        continue;
    end
    % find occurances of time travel
    tmin = find(idx_timetravel);
    % find all time points that lie in the past
    tpast = [];
    for i = 1:length(tmin)
        tmp = find(data.time{r}(tmin(i)+1:end)<=data.time{r}(tmin(i)));
        tpast = [tpast (tmin(i) + tmp)];
    end
    idx = zeros(size(data.time{r}));
    idx(tpast) = 1;
    % convolve to also delete the surrounding samples
    idx = logical(conv(idx,ones(1,5),'same'));
    % continue if all time points travel in the right direction in time
    if ~any(idx),   continue;   end
    
    % delete obnoxious data points
    if isfield(data,'samp') && all(size(data.time{r}) == size(data.samp{r}))
        data.samp{r} = data.samp{r}(~idx);
    end
    if isfield(data,'event') && all(size(data.time{r}) == size(data.event{r}))
        data.event{r} = data.event{r}(~idx);
    end
    if isfield(data,'pos')
        data.pos{r} = data.pos{r}(:,:,~idx);
    end
    if isfield(data,'ori')
        data.ori{r} = data.ori{r}(:,:,~idx);
    end
    data.time{r} = data.time{r}(~idx);
    
end

% convert event stamps to event codes using only the onset of an event
if isfield(data,'event')
    for r = 1:length(data.event)
        % if event is the direct read-out of a event channel
        if numel(data.event{r}) == numel(data.time{r})
            idx = data.event{r}>0;
            evts = unique(data.event{r}(idx));
            for e = 1:length(evts)
                % get event
                evt = data.event{r}==evts(e);
                % convolve to debounce the event channel
                evt = logical(conv(single(evt),[0 1 1],'same'));
                % convert to time domain (on- and off-sets)
                evt = km_logic2time(evt,data.time{r});
                % store only event onsets
                data.event{r} = [repmat(evts(e),size(evt,1),1) evt(:,1)];
            end
            % sort event in time
            if ~any(idx)
                data.event{r} = nan(0,2);
            else
                data.event{r} = sortrows(data.event{r},2);
            end
        end
    end
end

% dimension order
if ~isfield(data,'dimord')
    if ndims(data.pos{1}) == 3
        data.dimord = 'marker_axis_time';
    elseif ndims(data.pos{1}) == 2
        data.dimord = 'axis_time';
    elseif ndims(data.pos{1}) == 1
        data.dimord = 'time';
    else
        error('unknown dimension order');
    end
end

% run number
if ~isfield(data,'run')
    data.run = 1:length(data.time);
end

% sampling frequency
if ~isfield(data,'fsample')
    for r = 1:length(data.time)
        twin = max(data.time{r}) - min(data.time{r});
        data.fsample(r) = (length(data.time{r})-1)/twin;
    end
end

% marker labels
if ~isfield(data,'label')
    d = finddim(data.dimord,'marker');
    nmarker = size(data.pos{1},d);
    for m = 1:nmarker
        data.label{m} = sprintf('m%d',m);
    end
end

% spatial axes labels
aname = {'x','y','z'};
if ~isfield(data,'axis')
    d = finddim(data.dimord,'axis');
    naxes = size(data.pos{1},d);
    data.axis = aname(1:naxes);
else
    aname = aname(1:length(data.axis));
end

% permute the axes to fit the default x, y, z
if ~isequal(aname,data.axis) && isequal(sort(aname),sort(data.axis))
    [~,idx] = ismember(aname,data.axis);
    data.axis = data.axis(idx);
    for r = 1:length(data.pos)
        data.pos{r} = data.pos{r}(:,idx,:);
    end
    if isfield(data,'ori')
        for r = 1:length(data.ori)
            data.ori{r} = data.ori{r}(:,idx,:);
        end
    end
end
