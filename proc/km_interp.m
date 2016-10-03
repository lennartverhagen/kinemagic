function [cfg,data] = km_interp(cfg,data)
%--------------------------------------------------------------------------
%
% FIXME
% Now the data in interpolated if data sampling is unstable above a 1e-9
% second accuracy. That might be a little strict. Maybe it should be based
% on the step size. Sometimes the time is rounded to the millisecond, than
% you could get a train of sampling durations of for example:
% [4 4 4 4 5 4 4 4 4 5 4 4 4 4 5] etc.
% In such cases interpolation is not neccesary, but the time points should
% be adjusted (slightly).
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


% return if no interpolation is requested
if isfalse(cfg.interp) || isfalse(cfg.interp.method)
    return
end

% if interpolation was not deemed necessary based on the time, return
cfg.interp = checkinterp(cfg.interp,data);
if isfalse(cfg.interp)
    return
end

% interpolate time (if needed)
data = interptime(cfg.interp,data);

% interpolate data (if needed)
data = interpdata(cfg.interp,data);

% remove the interpolated time field from the data structure
data.time = data.interptime;
data = rmfield(data,'interptime');
if isfield(data,'samp')
    data = rmfield(data,'samp');
end


%% function checkinterp
%----------------------------------------
function cfg = checkinterp(cfg,data)

% get time and loop over trials
check = false;
for r = 1:length(data.time)
    
    x = data.time{r};
    
    % determine time step
    diffx = diff(x);
    if isnumeric(cfg.freq)
        step = 1/cfg.freq;
    elseif any(strcmpi(cfg.freq,{'max','orig','task'}))
        step = mode(diffx);
    end
    
    % determine if interpolation is needed
    if max(abs(diffx-step)) > 1e-9
        check = true;
    end
end

% return if no interpolation is needed
if ~check,   cfg = 'no'; end


%% function interptime
%----------------------------------------
function data = interptime(cfg,data)

% get time and loop over trials
for r = 1:length(data.time)
    
    x = data.time{r};
    if length(x) < 4
        error('number of samples is too small for sensible interpolation');
    end
    
    % determine time step
    diffx = diff(x);
    if isnumeric(cfg.freq)
        step = 1/cfg.freq;
    elseif any(strcmpi(cfg.freq,{'max','orig','task'}))
        step = mode(diffx);
    end
    step = round(step*1e9)/1e9;
    
    % get new time array
    xi = min(x) : step : max(x);
    
    % store new time and sampling frequency
    data.interptime{r} = xi;
    data.fsample(r) = 1/step;
    data.fsample(r) = round(data.fsample(r)*1e9)/1e9;
end


%% function interpdata
%----------------------------------------
function data = interpdata(cfg,data)

param = intersect(fieldnames(data),{'pos','ori'});
for p = 1:length(param)
    dat = data.(param{p});
    for r = 1:length(dat)
        % original time
        x = data.time{r};
        % time to interpolate to
        xi = data.interptime{r}(:);
        % original data
        dims = size(dat{r});
        y = permute(dat{r},[3 2 1]);
        y = reshape(y,dims([3 2 1]));
        % interpolate
        yi = interp1(x,y,xi,cfg.method,'extrap');
        % shape back to original size
        dims(3) = size(yi,1);
        yi = permute(yi,[3 2 1]);
        dat{r} = reshape(yi,dims);
    end
    data.(param{p}) = dat;
end
