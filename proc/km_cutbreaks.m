function [cfg,data] = km_cutbreaks(cfg,data)
%--------------------------------------------------------------------------
%
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

% find on and offsets of data chunks
idx_on = cell(1,length(data.time));
idx_off = cell(1,length(data.time));
for r = 1:length(data.time)
    
    % determine maximum number of missed samples (using data.samp) or
    % maximum duration of a break in seconds (using data.time)
    if ischar(cfg.cutbreaks.maxbreak)
        if isfield(data,'samp')
            dsamp = diff(data.samp{r});
        else
            dsamp = diff(data.time{r});
        end
        
        % step size in sample numbers or seconds
        samp = mode(dsamp);	
        try
            maxbreak = eval(cfg.cutbreaks.maxbreak);
        catch
            maxbreak = str2double(regexp(cfg.cutbreaks.maxbreak,'\d*','match','once')) * step;
        end
    else
        maxbreak = cfg.cutbreaks.maxbreak;
        if maxbreak < 1
            dsamp = diff(data.time{r});
        else
            if isfield(data,'samp')
                dsamp = diff(data.samp{r});
            else
                warning('KM:Configuration','You seem to be allowing breaks longer than 1 second. If you meant to specify a maximum number of samples instead of time in seconds, you should include a ''samp'' field in your data');
                dsamp = diff(data.time{r});
            end
        end
    end
    
    idx_on{r} = find([inf dsamp] > maxbreak);
    idx_off{r} = find([dsamp inf] > maxbreak);
end

% if more than one chunk is found in a run, split that run
nchunks = cellfun(@(x) length(x),idx_on);
if any(nchunks > 1)
    
    % add run index number to data
    if ~isfield(data,'run')
        data.run = 1:length(data.time);
    end
    
    % set desired dimension order
    data = km_setdimord(data,data.dimord,'marker_axis_time');
    
    % loop over old runs/repetitions
    newr = 1;
    olddata = data;
    for oldr = 1:length(olddata.time)
        
        % loop over chunks per run
        for c = 1:nchunks(oldr)
            idx = idx_on{oldr}(c):idx_off{oldr}(c);
            
            % check if the new chunk is long enough
            if ischar(cfg.cutbreaks.minchunk)
                minchunk = str2double(regexp(cfg.cutbreaks.minchunk,'\d*','match','once'));
                if sum(idx) < minchunk, continue;   end
            else
                minchunk = cfg.cutbreaks.minchunk;
                mintim = olddata.time{oldr}(idx(1));
                maxtim = olddata.time{oldr}(idx(end));
                if maxtim-mintim < minchunk, continue;   end
            end
            
            % assign correct run number
            data.run(newr) = olddata.run(oldr);
            
            % position
            if isfield(data,'pos')
                data.pos{newr} = olddata.pos{oldr}(:,:,idx);
            end
            
            % orientation
            if isfield(data,'ori')
                data.ori{newr} = olddata.ori{oldr}(:,:,idx);
            end
            
            % time
            data.time{newr} = olddata.time{oldr}(idx);
            
            % sample number
            if isfield(data,'samp')
                data.samp{newr} = olddata.samp{oldr}(idx);
            end
            
            % event
            if isfield(data,'event')
                mintim = min(data.time{newr});
                maxtim = max(data.time{newr});
                evtidx = olddata.event{oldr}(:,2) >= mintim & olddata.event{oldr}(:,2) <= maxtim;
                data.event{newr} = olddata.event{oldr}(evtidx,:);
            end
            
            % increase new run number
            newr = newr + 1;
            
        end
    end
    
    % set dimension order back
    data = km_setdimord(data,'marker_axis_time',data.dimord);
    
    % update sampling frequency
    if isfield(data,'fsample')
        for r = 1:length(data.time)
            twin = max(data.time{r}) - min(data.time{r});
            data.fsample(r) = (length(data.time{r})-1)/twin;
        end
    end
    
end