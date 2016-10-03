function [offset trg trgchk evtchk] = km_aligntrgevt(cfg,trg,evt)
%--------------------------------------------------------------------------
%
% See also KM_ALIGNLOG, KM_CUTTRIALS
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% set configuration
if ~isfield(cfg,'initidx'),     cfg.initidx     = 1:2;     	end
if ~isfield(cfg,'err'),         cfg.err         = 0.01;     end


% settings
initidx	= cfg.initidx;
err     = cfg.err;

% return quickly if nothing to do
if length(trg)<length(initidx) || length(evt)<length(initidx)
    warning('KM:cuttrials','not enough trials in this run to align the logfile with the data');
    offset = nan;
    return
end

% check if trial timing is jittered or constant
trgstep = median(diff(trg));
if mean(abs(diff(trg) - trgstep) < err) > 0.5
    warning('KM:cuttrials','The trial trigger is not jittered, this can be detrimental for aligning the data with the trigger');
end

% get initial estimation of the offset, this only works properly if you
% have jittered trigger timing
offset = [];
trgidx = initidx;
while isempty(offset) && max(trgidx) <= length(trg)
    evtidx = initidx;
    while isempty(offset) && max(evtidx) <= length(evt)
        offset = trg(trgidx)-evt(evtidx);
        if any(abs(offset - mean(offset)) > err)
            offset = [];
            evtidx = evtidx + 1;
        end
    end
    trgidx = trgidx + 1;
end
offset = mean(offset);

% check if an offset could be determined
if isempty(offset) || isnan(offset)
    warning('KM:cuttrials','no initial offset could be detected to align the data with the trigger');
    offset = nan;
    return
end

% using the trigger times as the gold standard, try to match the triggers
% with the events
trgchk = false(1,length(trg));
evtchk = false(1,length(evt));
t = 1;
e = 1;
while t <= length(trg) && e <= length(evt)
    % calculate new offset
    newoff = trg(t)-evt(e);
    
    % check if it falls within the error margin
    if abs(mean(offset)-newoff) <= err
        trgchk(t) = true;
        evtchk(e) = true;
        % update offset
        if t==1 || e==1
            offset = newoff;
        else
            offset = [offset newoff];
        end
        % update trigger and event indeces
        t = t + 1;
        e = e + 1;
    elseif newoff < mean(offset)-err
        t = t + 1;
    elseif newoff > mean(offset)+err
        e = e + 1;
    else
        % if mean(offset) results in a nan, skip one step ahead
        t = t + 1;
        e = e + 1;
    end
end

% calculate the offset between the selected triggers and events
offset = trg(trgchk) - evt(evtchk);
% cut outliers
offset = offset(~cutoutlier(offset,2,'iqr'));
% calculate the mean
offset = mean(offset);

% keep only unused triggers
idx = find(trgchk,1,'last');
trg = trg(idx+1:end);

% if no offset could be found, make it not a number to ensure that no
% triggers will be found
if isempty(offset) || isnan(offset)
    warning('KM:cuttrials','no offset could be found to align the data with the trigger');
    offset = nan;
end