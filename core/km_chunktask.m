function [cfg,data,rtrnflg] = km_chunktask(cfg,data)
%--------------------------------------------------------------------------
%
% See also KM_DOTASK
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% let's see what there is to do
rtrnflg = false;
if isfield(cfg,'task') && ~isempty(cfg.task)
    task = cfg.task;
else
    return
end

% get the correct naming conventions for the different tasks
task = km_settask(task,'strict');
if ~iscell(task);   task = {task};	end

% find unique tasks
b = unique(task);

% check if a task in b does not want to be chunked
tmp_b = {};
for i = 1:length(b)
    if isfield(cfg.(b{i}),'chunk') && isfalse(cfg.(b{i}).chunk)
        tmp = cellfun(@(n) sprintf('%s_%i',b{i},n),{1,2,3},'UniformOutput',false);
        tmp_b = [tmp_b tmp];
        task(ismember(task,b{i})) = tmp;
    else
        tmp_b = [tmp_b b{i}];
    end
end
b = tmp_b;

% return if no chunking is needed
if length(task) == length(b);
    return
end

% loop over unique elements and number the repetitions
c = zeros(size(task));
for i = 1:length(b)
    idx = ismember(task,b{i});
    if sum(idx) > 1
        c(idx) = 1:sum(idx);
    end
end

% identify task chunks
chunknr = 1;
i = 1;
j = 2;
d = zeros(size(task));
while j <= length(task)
    if length(unique(task(i:j))) ~= j-i+1
        d(i:j-1) = chunknr;
        chunknr = chunknr + 1;
        i = j;
    end
    j = j + 1;
end
d(i:end) = chunknr;

% return if no chunking is needed
if max(d) < 2
    return
end

% store original configuration
cfg.orig = cfg;

% loop over chunks
for i = 1:max(d)
    cfg.task = cfg.orig.task(d==i);
    
    % get correct fields for this task chunk
    for j = 1:length(cfg.task)
        task = cfg.task{j};
        idx = find(d==i,1,'first') + j - 1;
        if c(idx)
            cfg.(task) = cfg.orig.(task){c(idx)};
        end
    end
    
    % execute task chunk
    cfg.save = 'no';
    [cfg,data] = km_dotask(cfg,data);
    
    % update original fields
    for j = 1:length(cfg.task)
        task = cfg.task{j};
        idx = find(d==i,1,'first') + j - 1;
        if c(idx) > 0
            cfg.orig.(task){c(idx)} = cfg.(task);
        end
    end
    
end

% restore fields
for j = 1:length(cfg.orig.task)
    task = cfg.orig.task{j};
    if c(j) > 0
        cfg.(task) = cfg.orig.(task);
    end
end
cfg.save = cfg.orig.save;
cfg.task = cfg.orig.task;
cfg = rmfield(cfg,'orig');

% save the data if so requested
km_save(cfg,data);

% prompt the caller to return
rtrnflg = true;