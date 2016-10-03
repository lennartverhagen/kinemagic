function varargout = km_movselect(cfg,data)
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

% get movement
if isfield(cfg,'movement') && ~isstruct(cfg.movement) && ~isempty(cfg.movement)
    movement = cfg.movement;
elseif isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
    movement = cfg.movdef.movement;
elseif isfield(data,'movement')
    movement = data.movement;
else
    error('no movement provided');
end
oldmovement = movement;

% make sure that all trials have the same number of movements
movement = km_setequalnmov(movement);

% use km_movparam to obtain the criterium for each movement
if strcmpi(cfg.movselect.param{1},'index')
    val = repmat(1:size(movement{1},1),[length(movement) 1]);
else
    val = nan(length(movement),size(movement{1},1));
    for m = 1:size(movement{1},1)
        tcfg = cfg;
        tcfg.movparam = cfg.movselect;
        tcfg.movparam.movnr = m;
        tcfg = km_movparam(tcfg,data);
        val(:,m) = tcfg.trl(:,ismember(tcfg.vars,tcfg.movparam.pname));
    end
end

% sort the data
val(isnan(val)) = -Inf;
[~,idx] = sort(val,2,cfg.movselect.sort);

% select the movements
if ~isfalse(cfg.movselect.movidx)
    idx = idx(:,cfg.movselect.movidx);
else
    idx = idx(:,1:max(cfg.movselect.nmov));
end
movement = cellfun(@(x,y) x(y,:),movement,mat2cell(idx,ones(size(idx,1),1))','UniformOutput',false);

% check for correct number of movements
nmov = cellfun(@(x) sum(any(~isnan(x),2)),movement);
idx_excl = nmov < min(cfg.movselect.nmov) | nmov > max(cfg.movselect.nmov);
movement(idx_excl) = {nan(size(movement{1}))};

% store
cfg.movselect.val = val;
cfg.movselect.idx = idx;
cfg.movselect.allmov = oldmovement;

% return movement definition or store in data and configuration structure
if nargout == 1
    varargout = {movement};
else
    if length(cfg.movdef) == 1
        cfg.movdef.movement = movement;
    end
    cfg.movement	= movement;
    data.movement	= movement;
    varargout       = {cfg,data};
end