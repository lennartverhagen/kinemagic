function varargout = km_artfctdef(cfg,data)
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


%% Determine artifact according to task
%----------------------------------------
% loop over artifact types
if ~iscell(cfg.artfctdef.type),   cfg.artfctdef.type = {cfg.artfctdef.type};  end
artifact = cell(1,length(cfg.artfctdef.type));
for t = 1:length(cfg.artfctdef.type)
    
    % initialize type field
    tname = cfg.artfctdef.type{t};
    if ~isfield(cfg.artfctdef,tname),	cfg.artfctdef.(tname) = [];     end
    if strcmpi(tname,'log'),
        if ~isfield(cfg,'log'),       	error('no log provided');       end
        cfg.artfctdef.log.log = cfg.log;
    end
    
    % evaluate appropriate function
    artname = ['artfctdef_' tname];
    cfg.artfctdef.(tname) = feval(['km_' artname],cfg.artfctdef.(tname),data);
    
    % store in configuration structure
    artifact{t} = cfg.artfctdef.(tname).artifact;
    
end

% combine artifacts if needed
if length(artifact) > 1
    oldart = artifact;
    artifact = cell(1,length(data.time));
    for t = 1:length(data.time);
        art = false(length(oldart),length(data.time{t}));
        for a = 1:length(oldart)
            art(a,:) = time2logic(oldart{a}{t},data.time{t});
        end
        artifact{t} = logic2time(any(art,1),data.time{t});
    end
else
    artifact = artifact{1};
end

% return artifact or store combined artifact in data and configuration structure
if nargout == 1
    varargout = {artifact};
else
    cfg.artfctdef.artifact = artifact;
    cfg.artifact	= artifact;
    data.artifact	= artifact;
    varargout       = {cfg,data};
end