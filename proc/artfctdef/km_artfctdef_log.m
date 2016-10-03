function cfg = km_artfctdef_log(cfg,data)
%--------------------------------------------------------------------------
%
% See also KM_ARTFCTDEF
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% set configuration
if isfield(cfg,'type'),             cfg             = rmfield(cfg,'type');  end
if ~isfield(cfg,'log'),             error('no log provided');               end

if ~isfield(cfg,'beg'),             cfg.beg         = 't_beg';              end
if ~isfield(cfg,'end'),             cfg.end         = cfg.beg;              end
if ~isfield(cfg,'prepad'),          cfg.prepad      = 0;                    end     % prepad in seconds relative to begin of artifact
if ~isfield(cfg,'postpad'),         cfg.postpad     = 0;                    end     % postpad in seconds relative to end of artifact

% check configuration
if isfield(cfg,'event')
    cfg.art = cfg.event;
    cfg = rmfield(cfg,'event');
end
if isfield(cfg,'stim')
    cfg.art = cfg.stim;
    cfg = rmfield(cfg,'stim');
end
if isfield(cfg,'art')
    cfg.beg = cfg.art;
    cfg.end = cfg.art;
end
if isfield(cfg,'offset')
    cfg.prepad = cfg.prepad + cfg.offset;
    cfg.postpad = cfg.postpad + cfg.offset;
end
    
% variable names
vars = cfg.log.vars;

% loop over chunks
cfg.artifact = cell(1,length(data.pos));
for r = 1:length(data.pos)
    
    % log data
    logdata     = cfg.log.data{cfg.log.runnr==data.run(r)};
    tim         = data.time{r};
    
    % define artifact begin and end
    idx_beg     = strcmpi(vars,cfg.beg);
    idx_end     = strcmpi(vars,cfg.end);
    timbeg      = logdata(:,idx_beg)+cfg.prepad;
    timend      = logdata(:,idx_end)+cfg.postpad;
    if isempty(timbeg), timbeg = timend;    end
    if isempty(timend), timend = timbeg;    end
    
    % keep only the time points within the data chunk
    if ~isempty(timbeg) && ~isempty(timend)
        timmin = min(tim);
        timmax = max(tim);
        idxbeg = timbeg>=timmin & timbeg<=timmax;
        idxend = timend>=timmin & timend<=timmax;
        timbeg = timbeg(all([idxbeg idxend],2));
        timend = timend(all([idxbeg idxend],2));
    end
    
    % store the artifact
    cfg.artifact{r} = [timbeg-cfg.prepad timend+cfg.postpad];
    
end
