function cfg = km_artfctdef_tms(cfg,data)
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
if ~isfield(cfg,'critvel'),         cfg.critvel     = 1.2;                  end     % velocity in m/s
if ~isfield(cfg,'critacc'),         cfg.critacc     = 0.08;                 end     % velocity in m/s²
if ~isfield(cfg,'prepad'),          cfg.prepad      = -0.005;               end     % prepad artifact in seconds
if ~isfield(cfg,'postpad'),         cfg.postpad     = 0.020;                end     % postpad artifact in seconds

% detect TMS artifacts
cfg.artifact = cell(1,length(data.pos));
for r = 1:length(data.pos)
    
    % calculate velocity and acceleration
    pos = data.pos{r};
    vel = nan(size(pos));
    vel(:,:,2:end) = abs(diff(pos,1,3))*data.fsample(r);
    acc = nan(size(pos));
    acc(:,:,2:end) = abs(diff(vel,1,3))*data.fsample(r);
    
    % apply criteria
    art = vel > cfg.critvel & acc > cfg.critacc;
    
    % combine markers and axes
    art = squeeze(any(any(art,1),2));
    
    % go from the sample to the time domain
    art = logic2time(art,data.time{r});
    
    % add padding
    art = [art(:,1)+cfg.prepad art(:,2)+cfg.postpad];
    
    % go to sample domain and back to combine overlapping artifacts
    art = time2logic(art,data.time{r});
    art = logic2time(art,data.time{r});
    
    % store the artifact
    cfg.artifact{r} = art;
    
end