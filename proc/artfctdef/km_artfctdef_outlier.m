function cfg = km_artfctdef_outlier(cfg,data)
%--------------------------------------------------------------------------
%
% See also KM_ARTFCTDEF
%
% FIXME: include option to use the inter-quartile range instead of the
% z-value to define outliers
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% set configuration
if isfield(cfg,'type'),             cfg             = rmfield(cfg,'type');  end
if ~isfield(cfg,'param'),           cfg.param       = {'pos','vel'};        end     % parameters to take into account
if ~iscell(cfg.param),              cfg.param       = {cfg.param};          end
if ~isfield(cfg,'critpos'),         cfg.critpos     = 1;                    end     % critical position in meters
if ~isfield(cfg,'critvel'),         cfg.critvel     = 10;                   end     % critical velocity in meters/second
if ~isfield(cfg,'critacc'),         cfg.critacc     = 2000;                 end     % critical acceleration in meters/second^2
if ~isfield(cfg,'critzval'),        cfg.critzval    = 20;                   end     % critical z-value
if ~isfield(cfg,'prepad'),          cfg.prepad      = -0.002;               end     % prepad artifact in seconds
if ~isfield(cfg,'postpad'),         cfg.postpad     = 0.002;                end     % postpad artifact in seconds

% detect z-value artifacts
cfg.artifact = cell(1,length(data.pos));
for r = 1:length(data.pos)
    
    % loop over parameters
    art = false(size(data.pos{r}));
    for p = 1:length(cfg.param)
        pname = cfg.param{p};
        
        % get data
        switch lower(pname)
            case 'pos'
                dat = data.pos{r};
            case 'vel'
                pos = data.pos{r};
                dat = nan(size(pos));
                dat(:,:,2:end) = abs(diff(pos,1,3))*data.fsample(r);
            case 'acc'
                pos = data.pos{r};
                vel = nan(size(pos));
                vel(:,:,2:end) = abs(diff(pos,1,3))*data.fsample(r);
                dat = nan(size(pos));
                dat(:,:,2:end) = abs(diff(vel,1,3))*data.fsample(r);
        end
    
        % apply extreme values criteria
        fldname = ['crit' pname];
        if ~isfalse(cfg.(fldname)),   art = art | abs(dat) > cfg.(fldname); end
    
        % plot for testing
        %if any(art(:))
        %    figure(1); hold off; plot(squeeze(pos)*20,'r'); hold on; plot(squeeze(vel)/5,'g'); plot(squeeze(acc)/1000,'b'); plot(squeeze(art),'k','linewidth',2); hold off;
        %end
    
        % z-transform the data
        mudat = repmat(nanmean(dat,3),[1 1 size(dat,3)]);
        stddat = repmat(nanstd(dat,0,3),[1 1 size(dat,3)]);
        dat = abs((dat-mudat)./stddat);
    
        % apply z-value criteria
        art = art | dat > cfg.critzval;
    end
    
    % combine markers and axes
    art = squeeze(any(any(art,1),2));
    
    % go from the sample to the time domain
    art = logic2time(art,data.time{r});
    
    % add padding
    art = [art(:,1)+cfg.prepad art(:,2)+cfg.postpad];
    
    % go to sample domain and back to combine overlapping artifacts
    art = time2logic(art,data.time{r});
    art = logic2time(art,data.time{r});
    
    % plot for testing
    %if ~isempty(art)
    %    figure(2); hold off; plot(squeeze(pos),'r'); hold on; plot(squeeze(vel),'g'); plot(squeeze(acc),'b'); plot(time2logic(art,data.time{r})*10,'k','linewidth',2); hold off;
    %    disp(r);
    %    A = 1;
    %end
    
    % store the artifact
    cfg.artifact{r} = art;
    
end