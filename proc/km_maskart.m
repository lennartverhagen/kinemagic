function [cfg,data] = km_maskart(cfg,data)
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


% get artifact
if isfield(cfg,'artifact') && ~isstruct(cfg.artifact) && ~isempty(cfg.artifact)
    artifact = cfg.artifact;
elseif isfield(cfg,'artfctdef') && isfield(cfg.artfctdef,'artifact')
    artifact = cfg.artfctdef.artifact;
elseif isfield(data,'artifact')
    artifact = data.artifact;
elseif isfield(cfg,'artfctdef')
    artifact = km_artfctdef(cfg,data);
else
    error('no artifact provided');
end

% check artifact
if ~iscell(artifact) || length(artifact) == length(data.pos)
    artifact = {artifact};
end

% store correct artifact back in configuration structure for output
nart = length(artifact);
if nart==1
    cfg.artifact = artifact{1};
else
    cfg.artifact = artifact;
end

% make sure the artifact masking configuration matches the number of artifacts
if ~iscell(cfg.maskart),    cfg.maskart = {cfg.maskart};    end
if length(cfg.maskart) == 1 && nart > 1
    cfg.maskart = repmat(cfg.maskart,1,nart);
end

% loop over artifacts to mask them
for a = 1:nart
    
    % set configuration
    tcfg = [];
    tcfg.maskart	= cfg.maskart{a};
    tcfg.artifact	= artifact{a};
    
    % mask artifact
    [~, data] = maskartifact(tcfg,data);
    
end

% update dataset
if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
cfg.dataset = [cfg.dataset '_MASK' cfg.maskart{1}];
%--------------------------------------------------------------------------


%% function maskartifact
%----------------------------------------
function [cfg,data] = maskartifact(cfg,data)

% loop over runs
param = intersect(fieldnames(data),{'pos','ori'});
for r = 1:length(data.time)
    
    % go from the time to the sample domain
    art = time2logic(cfg.artifact{r},data.time{r});
    
    % continue if no artifact
    if ~any(art(:)),    continue;   end
    
    % mask the artifact using the specified mode
    switch lower(cfg.maskart)
        
        case 'no'
            
        case 'cut'
            
            % cut artefact out of data and time
            for p = 1:length(param)
                data.(param{p}){r} = data.(param{p}){r}(:,~art);
            end
            data.time{r} = data.time{r}(~art);
        
        case 'nan'
            
            % replace data with NaNs
            for p = 1:length(param)
                data.(param{p}){r}(:,:,art) = NaN;
            end
            
        case {'nearest','linear','spline','pchip','cubic'}
            
            for p = 1:length(param)
                
                % padded artifact
                npadsamp = 100;
                idx = conv(single(art),ones(1,npadsamp),'same')>0;
                
                % original time
                xi = data.time{r}(idx);
                xi = xi(:);
                
                % time without artifact
                x = xi(~art(idx));
                
                % data without artifact
                y = data.(param{p}){r}(:,:,idx);
                y = y(:,:,~art(idx));
                dims = size(y);
                y = permute(y,[3 2 1]);
                y = reshape(y,dims([3 2 1]));
                
                % interpolate
                yi = interp1(x,y,xi,lower(cfg.maskart),'extrap');
                
                % reshape data to fit original
                dims(3) = size(yi,1);
                yi = permute(yi,[3 2 1]);
                yi = reshape(yi,dims);
                data.(param{p}){r}(:,:,idx) = yi;
                
                % old code (interpolation on all time points)
                % % original time
                % xi = data.time{r};
                % 
                % % time without artifact
                % x = data.time{r}(~art);
                % 
                % % data without artifact
                % y = data.(param{p}){r}(:,:,~art);
                % y = permute(y,[3 2 1]);
                % 
                % interpolate
                % yi = interp1(x,y,xi,lower(cfg.maskart),'extrap');
                % yi = permute(yi,[3 2 1]);
                % data.(param{p}){r} = yi;
            end
            
        otherwise
            warning('KM:NotSupported','%s: this option (%s) is not supported.',mfilename,cfg.maskart);
            
    end
    
end
