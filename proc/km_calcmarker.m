function [cfg,data] = km_calcmarker(cfg,data)
%--------------------------------------------------------------------------
%
%
% example configuration
% cfg.calcmarker(1).name = 'mean_AB';
% cfg.calcmarker(1).type = 'mean';
% cfg.calcmarker(1).markers = {'A','B'};
% cfg.calcmarker(2).name = 'grip_CD';
% cfg.calcmarker(2).type = 'grip';
% cfg.calcmarker(2).markers = {'C','D'};
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


% loop over requested calculated markers
for c = 1:length(cfg.calcmarker)
    switch lower(cfg.calcmarker(c).type)
        
        
        % do calculate a new marker
        %----------------------------------------
        case {'no','none'}
            
            
        % calculate a grip marker
        %----------------------------------------
        case 'grip'
            
            error('Use km_grip to calculate grip markers.')
            
            
        % calculate a new marker based on the mean of others
        %----------------------------------------
        case 'mean'
            
            % select markers
            nmarkers = length(data.label);
            markersel = km_labelselection(cfg.calcmarker(c).marker,data.label);
            markersel = match_str(data.label,markersel);
            
            % position data
            if isfield(data,'pos')
                for r = 1:length(data.pos)
                    data.pos{r}(nmarkers+1,:,:) = squeeze(mean(data.pos{r}(markersel,:,:),1));
                end
            end
            
            % orientation data
            if isfield(data,'ori')
                for r = 1:length(data.ori)
                    data.ori{r}(nmarkers+1,:,:) = squeeze(mean(data.ori{r}(markersel,:,:),1));
                end
            end
            
            % update labels
            if isfield(cfg.calcmarker(c),'name')
                data.label{nmarkers+1} = cfg.calcmarker(c).name;
            else
                data.label{nmarkers+1} = strcat('mean_',data.label{markersel});
            end
        
            
        % not recognized
        %----------------------------------------
        otherwise
            error('calculation type %s not supported yet',cfg.calcmarker(c).type)
    end
end
