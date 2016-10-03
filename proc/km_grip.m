function [cfg,data] = km_grip(cfg,data)
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
if length(data.axis) < 2
    cfg.grip.ori = false;
    if ~isfield(cfg.grip,'apt'),    cfg.grip.apt = false;   end
end
if isfalse(cfg.(task)) || (isfalse(cfg.grip.apt) && isfalse(cfg.grip.ori))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end


% remove grip field if present
if isfield(data,'grip')
    warning('KM:RemoveField','The field ''grip'' is present in the data. It will be replaced.');
    data = rmfield(data,'grip');
end

% initialize grip field
data.grip.label     = cfg.grip.label;
data.grip.dimord    = 'marker_axis_time';

% select aperture axes
if ~isfalse(cfg.grip.apt)
    if ischar(cfg.grip.apt.axis{1})
        aptaxis = {'xyz','xy','xz','yz','x','y','z'};
        aptaxis = km_labelselection(cfg.grip.apt.axis,aptaxis);
    else
        aptaxis = cfg.grip.apt.axis;
    end
    data.grip.aptaxis	= aptaxis;
else
    aptaxis = [];
end

% select orientation axes
if ~isfalse(cfg.grip.ori)
    if ischar(cfg.grip.ori.axis{1})
        oriaxis = {'xy','xz','yz','x','y','z'};
        oriaxis = km_labelselection(cfg.grip.ori.axis,oriaxis);
    else
        oriaxis = cfg.grip.ori.axis;
    end
    if length(oriaxis) > length(cfg.grip.ori.dir)
        cfg.grip.ori.dir = repmat({cfg.grip.ori.dir{1}},1,length(oriaxis));
    end
    if length(oriaxis) > length(cfg.grip.ori.offset)
        cfg.grip.ori.offset = repmat(cfg.grip.ori.offset(1),1,length(oriaxis));
    end
    data.grip.oriaxis	= oriaxis;
else
    oriaxis = [];
end

% loop over runs
for r = 1:length(data.pos)
    
    % intitialize data
    if ~isfalse(cfg.grip.apt)
        data.grip.apt{r}	= [];
    end
    if ~isfalse(cfg.grip.ori)
        data.grip.ori{r}	= [];
    end
    
    % loop over grip labels
    for g = 1:length(cfg.grip.label)
        
        % select markers
        markersel = km_labelselection(cfg.grip.marker{g},data.label);
        markersel = match_str(data.label,markersel);
    
        % check number of selected
        if length(markersel) ~= 2
            error('To calculate grip markers, two markers need to be selected, not %d. Maybe you forgot to put the marker selection in a cell array?',length(markersel))
        end
        
        % loop over aperture axes
        for aa = 1:length(aptaxis)
            
            % select aperture axes
            axissel = km_axisselection(aptaxis{aa},data.axis);
            axissel = ismember(data.axis,axissel);
            
            % get aperture
            apt = calc_apt(data.pos{r}(markersel,axissel,:),cfg.grip.apt.abs);
            
            data.grip.apt{r}(g,aa,:) = apt;
            
        end
        
        % loop over orientation axes
        for oa = 1:length(oriaxis)
            
            % get orientation direction
            if isnumeric(cfg.grip.ori.dir{oa})
                oridir = cfg.grip.ori.dir{oa};
            elseif strcmpi(cfg.grip.ori.dir{oa},'ccw')
                oridir = 1;
            elseif strcmpi(cfg.grip.ori.dir{oa},'cw')
                oridir = -1;
            end
            
            % select orientation plane (e.g. 'yz' orientation IN the yz-plane,
            % 'z' orientation relative to the z-axis/xy-plane
            axissel = km_axisselection(oriaxis{oa},data.axis);
            [~,axisidx] = ismember(axissel,data.axis);
            axissel = ismember(data.axis,axissel);
            
            % get orientation
            ori = calc_ori(data.pos{r}(markersel,:,:),axisidx,axissel);

            % apply orientation direction and offset
            ori = oridir*ori + cfg.grip.ori.offset(oa);
            data.grip.ori{r}(g,oa,:) = ori;
            
        end
        
    end
    
end


% function calc_apt
%----------------------------------------
function apt = calc_apt(dat,absmethod)

if nargin < 2,  absmethod = 'no';  end

apt = squeeze(diff(dat,1));
if ndims(apt) > 1
    % calculate grip aperture in Euclidian distance (using Pythagoras)
    apt = squeeze(sqrt(sum(apt.^2,1)));
elseif istrue(absmethod,'free')
	apt = abs(apt);
end
    

% function calc_ori
%----------------------------------------
function ori = calc_ori(data,axisidx,axissel)

% get adjacent (b) and opposide (a) sides of a right triangle
% to calculate grip orientation
if length(axisidx) == 2
    % a is the second selected axis of the rotation plane
    % b is the first selected axis of the rotation plane
    dat = data(:,axisidx,:);
    a = squeeze(diff(dat(:,2,:)));
    b = squeeze(diff(dat(:,1,:)));
elseif length(axisidx) == 1
    % a is the reference axis
    % b is the distance in the reference plane
    dat = data(:,axisidx,:);
    a = squeeze(diff(dat));
    dat = data(:,~axissel,:);
    b = squeeze(sqrt(sum(diff(dat).^2,2)));
else
    error('An orientation plane can only be described by one or two axes, not %d',length(axisidx))
end

% calculate grip orientation using invert tangens (SOH-CAH-TOA)
% and convert radians to degrees
ori = atan(a./b);
ori = ori/pi*180;

% where the adjacent side is negative, please add 180 degrees to get a full
% 360 degree range
ori(b<0) = ori(b<0)+180;