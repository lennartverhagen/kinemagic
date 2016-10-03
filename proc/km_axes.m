function [cfg,data] = km_axes(cfg,data)
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


% re-order the axes
[cfg,data] = orderaxes(cfg,data);

% mirror axes
[cfg,data] = mirroraxes(cfg,data);

% rotate axes
[cfg,data] = rotateaxes(cfg,data);


%% function orderaxes
%----------------------------------------
function [cfg,data] = orderaxes(cfg,data)
% re-order axes

% return if nothing to do
if isfalse(cfg.axes.order)
    return
end

% check configuration
if ischar(cfg.axes.order)
    cfg.axes.order = num2cell(cfg.axes.order);
elseif ~iscell(cfg.axes.order)
    cfg.axes.order = data.axis(cfg.axes.order);
end

order = finddim(data.axis,cfg.axes.order);
if ~isequal(order,1:length(data.axis))
    for r = 1:length(data.pos)
        data.pos{r} = data.pos{r}(:,order,:);
    end
end


%% function mirroraxes
%----------------------------------------
function [cfg,data] = mirroraxes(cfg,data)
% Change signs

% return if nothing to do
if isfalse(cfg.axes.mirror)
    return
end

% check configuration
if ischar(cfg.axes.mirror)
    cfg.axes.mirror = num2cell(cfg.axes.mirror);
elseif ~iscell(cfg.axes.order)
    cfg.axes.mirror = data.axis(cfg.axes.mirror);
end

mirror = ismember(data.axis,cfg.axes.mirror);
if any(mirror)
    for r = 1:length(data.pos)
        data.pos{r}(:,mirror,:) = -1*data.pos{r}(:,mirror,:);
    end
end


%% function rotateaxes
%----------------------------------------
function [cfg,data] = rotateaxes(cfg,data)
% Rotate reference frame
% Rotates datapoints around axis CLOCKWISE when looking to origin
% So the reference frame is rotated counter clockwise!
% in degrees around [azimuth, elevation, tilt]
% 1) tilt      (yz), angle to y-axis, round X
% 2) elevation (xz), angle to x-axis, round Y 
% 3) azimuth   (xy), angle to z-axis, round Z
%
%http://mathworld.wolfram.com/RotationMatrix.html

% return if nothing to do
if isfalse(cfg.axes.rotate) || length(cfg.axes.rotate) ~= length(data.axis)
    cfg.axes.rotate = zeros(1,length(data.axis));
    return
end

if any(cfg.axes.rotate)
    
    rx = cfg.axes.rotate(1)/180*pi;
    ry = cfg.axes.rotate(2)/180*pi;
    rz = cfg.axes.rotate(3)/180*pi;

    r1 = [  1       0       0
            0       cos(rx) sin(rx)
            0      -sin(rx) cos(rx) ];
        
    r2 = [  cos(ry) 0      -sin(ry)
            0       1       0
            sin(ry) 0       cos(ry) ];
        
    r3 = [  cos(rz) sin(rz) 0
           -sin(rz) cos(rz) 0
            0       0       1       ];
        
    rot_matrix = r3*r2*r1;       

    for r = 1:length(data.pos)
        for m = 1:data.marker
            data.pos{r}(m,:,:) = data.pos{r}(m,:,:) * rot_matrix;
        end
    end
    
end