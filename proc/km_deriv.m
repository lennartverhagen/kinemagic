function [cfg,data] = km_deriv(cfg,data)
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

% select position & orientation axes
derivaxis = [];
if isfield(data,'pos') || isfield(data,'ori')
    if ischar(cfg.deriv.axis{1})
        derivaxis = {'xyz','xy','xz','yz','x','y','z'};
        derivaxis = km_labelselection(cfg.deriv.axis,derivaxis);
    else
        derivaxis = cfg.deriv.axis;
    end
    data.derivaxis	= derivaxis;
end

% select grip axes
grpaptaxis = [];
grporiaxis = [];
if isfield(data,'grip')
    % grip aperture
    if isfield(data.grip,'apt')
        grpaptaxis = km_labelselection('all',data.grip.aptaxis);
    end
    % grip orientation
    if isfield(data.grip,'apt')
        grporiaxis = km_labelselection('all',data.grip.oriaxis);
    end
end


% loop over runs
for r = 1:length(data.time)
    
    % allocate position and orientation derivatives matrixes
    nmarkers = length(data.label);
    naxes = length(derivaxis);
    nsamp = size(data.pos{r},3);
    if isfield(data,'pos')
        if istrue(cfg.deriv.vel),   data.vel{r} = nan(nmarkers,naxes,nsamp);    end
        if istrue(cfg.deriv.acc),   data.acc{r} = nan(nmarkers,naxes,nsamp);	end
    end
    if isfield(data,'ori')
        if istrue(cfg.deriv.vel),   data.orivel{r} = nan(nmarkers,naxes,nsamp);	end
        if istrue(cfg.deriv.acc),   data.oriacc{r} = nan(nmarkers,naxes,nsamp);	end
    end
    
    % loop over position and orientation axes
    for da = 1:length(derivaxis)
        axissel = km_axisselection(derivaxis{da},data.axis);
        axissel = ismember(data.axis,axissel);
        if isfield(data,'pos')
            dat = data.pos{r}(:,axissel,:);
            vel = deriv(cfg.deriv,dat,data.fsample(r));
            if istrue(cfg.deriv.vel),   data.vel{r}(:,da,:) = vel;  end
            if istrue(cfg.deriv.acc)
                acc = deriv(cfg.deriv,vel,data.fsample(r));
                data.acc{r}(:,da,:) = acc;
            end
        end
        if isfield(data,'ori')
            dat = data.ori{r}(:,axissel,:);
            vel = deriv(cfg.deriv,dat,data.fsample(r));
            if istrue(cfg.deriv.vel),   data.orivel{r}(:,da,:) = vel;	end
            if istrue(cfg.deriv.acc)
                acc = deriv(cfg.deriv,vel,data.fsample(r));
                data.oriacc{r}(:,da,:) = acc;
            end
        end
    end
    
    % loop over grip aperture axes
    if ~isempty(grpaptaxis)
        if istrue(cfg.deriv.vel),   data.grip.aptvel{r} = nan(size(data.grip.apt{r}));	end
        if istrue(cfg.deriv.acc),   data.grip.aptacc{r} = nan(size(data.grip.apt{r}));	end
    end
    for aa = 1:length(grpaptaxis)
        axissel = ismember(data.grip.aptaxis,grpaptaxis{aa});
        dat = data.grip.apt{r}(:,axissel,:);
        vel = deriv(cfg.deriv,dat,data.fsample(r));
        if istrue(cfg.deriv.vel),   data.grip.aptvel{r}(:,aa,:) = vel;  end
        if istrue(cfg.deriv.acc)
            acc = deriv(cfg.deriv,vel,data.fsample(r));
            data.grip.aptacc{r}(:,aa,:) = acc;
        end
    end
    
    % loop over grip orientation axes
    if ~isempty(grporiaxis)
        if istrue(cfg.deriv.vel),   data.grip.orivel{r} = nan(size(data.grip.ori{r}));	end
        if istrue(cfg.deriv.acc),   data.grip.oriacc{r} = nan(size(data.grip.ori{r}));	end
    end
    for oa = 1:length(grporiaxis)
        axissel = ismember(data.grip.oriaxis,grporiaxis{oa});
        dat = data.grip.ori{r}(:,axissel,:);
        vel = deriv(cfg.deriv,dat,data.fsample(r));
        if istrue(cfg.deriv.vel),   data.grip.orivel{r}(:,oa,:) = vel;  end
        if istrue(cfg.deriv.acc)
            acc = deriv(cfg.deriv,vel,data.fsample(r));
            data.grip.oriacc{r}(:,oa,:) = acc;
        end
    end 

end


% function deriv
%----------------------------------------
function ddat = deriv(cfg,dat,fsample)

% return if no data
if isempty(dat)
    ddat = zeros(size(dat));
    warning('KM:Return','%s: Nothing to do...',task);
    return
end

switch lower(cfg.method)
    case 'gradient'
        [~, ~, gra] = gradient(dat, 1/fsample);
        if ndims(gra) == 3
            if size(gra,2) > 1
                % Euclidian (tangential) derivative
                ddat = sqrt( sum(gra.^2, 2) );
            else
                ddat = gra;
            end
        elseif ndims(gra) == 2
            ddat(:,1,:) = permute(gra,[2 1]);
        elseif ndims(gra) < 2
            ddat(1,1,:) = gra;
        end
        if istrue(cfg.abs,'free')
            ddat = abs(ddat);
        end
    case 'diff'
        % dat must have 3 dimensions
        ddat = nan(size(dat,1),1,size(dat,3));
        ddat(:,:,2:end) = diff(dat,1,3) * fsample;
        if size(dat,2) > 1
            % Euclidian (tangential) derivative
            ddat = sqrt( sum(( ddat .^2),2) );
        end
        if istrue(cfg.abs,'free')
            ddat = abs(ddat);
        end
    otherwise
        warning('KM:NotSupported','%s: This option (%s) is not supported',mfilename,method);
end
