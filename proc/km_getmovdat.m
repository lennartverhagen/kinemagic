function varargout= km_getmovdat(cfg,data,mov,tseries,reshapeflg)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------
% reshape the data output by default
if nargin < 5
    reshapeflg = true;
end

% FIXME: legacy reasons for app > apt
if ~isempty(regexp(tseries,'grip','once')) && isfield(data,'grip')
    tseries = regexprep(tseries,'app','apt');
    if isfield(data.grip,'appaxis')
        data.grip.aptaxis = data.grip.appaxis;
        data.grip = rmfield(data.grip,'appaxis');
    end
    if isfield(data.grip,'app')
        data.grip.apt = data.grip.app;
        data.grip = rmfield(data.grip,'app');
    end
    if isfield(data.grip,'appvel')
        data.grip.aptvel = data.grip.appvel;
        data.grip = rmfield(data.grip,'appvel');
    end
    if isfield(data.grip,'appacc')
        data.grip.aptacc = data.grip.appacc;
        data.grip = rmfield(data.grip,'appacc');
    end
end

% get indeces
[movnr marker axs] = km_getmovidx(cfg,tseries);

datalabel	= data.label;
datatime	= data.time;
switch lower(tseries)
    case {'pos','ori'}
        dataaxis	= data.axis;
    case {'vel','acc'}
        dataaxis	= data.derivaxis;
    case {'gripapt','gripaptvel','gripaptacc'}
        data        = data.grip;
        datalabel	= data.label;
        dataaxis	= data.aptaxis;
        tseries     = tseries(5:end);
    case {'gripori','griporivel'}
        data        = data.grip;
        datalabel	= data.label;
        dataaxis	= data.oriaxis;
        tseries     = tseries(5:end);
end

% select marker (try to keep the order as specified)
if isempty(marker), marker = datalabel; end;
if length(marker) > 1 && ~all(ismember(marker,datalabel))
    warning('KM:GetMovDat:MarkerSel','Marker order can not be guaranteed');
end
im = km_labelselection(marker,datalabel);
if all(ismember(im,marker)) && all(ismember(marker,im))
    [~,im] = ismember(marker,datalabel);
    nm = length(im);
else
    im = ismember(datalabel,im);
    nm = sum(im);
end

% select axes (try to keep the order as specified)
if isempty(axs), axs = dataaxis; end;
if length(axs) > 1 && ~all(ismember(axs,dataaxis))
    warning('KM:GetMovDat:AxesSel','Axes order can not be guaranteed');
end
ia = km_labelselection(axs,dataaxis);
if all(ismember(ia,axs)) && all(ismember(axs,ia))
    [~,ia] = ismember(axs,dataaxis);
    na = length(ia);
else
    ia = ismember(dataaxis,ia);
    na = sum(ia);
end

% check if markers and axis exist
if nm < 1
    error('the requested marker ''%s'' was not present in the data.',marker{1});
end
if na < 1
    error('the requested axis ''%s'' was not present in the data.',axs{1});
end

% make sure that all trials have the same number of movements
mov = km_setequalnmov(mov);
tim = cellfun(@(x) x(movnr,:),mov,'UniformOutput',false);

% first find the time indeces, this is needed to determine how many NaNs
% should be filled in when no movement is present
it = cell(1,length(tim));
for t = 1:length(tim)
    it{t} = time2logic(tim{t},datatime{t});
end
nit = cellfun(@sum,it);
if all(nit(nit>0)==1)
    nt = 1;
else
    nt = 101;
end

dat = cell(1,length(tim));
for t = 1:length(tim)
    if ~any(it{t})
        tmp = nan(nm,na,nt);
        %tmp = nan(nm,na,101);
        %tmp = nan(nm,na,1);
    else
        tmp = data.(tseries){t}(im,ia,it{t});
    end
    dims = size(tmp);
    if length(dims)==3 && dims(1) == 1 && reshapeflg
        dat{t} = reshape(tmp,dims(2:3));
    elseif length(dims)==2 && nm > 1 && nt == 101
        % FIXME: I see no other solution then to fake a third dimension by
        % adding a bunch of NaNs.
        %dat{t} = repmat(tmp,[1 1 2]);
        tmp(:,:,2) = nan(size(tmp));
        dat{t} = tmp;
    else
        dat{t} = tmp;
    end
end


% sort output
if nargout > 1
    varargout = {dat,dataaxis(ia),datalabel(im)};
else
    varargout = {dat};
end