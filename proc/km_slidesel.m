function select = km_slidesel(data,select,crit,consecflg)
%--------------------------------------------------------------------------
% This algorithm does some wicked stuff!
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% check input
if nargin < 3,  crit = 0;               end
if ~isnumeric(crit) && istrue(crit,'free'), crit = 0; end
if nargin < 4,  consecflg = 'onset';	end
if all(isnan(data(:))), return;         end
n = numel(data);
origdims = size(data);
if ~isequal(size(select),origdims)
    error('KM:SlideSel:InputDims','The size of the data does not match the size of the selection');
end
if n==max(origdims)
    data = reshape(data,[1 1 n]);
    select = reshape(select,[1 1 n]);
end
if ndims(data) ~= 3
    error('KM:SlideSel:InputDims','The data must be of size Nx1, 1xN or IxJxN');
end    

% loop over every dimension
for m = 1:size(data,1)
    for a = 1:size(data,2)
        
        % select data
        dat = squeeze(data(m,a,:));
        sel = squeeze(select(m,a,:));
                
        % determine the slope of the data
        slope = sign(gradient(dat));
        
        % apply criterium
        critsel = false(size(dat,1),2);
        critsel(dat>crit,1) = true;
        critsel(dat<-crit,2) = true;
        
        % combine selections on the positive slope
        pos = slope > 0;
        sel = selcomb(sel,pos,critsel,consecflg);
                
        % combine selections on the negative slope
        critsel = critsel(:,[2 1]);
        neg = slope < 0;
        sel = selcomb(sel,neg,critsel,consecflg);
        
        % split selections that cross signs if requested
        sel = signsplitsel(sel,dat,consecflg);
        
        % store selection
        select(m,a,:) = sel;
        
    end
end

% reshape select back
if n==max(origdims)
    select = reshape(select,origdims);
end


% function selcomb
%----------------------------------------
function sel = selcomb(sel,slope,crit,flg)

% transform selections to on- and off-sets
nsamp = length(sel);
sel = logic2idx(sel);
slopesel = logic2idx(slope);
oncrit = logic2idx(crit(:,1));
offcrit = logic2idx(crit(:,2));

% loop over the selections
for i = 1:size(sel,1)
    % movement onset > slide to the left
    if slope(sel(i,1))
        % slide to the left
        tmpslope = min(slopesel(sel(i,1) <= slopesel(:,2),1));
        % apply the onset criterium
        tmpcrit = min([sel(i,1) min(oncrit(sel(i,1) < oncrit(:,2),1))]);
        % combine
        sel(i,1) = max([tmpslope tmpcrit]);
    end
    % movement offset > slide to the right
    if slope(sel(i,2))
        % slide to the right
        tmpslope = max(slopesel(sel(i,2) >= slopesel(:,1),2));
        % apply the offset criterium
        tmpcrit = max([sel(i,2) max(offcrit(sel(i,2) > offcrit(:,1),2))]);
        % combine
        sel(i,2) = min([tmpslope tmpcrit]);
    end
end

% force separation of consecutive movements if requested
if ~isfalse(flg)
    % selections that share the exact same point (as on- and off-sets)
    idx = find(sel(1:end-1,2) == sel(2:end,1));
    sel = shiftselidx(sel,idx,2,flg);
    % selections that share consecutive points (as on- and off-sets)
    idx = find(sel(1:end-1,2)+1 == sel(2:end,1));
    sel = shiftselidx(sel,idx,1,flg);
end

% transform selection to indeces
sel = idx2logic(sel,nsamp)';

% function shiftselidx
%----------------------------------------
function sel = shiftselidx(sel,idx,n,flg)
if isempty(idx),    return; end
switch lower(flg)
    case 'onset'
        sel(idx+1,1) = sel(idx,2) + n;
    case 'offset'
        sel(idx,2) = sel(idx+1,1) - n;
    case 'both'
        sel(idx+1,1) = sel(idx,2) + n;
        sel(idx,2) = sel(idx+1,1) - n - 1;
    otherwise
end

% function signsplitsel
%----------------------------------------
function sel = signsplitsel(sel,dat,flg)
% return quickly
if isempty(regexp(flg,'sign','once'))
    return;
end
% transform selections to indeces
nsamp = length(sel);
sel = logic2idx(sel);
% find negative data sections
negdat = dat < 0;

% loop over the selections
for i = 1:size(sel,1)
    % check data sign within selections
    tmpsign = negdat(sel(i,1):sel(i,2));
    % if sign changes within selection, split sel
    if any(diff(tmpsign))
        % split sel
        tmp = find(diff([-2; tmpsign])~=0);
        tmpsel = [tmp [tmp(2:end)-1; length(tmpsign)]] - 1 + sel(i,1);
        % move onset, offset, or both
        if ~isempty(regexp(flg,'offset','once'))
            % move only offset
            tmpsel(1:end-1,1) = tmpsel(1:end-1,1) - 1;
        elseif ~isempty(regexp(flg,'both','once'))
            % move onset and offset
            tmpsel(2:end,1) = tmpsel(2:end,1) + 1;
            tmpsel(1:end-1,1) = tmpsel(1:end-1,1) - 1;
        else %~isempty(regexp(flg,'onset','once'))
            % move only onset
            tmpsel(2:end,1) = tmpsel(2:end,1) + 1;
        end
        % store tmpsel in sel
        sel(i,:) = tmpsel(1,:);
        sel = [sel; tmpsel(2:end,:)];
    end
end

% transform selection to indeces
sel = idx2logic(sel,nsamp)';