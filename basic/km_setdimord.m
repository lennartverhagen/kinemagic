function data = km_setdimord(data,dimord_data,dimord_goal)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% return if nothing to do
if isempty(data) || (isstruct(data)&& isequal(fieldnames(data),{'dimord'}))
    return
end

% set default desired dimension order
if nargin < 3
    dimord_goal = 'marker_axis_time';
end

% use a cell instead of a character array
if ~iscell(dimord_data)
    dimord_data = regexpi(dimord_data,'_','split');
end
if ~iscell(dimord_goal)
    dimord_goal = regexpi(dimord_goal,'_','split');
end

% check if the same dimensions are used
dimdiff = setdiff(dimord_goal,dimord_data);
if ~isempty(dimdiff)
    error('the requiered dimension ''%s'' was not present in the data',dimdiff{1})
end
if length(dimord_data) ~= length(dimord_goal)
    error('the number of dimensions in the data (%d) does not match the requiered number of dimensions (%d)',length(dimord_data),length(dimord_goal));
end
    
% return if nothing to do
if isequal(dimord_data,dimord_goal)
    return
end

% set desired dimension order
if isstruct(data)
    for t = 1:length(data.pos)
        data.pos{t} = permute(data.pos{t},finddim(dimord_data,dimord_goal));
    end
elseif iscell(data)
    for t = 1:length(data)
        data{t} = permute(data{t},finddim(dimord_data,dimord_goal));
    end
else
    data = permute(data,finddim(dimord_data,dimord_goal));
end
