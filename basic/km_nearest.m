function [idx] = km_nearest(array, val)

% NEAREST return the index of an array nearest to a scalar
% 
% [indx] = nearest(array, val)

% Copyright (C) 2002, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% % try to use the standard version supplied by FieldTrip, otherwise use the
% % script below.
% funname   = mfilename;
% funname   = funname(4:end);
% if exist(funname,'file') == 2
%     funhandle = str2func(funname);
%     [i] = funhandle(array, val);
%     return
% end

if isempty(array)
    idx = [];
    return
end

mbreal(array);
mbreal(val);

mbvector(array);
%mbscalar(val);
mbvector(val);

% ensure that it is a column vector
array = array(:);
val = val(:);

if any(isnan(val))
  error('incorrect value')
end

% loop over val elements
idx = nan(size(val));
for i = 1:length(val)
    if val>max(array)
        % return the last occurence of the nearest number
        [~, idx(i)] = max(flipud(array));
        idx(i) = length(array) + 1 - idx(i);
    else
        % return the first occurence of the nearest number
        [~, idx(i)] = min(abs(array(:) - val(i)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbreal(a)
if ~isreal(a)
  error('Argument to mbreal must be real');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbscalar(a)
if ~all(size(a)==1)
  error('Argument to mbscalar must be scalar');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mbvector(a)
if ndims(a) > 2 || (size(a, 1) > 1 && size(a, 2) > 1)
  error('Argument to mbvector must be a vector');
end

