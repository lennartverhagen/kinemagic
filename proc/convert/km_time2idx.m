function varargout = km_time2idx(varargin)
%--------------------------------------------------------------------------
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

funname   = mfilename;
funhandle = str2func(funname(4:end));
[varargout{1:nargout}] = funhandle(varargin{:});


%% OLD CODE
% function idx = km_time2idx(time_points,time_array)
% %--------------------------------------------------------------------------
% %
% % This file is part of the KineMagic toolbox
% % Copyright (C) 2010, Lennart Verhagen
% % L.Verhagen@donders.ru.nl
% % version 2010-01-01
% %--------------------------------------------------------------------------
% 
% % transpose time_points if arranged in the second dimension
% if size(time_points,2) > 2
%     time_points = time_points';
% end
% 
% idx = nan(size(time_points));
% if isempty(time_points)
%     return
% elseif size(time_points,2) == 1
%     for t = 1:size(time_points,1)
%         idx(t) = km_nearest(time_array,time_points(t));        
%     end
% elseif size(time_points,2) == 2
%     for t = 1:size(time_points,1)
%         idx_on = km_nearest(time_array,time_points(t,1));
%         idx_off = km_nearest(time_array,time_points(t,2));
%         idx(t,:) = [idx_on idx_off];
%     end
% else
%     error('Unsupported time_point matrix size. Size should be nx1 or nx2.')
% end