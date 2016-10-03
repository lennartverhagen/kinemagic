function varargout = km_idx2time(varargin)
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


% %% OLD CODE
% function time_points = km_idx2time(idx,time_array)
% %--------------------------------------------------------------------------
% %
% % This file is part of the KineMagic toolbox
% % Copyright (C) 2010, Lennart Verhagen
% % L.Verhagen@donders.ru.nl
% % version 2010-01-01
% %--------------------------------------------------------------------------
% 
% % transpose indeces if arranged in the second dimension
% if size(idx,2) > 2
%     idx = idx';
% end
% 
% % loop over indeces
% time_points = nan(size(idx));
% for i = 1:size(idx,1)
%     time_points(i,:) = time_array(idx(i,:));
% end