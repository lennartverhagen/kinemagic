function varargout = km_idx2logic(varargin)
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
% function logic_array = km_idx2logic(idx,nsamp)
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
% logic_array = false(1,nsamp);
% if isempty(idx)
%     return
% elseif size(idx,2) == 1
%     logic_array(idx) = true;
% elseif size(idx,2) == 2
%     for i = 1:size(idx,1)
%         logic_array(idx(i,1):idx(i,2)) = true;
%     end
% else
%     error('Unsupported indeces matrix size. Size should be nx1 or nx2.')
% end