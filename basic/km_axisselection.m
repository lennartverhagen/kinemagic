function axissel = km_axisselection(axis,dataaxis)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------


if ischar(axis)
    axissel = num2cell(axis);
elseif ~iscell(appaxis)
    axissel = dataaxis(axis);
else
    axissel = axis;
end