function par = km_gettrlpar(trl,vars,param)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

if strcmpi(param,'offset')
    ib = find(ismember(vars,'trlbeg'),1,'last');
    io = find(ismember(vars,'trloff'),1,'last');
    par = trl(:,io) - trl(:,ib);
    return
end

i = find(ismember(vars,param),1,'last');
par = trl(:,i);