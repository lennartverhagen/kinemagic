function mov = km_setequalnmov(mov)
%--------------------------------------------------------------------------
% make sure that all trials have the same number of movements, fill with
% NaNs if needed
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

nmov = cellfun(@(x) size(x,1),mov');
idx = find(nmov<max(nmov));
for i = 1:length(idx)
    n = size(mov{idx(i)},1);
    mov{idx(i)} = [mov{idx(i)}; nan(max(nmov)-n,2)];
end