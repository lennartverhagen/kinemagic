function tim = km_gettim(dattim,trl,vars,idx)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

if isstruct(dattim),	dattim = dattim.time;	end

idx_nan = isnan(idx);
idx(idx_nan) = 1;

off = km_gettrlpar(trl,vars,'offset');

ntrl = length(dattim);
if size(idx,2) == ntrl, idx = idx'; end


dims = size(idx);
tim = nan(dims);
for i = 1:prod(dims(2:end))
    tim(:,i) = cellfun(@(x,j) x(j),dattim',num2cell(idx(:,i))) + off;
end

tim(idx_nan) = NaN;