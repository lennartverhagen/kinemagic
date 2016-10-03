function datmat = km_datcell2mat(datcell,varargin)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

if isempty(datcell)
    datmat = [];
    return;
end

if ~iscell(datcell)
    datmat = datcell;
    return;
end

% check for consistency of the number of dimensions of data matrices
ndmat = unique(cellfun(@ndims,datcell));
if length(ndmat) > 1
    error('KM:DatCell2Mat:ndims','The number of dimensions is not equal in all data cells');
end

% check for consistency of the size of the data matrices
siz = cellfun(@size,datcell,'UniformOutput',false);
siz = unique(cat(1,siz{:}),'rows');
if size(siz,1) > 1
    error('KM:DatCell2Mat:size','The size of the matrix is not equal in all data cells');
end

% check if siz matches the given dimensions
if nargin > 1
    sizgiven = [varargin{:}];
    if ~isequal(sizgiven,siz(1:length(sizgiven)))
        error('KM:DatCell2Mat:GivenSize','The size of the data does not match the entered dimensions');
    end
end

% check number of dimensions of the data cell
sizcell = size(datcell);
if sum(sizcell>1)<2
    ndcell = 1;
    datcell = datcell(:);
else
    ndcell = length(sizcell);
end

% prepad dimensions to allow cell2mat
sizmat = [ones(1,ndcell) siz];
datcell = cellfun(@(x) reshape(x,sizmat),datcell,'UniformOutput',false);

% convert cell to mat
datmat = cell2mat(datcell);