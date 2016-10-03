function Cc = cellcont2str(C,flg)
%--------------------------------------------------------------------------
% CELLCONT2STR describes cell array contents in a string.
% FORMAT:   Cc      = cellcont2str(C,flg);
% INPUT:    C       - a cell array
%           flg     - flg determining if the output is a cell array of
%                     strings (default: 'cell'), or a concatenated string
%                     ('str').
%
% Copyright (C) 2010-2011, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2011-09-01
%--------------------------------------------------------------------------


%% Check input
%----------------------------------------
error(nargchk(1,2,nargin));
if ~iscell(C)
    error('LV:arg2str:BadArg','The input variable arg must be a cell array.');
end
%if length(C) ~= numel(C)
%    error('LV:arg2str:BadArg','The input variable arg must be a 1 x n cell array.');
%end
%C = C(:)';
if nargin < 2,	flg = 'cell';   end
flg = lower(flg);
if ~ismember(flg,{'cell','str','char','chararray'})
    error('LV:arg2str:BadFlag','The input variable flg was not recognized.');
end

% determine to print loose or tight
if strcmp(get(0,'formatspacing'),'loose');
    n = 24;
else
    n = 12;
end

% determine to add spaces
str_out = ismember(flg,{'str','char','chararray'});

% loop over cell contents and describe cell contents in a string
Cc = cell(size(C));
for i = 1:numel(C)
    describe = false;
    if ischar(C{i})
        if length(C{i}) > (n-1)
            describe = true;
        else
            Cc{i} = C{i};
        end
    elseif isnumeric(C{i}) || islogical(C{i})
        Cc{i} = num2str(C{i});
        if length(Cc{i}) > (n-1)
            describe = true;
        end
    elseif isa(C{i},'function_handle')
        Cc{i} = sprintf('@%s',func2str(C{i}));
    else
        describe = true;
    end
    % describe class
    if describe
        Cc{i} = sprintf('%0.f-by-%0.f %s',m,n,class(C{i}));
        if length(Cc{i}) > (n-1)
            Cc{i} = sprintf('%s',class(C{i}));
        end
    end
    % cut off if too long
    if length(Cc{i}) > (n-1)
        Cc{i} = Cc{i}(1:(n-1));
    end
    % add spaces (if requested)
    if str_out
        Cc{i} = [Cc{i} repmat(' ',1,max(0,n-length(Cc{i})))];
    end
end

% convert cell of character arrays to a single character array if requested
if str_out
    Cc = cell2mat(Cc);
end
