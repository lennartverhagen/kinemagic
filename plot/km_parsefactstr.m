function fact = km_parsefactstr(str)
%--------------------------------------------------------------------------
% parse the factor string, used by KM_PLOT
%
% See also KM_PLOT
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% make str a cell structure
if iscell(str)
    tstr = str;
else
    tstr = {str};
end

% loop over strings if str is a cell structure
fact = cell(1,length(tstr));
for f = 1:length(tstr)
    
    % get factor string
    factstr = tstr{f};
    
    % continue if fact str is already in the right format
    if iscell(factstr)
        fact{f} = factstr;
        continue
    end
    
    % remove spaces
    factstr = regexprep(factstr,'\s','');
    
    % split factors
    fact{f} = regexp(factstr,'x','split');
    
end
