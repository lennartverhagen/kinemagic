function idx = km_combidxr(idx,rejectparam)
%--------------------------------------------------------------------------
% combine the trial rejection as determined previously
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% use only reject parameters that are fields of idx
rejectparam = rejectparam(ismember(rejectparam,fieldnames(idx)));

% return if nothing should be rejected
if isempty(rejectparam)
    if isfield(idx,'zeros')
        idx.reject = false(size(idx.zeros));
    elseif isfield(idx,'apriori')
        idx.reject = false(size(idx.apriori));
    elseif isfield(idx,'nan')
        idx.reject = false(size(idx.nan));
    elseif isfield(idx,'behav')
        idx.reject = false(size(idx.behav));
    else
        error('number of trials unknown');
    end

% initialize idx.reject based on first rejection criterium
else
    idx.reject = getsubfield(idx,rejectparam{1});
    %idx.reject = idx.(rejectparam{1});
end

% continue selecting trials based on remaining rejection criteria
for i = 2:length(rejectparam)
    
    idx.reject = idx.reject | getsubfield(idx,rejectparam{i});
    %idx.reject = idx.reject | idx.(rejectparam{i});
    
end