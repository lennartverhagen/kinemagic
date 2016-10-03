function idxr = km_getidxr(cfg,datacfg)
%--------------------------------------------------------------------------
% get reject trial indeces
%
% See also KM_PLOT
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% check input
if nargin > 1
    % datacfg and cfg are combined, with cfg fields overwriting datacfg
    % fields in case of overlap
    cfg = combstruct(datacfg,cfg);
end

% get trial rejection indeces from data configuration
if isfield(cfg,'idxr')
    idxr = cfg.idxr;
    
% load trial rejection indeces
elseif ~isfalse(cfg.trlrej)
    % load IDXR
    cfg = km_setcfg(cfg,'sess');
    fname = fullfile(cfg.dir.proc,sprintf('%s%s%s',cfg.subj,cfg.sess,cfg.trlrej.dataset));
    idx = load(fname);
    
    % select parameters to run trlrej on
    listwise = cfg.trlrej.listwise;
    %if isfield(cfg,'movparam') && isfield(cfg.movparam,'pname')
    %    pnamelist = cfg.movparam.pname;
    %elseif isfield(cfg,'vars')
    %    pnamelist = cfg.vars;
    %else
    %    pnamelist = cfg.trlrej.listwise;
    %end
    %if ~ismember('apriori',pnamelist), pnamelist = [{'apriori'} pnamelist];
    %listwisesel = match_str(pnamelist,km_labelselection(cfg.trlrej.listwise,pnamelist));
    %listwise = pnamelist(listwisesel);
    
    idxr = km_combidxr(idx.cfg,listwise); clear idx;

% do not reject any trials
else
    nrpt = size(cfg.trl,1);
    idxr.reject = false(nrpt,1);
end
