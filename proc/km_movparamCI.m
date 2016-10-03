function cfg = km_movparamCI(cfg,data)
%--------------------------------------------------------------------------
%
% See also KM_MOVPARAM
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end

% get trial variables and variable names
if isfield(cfg,'trl')
    trl = cfg.trl;
else
    trl = [];
end
if isfield(cfg,'vars')
    vars = cfg.vars;
else
    vars = {};
end

% select parameters to run trlrej on
mpcfg = cfg.movparamCI;
if isfield(cfg,'movparam') && isfield(cfg.movparam,'pname')
    pnamelist = cfg.movparam.pname;
else
    pnamelist = cfg.vars;
end
paramsel = match_str(pnamelist,km_labelselection(mpcfg.param,pnamelist));
mpcfg.param = pnamelist(paramsel);

% get indeces of trials to reject if requested
if ~isfalse(mpcfg.trlrej)
    cfg = km_setcfg(cfg,'sess');
    fname_idx = fullfile(cfg.dir.ana,sprintf('%s%s%s',cfg.subj,cfg.sess,mpcfg.trlrej.dataset));
    idxr = load(fname_idx);
    % reject trials listwise
    listwisesel = match_str(pnamelist,km_labelselection(mpcfg.trlrej.listwise,pnamelist));
    mpcfg.trlrej.listwise = pnamelist(listwisesel);
    idxr = km_combidxr(idxr.cfg,mpcfg.trlrej.listwise);
    % select pairwise
    pairwisesel = match_str(pnamelist,km_labelselection(mpcfg.trlrej.pairwise,pnamelist));
    mpcfg.trlrej.pairwise = pnamelist(pairwisesel);
else
    idxr = [];
end

% loop over factors
if ~iscell(mpcfg.fact), mpcfg.fact = {mpcfg.fact};	end
for f = 1:length(mpcfg.fact)
    
    % get factor(s)
    fact = mpcfg.fact{f};
    fname = mpcfg.fname{f};
    
    % get levels
    factlvl = cell(1,length(fact));
    for ff = 1:length(fact)
        factlvl{ff} = mpcfg.lvl.(fact{ff});
    end
    
    % combine levels of all factors
    lvl = arraycomb(factlvl{:});
    
    % determine dimensions
    dims = cellfun(@(x) length(x),factlvl);
    if length(dims) == 1
        dims = [dims 1];
    end
    
    % loop over parameters
    for p = 1:length(mpcfg.param)
        
        % obtain parameter values
        pname = lower(mpcfg.param{p});
        par = km_gettrlpar(trl,vars,pname);
        
        % reject trials listwise if requested
        if ~isempty(idxr.reject)
            par(idxr.reject) = NaN;
        end
        
        % reject trials pairwise if requested
        if ismember(pname,mpcfg.trlrej.pairwise) && isfield(idxr,pname)
            par(idxr.(pname)) = NaN;
        end
        
        % initialize structure
        tile = repmat({[]},prod(dims),1);
        CI = struct('mean',tile,'var',tile,'CI',tile,'ub',tile,'lb',tile);
        
        % devide end positions over factor x level combinations and calculate
        % multi-dimensional confidence intervals
        [~,iv] = ismember(fact,vars);
        for i = 1:size(lvl,1)
            idx = all(trl(:,iv) == repmat(lvl(i,:),size(trl,1),1),2);
            %CI(i) = km_CI(par(idx));
            CI(i) = orderfields(km_CI(par(idx)),CI);
        end
        
        % store end point error information in configuration structure
        cfg.movparamCI.(fname).(pname) = CI;
        
    end
    
    % store levels in configuration structure
    cfg.movparam.(fname).lvl	= lvl;
    
end


% update dataset
if isempty(regexpi(cfg.dataset,'_ana$','once'))
    if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
    cfg.dataset = [cfg.dataset '_ANA'];
end

% update processing directory
cfg = km_setcfg(cfg,'dirana');
cfg.dir.proc = cfg.dir.ana;
