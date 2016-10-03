function idx = km_trlrej(cfg,data)
%--------------------------------------------------------------------------
% select trials to be rejected. Automated or by hand
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('FTW:Return','FTW %s: Nothing to do...',task);
    return
end


%% Initiate index structure
%----------------------------------------
% initiate trial indeces
if isfield(data,'time')
    ntrl = length(data.time);
else
    ntrl = size(cfg.trl,1);
end
idx.zeros = false(ntrl,1);
% select parameters to run trlrej on
if isfield(cfg,'vars');
    pnamelist = cfg.vars;
elseif isfield(cfg,'movparam') && isfield(cfg.movparam,'pname')
    pnamelist = cfg.movparam.pname;
else
    pnamelist = cfg.trlrej.pname;
end
aprioriparamsel = match_str(pnamelist,km_labelselection(cfg.trlrej.apriori.param,pnamelist));
cfg.trlrej.apriori.param = pnamelist(aprioriparamsel);
paramsel = match_str(pnamelist,km_labelselection(cfg.trlrej.param,pnamelist));
cfg.trlrej.param = pnamelist(paramsel);
% recheck/reset cfg
cfg = km_setcfg(cfg,task);

%% Apriori trial rejection
%----------------------------------------
% reject trials apriori based on a user specified function
if ~isfalse(cfg.trlrej.apriori.fun)
	idx = feval(cfg.trlrej.apriori.fun,cfg,data,idx);
else
    idx.apriori = false(ntrl,1);
end
% reject trials apriori based on parameter cutoffs
if ~isfalse(cfg.trlrej.apriori.param)
    for p = 1:length(cfg.trlrej.apriori.param)
        % param name
        pname = cfg.trlrej.apriori.param{p};
        if strcmpi(pname,'apriori'),	continue;	end
        if ~any(strcmp(cfg.vars,pname))
            warning('KM:NotSupported','KM %s: Parameter ''%s'' is not supported.',mfilename,pname);
        	continue
        end
        % get parameter and apply cutoffs
        dat = cfg.trl(:,strcmp(cfg.vars,pname));
        if ~isempty(regexp(cfg.trlrej.apriori.fun,'dat','once'))
            idx.apriori = idx.apriori || eval(cfg.trlrej.apriori.fun);
        else
            [~, idx.apriori] = cutoutlier(dat,NaN,'none',cfg.trlrej.apriori.cutmode,cfg.trlrej.(pname),idx.apriori);
        end
    end
end

%% Manual trial rejection
%----------------------------------------
if ~strcmpi(cfg.trlrej.fun,'auto')
    
    % evaluate function
	idx = feval(cfg.trlrej.fun,cfg,data,idx);
    
    
%% Automatic trial rejection
%----------------------------------------
else
    
    % get trial data and variable names
    vars = cfg.vars;
    trl  = cfg.trl;
    
    % loop over automatic trial rejection parameters
    for p = 1:length(cfg.trlrej.param)
        
        % param name
        pname = cfg.trlrej.param{p};
        if strcmpi(pname,'apriori'),	continue;	end
        if ~any(strcmp(vars,pname))
            warning('KM:NotSupported','KM %s: Parameter ''%s'' is not supported.',mfilename,pname);
        	continue
        end
        
        % get parameter
        dat = trl(:,strcmp(vars,pname));
        
        % get initial cutoffs
        if isfield(cfg.trlrej,pname)
            initcut = cfg.trlrej.(pname);
        else
            initcut = [-Inf Inf];
        end
        
        % select trials with too short/low or long/high parameter values
        [idx.describ.(pname).low, idxextrlow, cutlow] 	= cutoutlier(dat,3,'iqr','lower',initcut,idx.apriori);
        [idx.describ.(pname).high, idxextrhigh, cuthigh]= cutoutlier(dat,3,'iqr','higher',initcut,idx.apriori);
        idx.describ.(pname).cutlow = cutlow;
        idx.describ.(pname).cuthigh = cuthigh;
        idxextreme = idxextrlow | idxextrhigh;
        
        % calculate descriptive metrics
        step = (cuthigh-cutlow)/25;
        edges = (cutlow:step:(cuthigh+step)) - step/2;
        [~,bin] = histc(dat,edges); bin(bin==0) = NaN; modbin = mode(bin);
        if ~isnan(modbin),	modbin = edges(modbin) + step/2;	end
        idx.describ.(pname).all.values	= dat;
        idx.describ.(pname).all.mean    = nanmean(dat);
        idx.describ.(pname).all.median  = nanmedian(dat);
        idx.describ.(pname).all.mode    = modbin;
        [~,bin] = histc(dat(~idxextreme),edges); bin(bin==0) = NaN; modbin = mode(bin);
        if ~isnan(modbin),	modbin = edges(modbin) + step/2;	end
        idx.describ.(pname).sel.values  = dat(~idxextreme);
        idx.describ.(pname).sel.mean    = nanmean(dat(~idxextreme));
        idx.describ.(pname).sel.median  = nanmedian(dat(~idxextreme));
        idx.describ.(pname).sel.mode    = modbin;
        
        % combine low-high rejections
        idx.(pname) = idx.describ.(pname).low | idx.describ.(pname).high;
        
    end
    
end

% report trial indeces
cfg.trlrej.idx = idx;
idx.report = km_trlrej_report(cfg);

% only the idx needs to be saved
if ~isfalse(cfg.save)
    idx.save = 'cfg';
else
    idx.save = 'no';
end
idx.subj = cfg.subj;
idx.sess = cfg.sess;
idx.dir  = cfg.dir;
idx.dir.proc = cfg.dir.ana;

% update dataset
idx.dataset = [cfg.dataset '_IDXR'];
%--------------------------------------------------------------------------
