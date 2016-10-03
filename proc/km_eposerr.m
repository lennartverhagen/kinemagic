function [cfg,data] = km_eposerr(cfg,data)
%--------------------------------------------------------------------------
%
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
    error('no trl provided');
end
if isfield(cfg,'vars')
    vars = cfg.vars;
else
    error('no vars provided');
end

% get movement
if isfield(cfg.eposerr,'movname');
    movname = cfg.eposerr.movname;
else
    movname = 'movment';
end
if isfield(cfg,movname) && ~isstruct(cfg.(movname)) && ~isempty(cfg.(movname))
    mov = cfg.(movname);
elseif isfield(cfg,'movdef') && isfield(cfg.movdef,movname)
    mov = cfg.movdef.(movname);
elseif isfield(data,movname)
    mov = data.(movname);
else
    error('no %s provided',movname);
end

% make sure that all trials have the same number of movements
mov = km_setequalnmov(mov);

% get end point position data
tcfg = [];
tcfg.pos.movnr	= cfg.eposerr.movnr;
tcfg.pos.marker = cfg.eposerr.marker;
tcfg.pos.axis   = km_axisselection(cfg.eposerr.axis,data.axis);
endoffset = repmat({cfg.eposerr.endoffset},size(mov));
endmovement = cellfun(@(x,y) [x(:,2)+y x(:,2)+y],mov,endoffset,'UniformOutput',false);
pos = km_getmovdat(tcfg,data,endmovement,'pos');
dat = km_datcell2mat(pos);

% reject trials if requested
if ~isfalse(cfg.eposerr.trlrej)
    idxr = km_getidxr(cfg,cfg.eposerr);
    dat(idxr.reject,:,:) = NaN;
end

% loop over factors
nmark = length(cfg.eposerr.marker);
for f = 1:length(cfg.eposerr.fact)
    flg_posdef = 'posdef';
    
    % get factor(s)
    fact = cfg.eposerr.fact{f};
    fname = cfg.eposerr.fname{f};
    
    % get levels
    factlvl = cell(1,length(fact));
    for ff = 1:length(fact)
        factlvl{ff} = cfg.eposerr.lvl.(fact{ff});
    end
    
    % combine levels of all factors
    lvl = arraycomb(factlvl{:});
    
    % determine dimensions
    dims = cellfun(@(x) length(x),factlvl);
    if length(dims) == 1
        dims = [dims 1];
    end
        
    % loop over markers
    for m = 1:nmark
        mname = cfg.eposerr.marker{m};
        
        % initialize structure
        tile = repmat({[]},prod(dims),1);
        naxes = size(dat,3);
        if naxes ==1
            tmpCI = struct('x',tile);
        elseif naxes ==2
            tmpCI = struct('x',tile,'y',tile,'xy',tile);
        elseif naxes == 3
            tmpCI = struct('x',tile,'y',tile,'z',tile,'xy',tile,'yz',tile,'zx',tile,'xyz',tile);
        else
            tmpCI = struct;
        end
        
        % devide end positions over factor x level combinations and calculate
        % multi-dimensional confidence intervals
        [~,iv] = ismember(fact,vars);
        for i = 1:size(lvl,1)
            idx = all(trl(:,iv) == repmat(lvl(i,:),size(trl,1),1),2);
            if ~any(idx), warning('KM:EPosErr:notrials','No trials were selected'); end
            epos = dat(idx,m,:);
            epos = reshape(epos,size(epos,1),size(epos,3));
            % test if enough data points are present
            if sum(~any(isnan(epos),2)) <= size(epos,2)
                warning('KM:eposerr:NoData','In subject: %s%s, factor: %s, less data points present than the number of dimensions: multi-dimensional covariance decomposition is not possible',cfg.subj,cfg.sess,fname);
                flg_posdef = 'nonposdef';
            end
            epos_CI = km_mdCI(epos,flg_posdef); epos_CI = rmfield(epos_CI,{'conf','kdist','k','mu'});
            tmpCI(i) = orderfields(epos_CI,tmpCI);
        end
        
        % reorganize CI
        CI = struct;
        axisname = fieldnames(tmpCI);
        for an = 1:length(axisname)
            flds = fieldnames(tmpCI(i).(axisname{an}));
            paramflds = flds(~ismember(flds,{'eigvec','eigval','xyz','d','idx_out'}));
            for pfn = 1:length(paramflds)
                tmp = nan(size(lvl,1),1);
                for i = 1:size(lvl,1)
                    tmp(i) = tmpCI(i).(axisname{an}).(paramflds{pfn});
                end
                CI.(axisname{an}).(paramflds{pfn}) = reshape(tmp,dims);
            end
            cellflds = flds(ismember(flds,{'eigvec','eigval','xyz','d','idx_out'}));
            for xfn = 1:length(cellflds)
                tmp = cell(size(lvl,1),1);
                for i = 1:size(lvl,1)
                    tmp{i} = tmpCI(i).(axisname{an}).(cellflds{xfn});
                end
                CI.(axisname{an}).(cellflds{xfn}) = reshape(tmp,dims);
            end
        end
        
        % store end point error information in data structure
        % cfg.eposerr.(fname).(mname) = CI; % outdated
        data.eposerr.(fname).(mname) = CI;
        
    end
    
    % store levels in data structure
    % cfg.eposerr.(fname).lvl = lvl; % outdated
    data.eposerr.(fname).lvl = lvl;
    
end

% evaluate experimentally specific function
if isfield(cfg.eposerr,'expfun') && ~isfalse(cfg.eposerr.expfun)
    % cfg.eposerr = feval(cfg.eposerr.expfun,cfg.eposerr); % outdated
    data.eposerr = feval(cfg.eposerr.expfun,data.eposerr);
end


% update dataset
if isempty(regexpi(cfg.dataset,'_ana$','once'))
    if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
    cfg.dataset = [cfg.dataset '_ANA'];
end

% update processing directory
cfg = km_setcfg(cfg,'dirana');
cfg.dir.proc = cfg.dir.ana;
%--------------------------------------------------------------------------
