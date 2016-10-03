function [cfg,Sparam,SCI,Sepe,Stseries] = km_plot_getdata(cfg)
%--------------------------------------------------------------------------
% get kinematic data (parameters and time-course)
%
% See also KM_PLOT
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% condition, subject and session names
nsubj = length(cfg.subj);
nsess = length(cfg.sess);
bfact  = cfg.plot.bfact;
nfact  = cfg.plot.nfact;
ntseries = length(cfg.plot.tseries);
if ~isfalse(cfg.plot.tseries)
    bcname = cfg.plot.bcname;
    nbcond = length(bcname);
end

% load group analysis dataset
if nsubj == 1 && ~isempty(strfind(lower(cfg.subj{1}),'group'))
    error('Preanalyzed group-analyses should not be reconstructed, but loaded directly. Use cfg.load = true;');
end


%% collect data
%----------------------------------------
% initiate progress bar
if nsubj > 1
    hsubj = waitbar(0,'Processing individual subjects...');
end

% loop over subjects
Sparam = cell(nsubj,1);
SCI = cell(nsubj,1);
Sepe = cell(nsubj,1);
Stseries = cell(nsubj,1);
for s = 1:nsubj
    
    % define subject analysis directory
    if ~isfalse(cfg.plot.fetchdatafromcfg)
        dir_proc = 'fetchdatafromcfg';
    else
        dir_proc = getsubjsubdir(cfg,s,'ana');
%     if ~isfield(cfg.dir,'ana') || isempty(cfg.dir.ana)
%         dir_proc = fullfile(cfg.dir.root,cfg.subj{s},'kinematics','ana');
%     elseif ~isempty(regexp(cfg.dir.ana,['^' regexptranslate('escape',cfg.dir.root)],'once'))
%         dir_proc = cfg.dir.ana;
%     elseif isempty(fileparts(cfg.dir.ana)) || ~exist(cfg.dir.ana,'dir')
%         dir_proc = fullfile(cfg.dir.root,cfg.subj{s},cfg.dir.ana);
%     else
%         dir_proc = cfg.dir.ana;
%     end
    end
    
    % initialize data structure
    data = struct(...
        'dcfg',{},...
        'param',{},...
        'CI',{},...
        'epe',{},...
        'tseries',{},...
        'idxc',{},...
        'idxr',{},...
        'idxe',{}...
        );
    
    % add requested arbitrary fields
    if isfield(cfg.plot,'fetchdata') && ~isfalse(cfg.plot.fetchdata)
        for f = 1:length(cfg.plot.fetchdata)
            data().(cfg.plot.fetchdata{f}) = [];
        end
    end
    
    %% get session data
    %----------------------------------------
    % initiate progress bar
    if s == 1 && nsess > 1
        hproc = waitbar(0,'Loading datasets...');
        % find all waitbar handles, use the current one
        h = findobj(allchild(0),'flat','Tag','TMWWaitbar');
        if length(h) > 1
            % calculate new position
            pos = get(h(2),'Position');
            % shift current waitbar one step down
            pointsPerPixel = 72/get(0,'ScreenPixelsPerInch');
            pos(2) = pos(2)-pos(4)-28*pointsPerPixel;
            set(h(1),'Position',pos);
        end
    end
    
    % loop over sessions
    for ss = 1:nsess
        
        % load data
        data(ss) = load_data(cfg,s,ss,dir_proc);
        
        % reject outliers based on km_trlrej
        data(ss) = reject_outliers(cfg,data(ss));
        
        % exclude conditions based on cfg.plot.factexcl
        data(ss) = exclude_cond(cfg,data(ss));
        
        % normalize tseries in time and concatenate trials
        data(ss).tseries = tseries_norm(cfg,data(ss).tseries);
        
        % update progress bar
        if nsess > 1
            waitbar(ss/nsess,hproc);
        end
        
    end
    
    % update the progressbar
    if nsess > 1
        waitbar(ss/nsess,hproc,'Merge datasets...');
    end
    
    %% combine sessions and select conditions and rejections
    %----------------------------------------
    % combine sessions if requested
    dcfg    = km_combsess({data(:).dcfg},'cfg');
    param	= km_combsess({data(:).param},'param');
    CI      = km_combsess({data(:).CI},'CI');
    epe     = km_combsess({data(:).epe},'epe');
    tseries = km_combsess({data(:).tseries},'tseries');
    idxc    = km_combsess({data(:).idxc},'idxc');
    idxr    = km_combsess({data(:).idxr},'idx');
    idxe    = km_combsess({data(:).idxe},'idx');
    
    % store the trial matrix in the configuration structure
    if s == 1,  cfg.vars = dcfg.vars;   end
    cfg.trl{s} = dcfg.trl;
    cfg.idxr{s} = idxr;
    cfg.idxe{s} = idxe;
    
    % fetch a designated field if requested
    if isfield(cfg.plot,'fetchcfg') && ~isfalse(cfg.plot.fetchcfg)
        for f = 1:length(cfg.plot.fetchcfg)
            fld = cfg.plot.fetchcfg{f};
            cfg.(fld){s} = dcfg.(fld);
        end
    end
    
    % dirty HACK to get traject from data
    for ss = 1:nsess
        if isfield(data(ss),'traject')
            cfg.traject{s}{ss} = data(ss).traject;
        end
    end
    
    % rename session factor
    [CI,epe] = rename_sessfact(CI,epe);
    
    if nsess > 1
        waitbar(0,hproc,'Processing factors...');
    end
    
    % select conditions specific parameters
    for f = 1:nfact
        % get factor(s)
        fact = cfg.plot.fact{f};
        fname = cfg.plot.fname{f};
        
        % get levels
        factlvl = cell(1,length(fact));
        for ff = 1:length(fact)
            factlvl{ff} = cfg.plot.lvl.(fact{ff});
        end
        
        % combine levels of all factors
        lvl = arraycomb(factlvl{:});
        
        % select params repetitions
        tmpparam = cell(size(lvl,1),1);
        tmptseries = struct;
        for i = 1:size(lvl,1)
            [~,idxcond] = ismember(fact,bfact);
            A = idxc(:,idxcond);
            B = repmat(lvl(i,:),size(idxc,1),1);
            idx = all(A==B,2);
            if ~isempty(param)
                tmpparam{i} = param(idx,:);
            end
            for ts = 1:ntseries
                tsname = cfg.plot.tseries{ts};
                switch ndims(tseries.(tsname))
                    case 4,     tmptseries.(tsname){i} = tseries.(tsname)(idx,:,:,:);
                    case 3,     tmptseries.(tsname){i} = tseries.(tsname)(idx,:,:);
                    otherwise,  tmptseries.(tsname){i} = tseries.(tsname)(idx,:);
                end
            end
        end
        
        % store subject params
        dims = cellfun(@(x) length(x),factlvl);
        if length(dims) == 1, dims = [1 dims];  end
        if ~isempty(param)
            Sparam{s}.(fname) = reshape(tmpparam,dims);
        end
        for ts = 1:ntseries
            tsname = cfg.plot.tseries{ts};
            Stseries{s}.(tsname).(fname) = reshape(tmptseries.(tsname),dims);
            Stseries{s}.(tsname).dimord = tseries.dimord.(tsname);
            Stseries{s}.(tsname).axis   = tseries.axis.(tsname);
            Stseries{s}.(tsname).marker = tseries.marker.(tsname);
        end
        
        % find CI factor
        % get CI data if requested
        if ~isempty(CI)
            tmpfact = CI.fact{cellfun(@(x) all(ismember(fact,x)),CI.fact)};
            tmpfname = sprintf('%sx',tmpfact{:}); tmpfname = tmpfname(1:end-1);
            tmpparam = regexprep(cfg.plot.CI.param,'CI$','');
            if isfield(CI,fname)
                for p = 1:length(tmpparam)
                    SCI{s}.(fname).(tmpparam{p}) = [CI.(fname).(tmpparam{p}).CI];
                end
            elseif isfield(CI,tmpfname)
                [~,J] = ismember(fact,tmpfact);
                for p = 1:length(tmpparam)
                    tmpCI = [CI.(tmpfname).(tmpparam{p}).CI];
                    SCI{s}.(fname).(tmpparam{p}) = permute(tmpCI,J);
                end
            end
        end
        
        % find epe factor
        % get epe data if requested
        if ~isempty(epe)
            tmpfact = epe.fact{cellfun(@(x) all(ismember(fact,x)),epe.fact)};
            tmpfname = sprintf('%sx',tmpfact{:}); tmpfname = tmpfname(1:end-1);
            marker = cfg.plot.epe.marker;
            axisname = cfg.plot.epe.axis;
            paramname = cfg.plot.epe.param;
            b = cfg.plot.epe.scale;
            if isfield(epe,fname)
                for m = 1:length(marker)
                    Sepe{s}.(fname).(marker{m}) = b * epe.(fname).(marker{m}).(axisname).(paramname);
                end
            elseif isfield(epe,tmpfname)
                [~,J] = ismember(fact,tmpfact);
                for m = 1:length(marker)
                    tmpepe = b * epe.(tmpfname).(marker{m}).(axisname).(paramname);
                    Sepe{s}.(fname).(marker{m}) = permute(tmpepe,J);
                end
            end
            if isempty(Sepe{s}), warning('KM:plot:NoEPEData','No end-point-error data is found (maybe you did not specify the correct factor).'); end
        end
        
        % update progress bar
        if nsess > 1
            waitbar(f/nfact,hproc);
        end
        
    end
    
    % select condition specific time series
    if ~isfalse(ntseries)
        %% FIXME: tseries are now handled as params, but this does not
        %% allow substraction and contrasts...
%         % get trial indeces for the requested conditions
%         idxc = getcondidx(cfg.plot.bcond,idxc);
%         
%         % FIXME: adjust code to tseries and param
%         %% select condition trials
%         %----------------------------------------
%         for c = 1:nbcond
%             if isfield(Stseries{s},bcname{c}),	continue;	end
%             
%             % select trials
%             idxcond = idxc.(bcname{c});
%             trials	= find(idxcond & ~idxr.reject);
%             
%             % % do time-lock baseline analyses
%             % cfg.tlana.trials = trials;
%             % [cfg,cdata] = ftw_tlana(cfg,tseries);
%             
%             % Stseries{s}.(bcname{c})	= cdata;
%             
%             Stseries{s}.(bcname{c})	= cdata;
%         end
    end
    
    % update progress bar
    if nsubj > 1
        waitbar(s/nsubj,hsubj);
    end
        
end

% now it is time to clean up some memory
clear tseries param;

% close the progressbar
if nsess > 1
    close(hproc);
end

% close the progressbar
if nsubj > 1
    close(hsubj);
end


%% remove conditions with too few repetitions
%----------------------------------------

if cfg.plot.reptcontrib
    % loop over subjects
    for s = 1:nsubj
        
        if ~isempty(Sparam{s})
            % loop over factors
            for f = 1:nfact
                dat = Sparam{s}.(cfg.plot.fname{f});
                for i = 1:numel(dat);
                    for p = 1:size(dat{i},2)
                        % keep only data with enough repetitions contributing
                        if sum(~isnan(dat{i}(:,p)),1) < cfg.plot.reptcontrib
                            Sparam{s}.(cfg.plot.fname{f}){i}(:,p) = NaN;
                        end
                    end
                end
            end
        end
        
        % loop over time-series
        for ts = 1:ntseries
            tsname = cfg.plot.tseries{ts};
            % loop over factors
            for f = 1:nfact
                fname = cfg.plot.fname{f};
                dat = Stseries{s}.(tsname).(fname);
                for i = 1:numel(dat);
                    % keep only data with enough repetitions contributing
                    dims = size(dat{i});
                    switch length(dims)
                        case 4
                            for j2 = 1:dims(2)
                                for j3 = 1:dims(3)
                                    if sum(~any(isnan(squeeze(dat{i}(:,j2,j3,:))),2)) < cfg.plot.reptcontrib
                                        dat{i}(:,j2,j3,:) = nan([dims(1) 1 1 dims(4)]);
                                    end
                                end
                            end
                        case 3
                            for j2 = 1:dims(2)
                                if sum(~any(isnan(squeeze(dat{i}(:,j2,:))),2)) < cfg.plot.reptcontrib
                                    dat{i}(:,j2,:) = nan([dims(1) 1 dims(2)]);
                                end
                            end
                        otherwise
                            if sum(~any(isnan(dat{i}),2)) < cfg.plot.reptcontrib
                                dat{i}(:,:) = nan(dims);
                            end
                            
                    end
                end
            end
            
            % TODO: adjust code to tseries
%             % loop over conditions
%             for c = 1:nbcond
%                 dat = Stseries{s}.(bcname{c}).powspctrm;
%                 nrept = size(dat,1);
%                 
%                 % TODO: check if the number of reptcontrib should really be
%                 % based on the number of NaNs. Shouldn't it be based on the
%                 % number of non-nans?
%                 dbstop
%                 
%                 % keep only data with enough repetitions contributing
%                 maxnrnans = nrept-cfg.plot.reptcontrib;
%                 idx_nan = repmat( sum(isnan(dat),1) > maxnrnans ,nrept,1);
%                 Stseries{s}.(bcname{c}).powspctrm(idx_nan) = NaN;
%             end
        end
        
    end
end


%% assign report directory
%----------------------------------------
if ~isfield(cfg.dir,'report')
    if nsubj > 1
        cfg.dir.report = fullfile(cfg.dir.root,'Group','kinematics','report');
    else
        cfg.dir.report = fullfile(cfg.dir.root,cfg.subj{1},'kinematics','report');
    end
end
%--------------------------------------------------------------------------


%% function load_data
%----------------------------------------
function data = load_data(cfg,s,ss,dir_proc)
ntseries = length(cfg.plot.tseries);

% allow functionality for names with and without a session string
sess = cfg.sess{ss};
if ~isempty(sess) && ~strcmpi(sess(1),'_')
    sess = ['_' sess];
end

% load configuration and parameters
if ~isfalse(cfg.plot.fetchdatafromcfg)
    dcfg = cfg.dcfg{s}{ss};
else
    tdcfg = load(fullfile(dir_proc,sprintf('%s%s%s_CFG',cfg.subj{s},sess,cfg.dataset)));
    dcfg = tdcfg.cfg;
end

% add variables to trl matrix if requested
% addvar.newvar	= 'newvar';
% addvar.oldvar	= {'oldvar1','oldvar2'};
% addvar.fun	= 'tmp(:,1)==1 & tmp(:,2)==2';
for i = 1:length(cfg.plot.addvar)
    idx = cellfun(@(x) find(strcmpi(dcfg.vars,x)),cfg.plot.addvar(i).oldvar);
    dcfg.vars = [dcfg.vars(1:max(idx)) {cfg.plot.addvar(i).newvar} dcfg.vars(max(idx)+1:end)];
    tmp = dcfg.trl(:,idx); % this variable will be used in the string to evaluate below
    tmp = eval(cfg.plot.addvar(i).fun);
    dcfg.trl = [dcfg.trl(:,1:max(idx)) tmp dcfg.trl(:,max(idx)+1:end)];
end

% check parameter availability
if ~isempty(cfg.plot.param) || ~isfalse(cfg.plot.epe) || ~isfalse(cfg.plot.CI)
    [tmp,idxvar] = ismember(cfg.plot.param,dcfg.vars);
    if ~all(tmp)
        tmpstr = sprintf('%s, ',cfg.plot.param{~tmp});
        tmpstr = tmpstr(1:end-2);
        error('not all parameters requested to plot are present in the dataset: %s',tmpstr);
    end
    param = dcfg.trl(:,idxvar);
else
    param = [];
end

% get trial condition indeces
[~,idxcond] = ismember(cfg.plot.bfact,dcfg.vars);
idxc = dcfg.trl(:,idxcond);
% get trial exclusion indeces
[~,idxexcl] = ismember(cfg.plot.factexcl,dcfg.vars);
idxe = dcfg.trl(:,idxexcl);

% load reject trial indeces
tcfg = [];
tcfg.subj       = cfg.subj{s};
tcfg.sess       = cfg.sess{ss};
tcfg.dir.proc   = dir_proc;
tcfg.trlrej     = cfg.plot.trlrej;
idxr        = km_getidxr(tcfg,dcfg);

% load parameter confidence interval data if requested
if ~isfalse(cfg.plot.CI)
    CI = dcfg.movparamCI;
else
    CI = [];
end

% load data
if ~isfalse(cfg.plot.tseries) || ~isfalse(cfg.plot.epe) || (isfield(cfg.plot,'fetchdata') && ~isfalse(cfg.plot.fetchdata))
    tdata = load(fullfile(dir_proc,sprintf('%s%s%s_DATA',cfg.subj{s},sess,cfg.dataset)));
end

% load end-point-error data if requested
if ~isfalse(cfg.plot.epe)
    if isfield(tdata.data,'eposerr')
        epe = tdata.data.eposerr;
    else
        epe = dcfg.eposerr;
    end
    if ~isfield(epe,'fact')
        epe.fact = cfg.plot.fact;
        epe.lvl = cfg.plot.lvl;
        epe.nfact = cfg.plot.nfact;
        epe.fname = cfg.plot.fname;
        epe.bfact = cfg.plot.bfact;
        epe.nbfact = cfg.plot.nbfact;
    end
else
    epe = [];
end

% process time series
if ~isfalse(cfg.plot.tseries)    
    % process data if requested
    if ~isfalse(cfg.tseries.proc) && ~isfalse(cfg.tseries.proc.task)
        tcfg            = cfg.tseries.proc;
        tcfg.subj       = cfg.subj{s};
        tcfg.sess       = cfg.sess{ss};
        tcfg.dataset	= cfg.dataset;
        [~,tdata.data]	= km_dotask(tcfg,tdata.data);
    end
    
    tcfg = cfg.tseries;
    for ts = 1:ntseries
        tsname = cfg.plot.tseries{ts};
        % get movement
        if isfield(cfg.tseries.(tsname),'movname');
            movname = cfg.tseries.(tsname).movname;
        else
            movname = 'movement';
        end
        if isfield(dcfg,movname) && ~isstruct(dcfg.(movname)) && ~isempty(dcfg.(movname))
            mov = dcfg.(movname);
        elseif isfield(dcfg,'movdef') && isfield(dcfg.movdef,movname)
            mov = dcfg.movdef.(movname);
        elseif isfield(tdata.data,movname)
            mov = tdata.data.(movname);
        else
            error('no %s provided for time-series [%s] plotting',movname,tsname);
        end
        
        % make sure that all trials have the same number of movements
        mov = km_setequalnmov(mov);
        
        % get time series
        tseries.(tsname) = km_getmovdat(tcfg,tdata.data,mov,tsname);
        
        % get time
        movnr = km_getmovidx(tcfg,tsname);
        tim = cellfun(@(x) x(movnr,:),mov,'UniformOutput',false);
        tseries.time = cell(1,length(tim));
        for t = 1:length(tim)
            it = km_time2logic(tim{t},tdata.data.time{t});
            tseries.time{t} = tdata.data.time{t}(it);
        end
        
        % get baseline window when object is held
        if any(strcmpi(tsname,{'pos','gripapt'})) && ~isfalse(cfg.tseries.(tsname).basecorr)
            basewin = cfg.tseries.(tsname).basecorr.win;
            if isa(basewin,'function_handle')
                basewin = cellfun(basewin,mov,'UniformOutput',false);
            elseif ~iscell(basewin)
                basewin = cellfun(@(x) basewin,mov,'UniformOutput',false);
            end
            tseries.basedata = km_getmovdat(tcfg,tdata.data,basewin,tsname);
            tseries.basecorr = cfg.tseries.(tsname).basecorr;
        end
        
        % store dimension order
        nd = ndims(tseries.(tsname){1});
        if nd == 2
            tseries.dimord.(tsname) = 'axis_time';
        elseif nd == 3
            tseries.dimord.(tsname) = 'sensor_axis_time';
        else
            tseries.dimord.(tsname) = 'unknown';
        end
        % store axis names
        tseries.axis.(tsname) = tcfg.(tsname).axis;
        tseries.marker.(tsname) = tcfg.(tsname).marker;
        
        % adjust sensor position from nail to skin
        if strcmpi(tsname,'pos') && ~isfalse(cfg.tseries.(tsname).basecorr)
            %tseries = basecorr_pos(cfg,tseries);
        end
        
        % adjust grip aperture to remove finger width
        if strcmpi(tsname,'gripapt') && ~isfalse(cfg.tseries.(tsname).basecorr)
            tseries = basecorr_gripapt(cfg,tseries);
        end
        
    end
else
    tseries = struct;
end

% collect output in structure
data.dcfg = dcfg;
data.param = param;
data.CI = CI;
data.epe = epe;
data.tseries = tseries;
data.idxc = idxc;
data.idxr = idxr;
data.idxe = idxe;

% fetch a designated field if requested
if isfield(cfg.plot,'fetchdata') && ~isfalse(cfg.plot.fetchdata)
    for f = 1:length(cfg.plot.fetchdata)
        fld = cfg.plot.fetchdata{f};
        data.(fld) = tdata.data.(fld);
    end
end

%--------------------------------------------------------------------------

%% function basecorr_pos
%----------------------------------------
function data = basecorr_pos(cfg,data)

descrip = data.basecorr.descrip;
if ~isa(descrip,'function_handle')
    if any(strcmpi(descrip,{'mean','median'}))
        descrip = ['nan' descrip];
    end
    descrip = str2func(lower(descrip));
end
% execute description function over the last dimension if appropriate
if ~isempty(regexp(func2str(descrip),'(mean|median|mode|sum|minus)$','once'))
    dim = repmat({length(size(data.basedata{1}))},size(data.basedata));
    val = cellfun(descrip,data.basedata,dim,'UniformOutput',false);
else
    val = cellfun(descrip,data.basedata,'UniformOutput',false);
end
corrval = data.basecorr.val;

switch lower(data.basecorr.corrflg)
    case {'sess','session'}
        % reshape val to get trials in the fourth dimension
        val = reshape([val{:}],[size(val{1}) length(val)]);
        % get descriptives for each sensors axis combination
        val = feval(descrip,val,length(size(val)));
        % FIXME: I assume that you would like to do correct the position to
        % on the basis of a fixed aperture, not position!
        if size(val,1) == 2
            % calculate the mean of the sensors
            m = repmat(mean(val,1),[2 1]);
            % calculate the apperture between the sensors
            app = sqrt(sum(diff(val,1).^2));
            % determine the new distance to the desired position based on
            % the ratio between the desired and observed aperture
            r = corrval/app;
            d = (r*val + (1-r)*m) - val;
            % Ouch, painful. I have been so smart right up to this point,
            % and now I realize that you need the orientation of each
            % individual sensor before you can correct the data!!!
        else
            error('FIXME:position can only be corrected on the basis of an aperture at the moment.');
        end
    case {'trl','trial'}
    otherwise
        error('Baseline correction flag [%s] is not supported',data.basecorr.corrflg);
end


%% function basecorr_gripapt
%----------------------------------------
function data = basecorr_gripapt(cfg,data)

descrip = data.basecorr.descrip;
if ~isa(descrip,'function_handle')
    if any(strcmpi(descrip,{'mean','median'}))
        descrip = ['nan' descrip];
    end
    descrip = str2func(lower(descrip));
end
% execute description function over the last dimension if appropriate
if ~isempty(regexp(func2str(descrip),'(mean|median|mode|sum|minus)$','once'))
    dim = repmat({length(size(data.basedata{1}))},size(data.basedata));
    val = cellfun(descrip,data.basedata,dim,'UniformOutput',false);
else
    val = cellfun(descrip,data.basedata,'UniformOutput',false);
end
corrval = data.basecorr.val;

switch lower(data.basecorr.corrflg)
    case {'sess','session'}
        % reshape val to get trials in the fourth dimension
        val = reshape([val{:}],[size(val{1}) length(val)]);
        % TODO: what if more than one dimension is provided?
        val = feval(descrip,val,length(size(val)));
        % correct data with base value
        data.gripapt = cellfun(@(x) bsxfun(@minus,x,val-corrval),data.gripapt,'UniformOutput',false);
    case {'trl','trial'}
        % correct data with base value
        data.gripapt = cellfun(@(x,y) bsxfun(@minus,x,y-corrval),data.gripapt,val,'UniformOutput',false);
    otherwise
        error('Baseline correction flag [%s] is not supported',data.basecorr.corrflg);
end


%% function reject_outliers
%----------------------------------------
function data = reject_outliers(cfg,data)
ntseries = length(cfg.plot.tseries);

% substitute whole trials with nans
if ~isempty(data.param)
    data.param(data.idxr.reject,:) = NaN;
end
for ts = 1:ntseries
    tsname = cfg.plot.tseries{ts};
    tmpidx = find(data.idxr.reject);
    for ti = 1:length(tmpidx)
        data.tseries.(tsname){tmpidx(ti)} = nan(size(data.tseries.(tsname){tmpidx(ti)}));
    end
end

% substitute a single parameter value (pairwise) with a nan
if ~isempty(data.param)
pairwisesel = match_str(cfg.plot.param,km_labelselection(cfg.plot.trlrej.pairwise,cfg.plot.param));
for p = pairwisesel(:)'
    pname = cfg.plot.param{p};
    if isfield(cfg.plot.trlrej,pname)
        tmpidxr = km_combidxr(data.idxr,cfg.plot.trlrej.(pname));
        data.param(tmpidxr.reject,p) = NaN;
    elseif isfield(data.idxr,pname)
        data.param(data.idxr.(pname),p) = NaN;
    end
end
end
% substitute time-series individually (pairwise) with nans
for ts = 1:ntseries
    tsname = cfg.plot.tseries{ts};
    if isfield(cfg.tseries.(tsname),'trlrej')
        tsidxr = km_combidxr(data.idxr,cfg.tseries.(tsname).trlrej);
        tmpidx = find(tsidxr.reject);
        for ti = 1:length(tmpidx)
            data.tseries.(tsname){tmpidx(ti)} = nan(size(data.tseries.(tsname){tmpidx(ti)}));
        end
    end
end
%--------------------------------------------------------------------------


%% function exclude_cond
%----------------------------------------
function data = exclude_cond(cfg,data)
ntseries = length(cfg.plot.tseries);

if ~isempty(cfg.plot.factexcl) && ~isempty(data.idxe)
    % combine excluded factor x levels
    tmpe = false(size(data.idxe,1),1);
    for fe = 1:length(cfg.plot.factexcl)
        fename = cfg.plot.factexcl{fe};
        ie = ismember(data.idxe(:,fe),cfg.plot.lvlexcl.(fename));
        tmpe(ie) = true;
    end
    data.idxe = tmpe;
    
    % substitute excluded trials with nans
    if ~isempty(data.param)
        data.param(data.idxe,:) = NaN;
    end
    for ts = 1:ntseries
        tsname = cfg.plot.tseries{ts};
        tmpidx = find(data.idxe);
        for ti = 1:length(tmpidx)
            data.tseries.(tsname){tmpidx(ti)} = nan(size(data.tseries.(tsname){tmpidx(ti)}));
        end
    end
end
%--------------------------------------------------------------------------


%% function tseries_norm
%----------------------------------------
function tseries = tseries_norm(cfg,tseries)
ntseries = length(cfg.plot.tseries);

for ts = 1:ntseries
    tsname = cfg.plot.tseries{ts};
    for t = 1:length(tseries.time)
        % original data
        nd = ndims(tseries.(tsname){t});
        dims = size(tseries.(tsname){t});
        y = permute(tseries.(tsname){t},nd:-1:1);
        y = reshape(y,dims(nd:-1:1));
        if any(isnan(y(:)))
            yi = nan([101 dims(nd-1:-1:1)]);
        else
            % original time
            x = tseries.time{t}';
            % time to interpolate to
            xi = linspace(min(x),max(x),101)';
            yi = interp1(x,y,xi,'pchip','extrap');
        end
        % shape back to original size
        dims(end) = size(yi,1);
        yi = permute(yi,nd:-1:1);
        tseries.(tsname){t} = reshape(yi,dims);
    end
    % prepend the data with a new first dimension
    tseries.(tsname) = cellfun(@(x) shiftdim(x,-1),tseries.(tsname),'UniformOutput',false);
    % concatenate trials in the first dimension
    tseries.(tsname) = vertcat(tseries.(tsname){:});
    % update dimension order
    tseries.dimord.(tsname) = ['rept_' tseries.dimord.(tsname)];
end
if ~isfalse(cfg.plot.tseries)
    % remove time field
    tseries = rmfield(tseries,'time');
end
%--------------------------------------------------------------------------


%% function rename_sessfact
%----------------------------------------
function [CI,epe] = rename_sessfact(CI,epe)

% rename session factor for parameter confidence intervals
if ~isempty(CI)
    %sessfact = cfg.plot.CI.sessfact;
    sessfact = 'site';
    % loop over factors
    for f = 1:CI.nfact
        % replace fact
        CI.fact{f}(ismember(CI.fact{f},'sess')) = {sessfact};
        % replace fname
        CI.fname{f} = regexprep(CI.fname{f},'sess',sessfact);
    end
    % replace bfact
    CI.bfact(ismember(CI.bfact,'sess')) = {sessfact};
    % replace lvl
    if ~isfield(CI.lvl,sessfact) && isfield(CI.lvl,'sess')
        CI.lvl.(sessfact) = CI.lvl.sess;
        CI.lvl = rmfield(CI.lvl,'sess');
    end
    % replace field names
    flds = fieldnames(CI);
    idx_flds = cellfun(@(x) ~isempty(strfind(x,'sess')),flds);
    flds = flds(idx_flds);
    for f = 1:length(flds)
        newfname = regexprep(flds{f},'sess',sessfact);
        CI.(newfname) = CI.(flds{f});
        CI = rmfield(CI,flds{f});
    end
end

% rename session factor for end point error estimations
if ~isempty(epe) && isfield(epe,'fact')
    %sessfact = cfg.plot.epe.sessfact;
    sessfact = 'site';
    % loop over factors
    for f = 1:epe.nfact
        % replace fact
        epe.fact{f}(ismember(epe.fact{f},'sess')) = {sessfact};
        % replace fname
        epe.fname{f} = regexprep(epe.fname{f},'sess',sessfact);
    end
    % replace bfact
    epe.bfact(ismember(epe.bfact,'sess')) = {sessfact};
    % replace lvl
    if ~isfield(epe.lvl,sessfact) && isfield(epe.lvl,'sess')
        epe.lvl.(sessfact) = epe.lvl.sess;
        epe.lvl = rmfield(epe.lvl,'sess');
    end
    % replace field names
    flds = fieldnames(epe);
    idx_flds = cellfun(@(x) ~isempty(strfind(x,'sess')),flds);
    flds = flds(idx_flds);
    for f = 1:length(flds)
        newfname = regexprep(flds{f},'sess',sessfact);
        epe.(newfname) = epe.(flds{f});
        epe = rmfield(epe,flds{f});
    end
end
%--------------------------------------------------------------------------
    