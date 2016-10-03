function cfg = km_movparam(cfg,data)
%--------------------------------------------------------------------------
%
% Movement parameters
% -------------------
% All movement parameters are denoted with an abbreviation: 'rt' stands for
% 'reaction time' and 'mv' stands for 'mean velocity'. Some parameters can
% depend on a selected movement phase, a sensor, or an axis. Please see
% below for more information.
%
% TIME related
% rt        - reaction time
% mt        - movement duration
% rmt       - relative movement phase duration, e.g. mt_T/mt
%
% VELOCITY related
% v         - instantaneous velocity (at onset of the denoted phase)
% mv        - mean velocity
% tl        - trajectory length (integral of velocity)
% pv        - peak velocity (positive peak velocity 'ppv' is also accepted)
% tpv       - time to peak velocity (from movement onset)
% rtpv      - relative time to peak velocity (to movement duration)
% apv, tapv, rtapv 	- peak of absolute velocity
% npv, tnpv, rtnpv 	- negative peak of velocity
%
% ACCELERATION related
% pa, tpa, rtpa	- peak acceleration (similar to 'pv')
% pd, tpa, rtpd - peak deceleration (similar to 'pa')
% nvc       - number of velocity changes
%
% GRIP related
% ga        - instantaneous grip aperture (at onset of the denoted phase)
% pga, tpga,rtpga       - peak grip aperture
% gv        - instantaneous grip velocity (at onset of the denoted phase)
% mgv       - mean grip velocity
% zgv       - integral of grip velocity (similar to 'tl')
% pgv, tpv, rtpgv       - peak grip velocity ('ppgv' is also accepted; similar to 'pv')
% npgv, tnpgv, rtnpgv   - negative peak grip velocity (as 'npv')
% go        - instantaneous grip orientation (at onset of the denote phase)
%
% POSITION related
% pos       - instantaneous position (at onset of the denoted phase)
% maxpos    - maximal position (maximal positive position 'pmaxpos' is also accepted)
% amaxpos   - maximal absolute position
% nmaxpos   - maximal negative position
% zpos      - integral of position (similar to 'tl')
% pzpos     - integral of positive position (similar to 'tl')
% nzpos     - integral of negative position (similar to 'tl')
%
%
% Movement phases
% ---------------
% The movement can be subdivided in several phases (e.g. transport and
% approach). When applicable, parameters can be calculated separately for
% these phases. A parameter specific for a single movement phase is
% appended with an underscore '_' and a capital (!) letter indicating the
% phase: 'mv_T' stands for the 'mean velocity of the transport phase'.
%
% common phases
% T - transport - from start to a defined end of the transport
% A - approach  - from end of transport to end of movement
% B - beginning - an epoch time-locked to the beginning of the movement
% E - end       - an epoch time-locked to the end of the movement
% W - window    - a fully specifyable window of interest
%
% special phases
% Tt - thumb transport - transport phase defined for the thumb specifically
% Ti - index transport - transport phase defined for the index specifically
%
% reference phases
% R  - reference 	- a fully specifyable reference window
% Rb - beginning	- an epoch time-locked to the beginning of the movement
% Re - end          - an epoch time-locked to the end of the movement
%
%
% Differential and relative parameters
% ---------------
% For almost all parameters a differential value can be calculated, either
% with respect to a reference window, or maybe another phase. The names of
% differential parameters are prepadded by a 'd'. The reference
% window/phase can be indicated by postpadding the name with an underscore
% '_' and the phase denotation. Example 1: 'mt_T' is the duration of the
% transport movement phase. Example 2: 'dmt_T' is the difference between
% the total movement duration and the transport phase duration. Example 3:
% 'dmt_T_A' is the difference between the transport and approach durations.
% There is one interesting addition: if the parameter name is prepadded by
% 'ad', than the absolute difference is taken, instead of the signed
% difference. This is not frequently used, but might be relevant to
% calculate the absolute differenctial grip orientation ('adgo').
%
%
% Sensors
% -------
% ATTENTION: Sensor selection is not yet implemented
% Some parameters can be calculated separately for individual sensors, if
% so, the sensor indicator is appended after the parameter abbreviation
% following a underscore '_' and a lower case (!) letter indicating the
% sensor: 'mv_T_t' stands for the 'mean velocity of the thumb during the
% transport phase'.
%
% common sensors
% t - thumb
% i - index
% j - first metacarpophalangeal joint
% w - wrist
% o - object
%
%
% Axes
% ----
% ATTENTION: Axis selection is not yet implemented
% Some parameters can depend on the axis dimension, 'position' for example.
% Suggestion: append '_ax' where ax denotes the axis in lower case.
%
%
% Movement number
% ---------------
% ATTENTION: Movement number selection is not yet implemented
% A trial can contain more than one movement (e.g. towards target and
% back). Parameters can be calculated per movement, not per trial.
% Suggestion: append '_#' where # denotes the movement number.
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2011-11-09
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end

% FIXME: keep for legacy reasons only: app > apt
if isfield(data,'grip')
    if isfield(data.grip,'appaxis')
        data.grip.aptaxis = data.grip.appaxis;
        data.grip = rmfield(data.grip,'appaxis');
    end
    if isfield(data.grip,'app')
        data.grip.apt = data.grip.app;
        data.grip = rmfield(data.grip,'app');
    end
    if isfield(data.grip,'appvel')
        data.grip.aptvel = data.grip.appvel;
        data.grip = rmfield(data.grip,'appvel');
    end
    if isfield(data.grip,'appacc')
        data.grip.aptacc = data.grip.appacc;
        data.grip = rmfield(data.grip,'appacc');
    end
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

% get movement
if isfield(cfg,'movement') && ~isstruct(cfg.movement) && ~isempty(cfg.movement)
    movement = cfg.movement;
elseif isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
    movement = cfg.movdef.movement;
elseif isfield(data,'movement')
    movement = data.movement;
elseif isfield(cfg,'movdef')
    movement = km_movdef(cfg,data);
else
    error('no movement provided');
end
% make sure that all trials have the same number of movements
movement = km_setequalnmov(movement);
% get the data on which the parameters are based
mpcfg = cfg.movparam;

% if movement parameters have been previously calculated, remove them
% completely if this is requested
if isfield(mpcfg,'pname') && ~isfalse(mpcfg.clean)
    idx_pname = ismember(cfg.vars,mpcfg.pname);
    vars = vars(~idx_pname);
    trl = trl(:,~idx_pname);
end

% execute study specific function if requested
if ~isfalse(cfg.movparam.expfun)
    if strcmpi(cfg.movparam.expfun,'exp')
        cfg = setcfg_exp(cfg);
        cfg.movparam.expfun = ['km_movparamfun_' cfg.exp];
    end
    [cfg,flg_return] = feval(cfg.movparam.expfun,cfg,data,'begin');
    if istrue(flg_return), return; end
end

%% movement (phase) definitions
%----------------------------------------
% movement number
if isfield(mpcfg,'movnr')
    movnr = mpcfg.movnr;
else
    movnr = km_getmovidx(mpcfg,'time');
end
nmov = length(movnr);
% keep only selected movement numbers
mov.total = cellfun(@(x) [x(movnr,1) x(movnr,2)],movement,'UniformOutput',false);
% movement begin and end
%[mov.beg.total,mov.end.total] = mov2begend(mov.total);

% identity absolute parameters
prmabs = regexp(mpcfg.param,'(?<=^d?)a[a-zA-Z]');
prmabs(cellfun(@isempty,prmabs)) = {0}; prmabs = [prmabs{:}];

% identify differential parameters
prmdiff = regexp(mpcfg.param,'(?<=^a?)d[a-zA-Z]');
prmdiff(cellfun(@isempty,prmdiff)) = {0}; prmdiff = [prmdiff{:}];

% check if absolute and differential tokens match
if any(prmabs>3) || any(prmdiff>3) || ...
   any(prmabs==0&prmdiff==2) || any(prmabs==2&prmdiff==0) || ...
   any(prmabs==1&prmdiff==1) || any(prmabs==2&prmdiff==2)
    error('KN:MovParam:Prefix','The parameter name prefix (with ''a'' and ''d'') could not be correctly identified.');
end

% isolate parameter names
[prmtok,prmlst] = list_tokens(mpcfg.param,'^[ad]{0,2}([a-zA-Z]+)');

% isolate phase indicators
[phstok,phslst,phsreflst] = list_tokens(mpcfg.param,'_([A-Z][a-zA-Z]*)');

% sort phases of interest and reference
tmp = (prmdiff>0) & cellfun(@isempty,phsreflst);
phsreflst(tmp) = phslst(tmp);
phslst(tmp) = {''};
phslst(cellfun(@isempty,phslst)) = {'total'};

% check if all differential parameters also have a reference phase
if any(cellfun(@isempty,phsreflst(prmdiff>0)))
    error('KN:MovParam:Prefix','For at least one differential parameter no reference phase (e.g. ''_Re'' is specified.');
end
if any(~cellfun(@isempty,phsreflst(prmdiff==0)))
    error('KN:MovParam:Prefix','For at least one parameter a reference phase (e.g. ''_Re'' is specified, but the parameter name is not prepadded with a ''d''. The reference phase will be ignored.');
end

% isolate sensor indicators
% FIXME: sensor indicators are not implemented yet
[snstok,snslst,snsreflst] = list_tokens(mpcfg.param,'_([a-z]+)');

% isolate axis indicators
% FIXME: axis indicators are not implemented yet

% isolate movement number indicators
% FIXME: movement indicators are not implemented yet
[movtok,movlst,movreflst] = list_tokens(mpcfg.param,'_([0-9]+)');

% transport and approach phases of the movement
if any(ismember(phstok,{'T','A'})) || any(ismember(prmtok,{'pa','tpa','rtpa','pd','tpd','rtpd'}))
    %if any(nmov>1), error('currently only one movement per trials is supported');   end
    tmpcfg = mpcfg;
    tmpcfg.vel = mpcfg.vel.T;
    tmpcfg.vel.T = mpcfg.vel.T;
    [mov.T,mov.A] = get_mov_TA(tmpcfg,data,mov.total,trl,vars);
end
% thumb transport movement
if any(ismember(phstok,{'Tt','At'}))
    tmpcfg = mpcfg;
    tmpcfg.vel = mpcfg.vel.Tt;
    tmpcfg.vel.T = mpcfg.vel.Tt;
    [mov.Tt,mov.At] = get_mov_TA(tmpcfg,data,mov.total,trl,vars);
end
% index transport movement
if any(ismember(phstok,{'Ti','Ai'}))
    tmpcfg = mpcfg;
    tmpcfg.vel = mpcfg.vel.Ti;
    tmpcfg.vel.T = mpcfg.vel.Ti;
    [mov.Ti,mov.Ai] = get_mov_TA(tmpcfg,data,mov.total,trl,vars);
end
% begin movement
if any(ismember(phstok,'B'))
    mov.B = offset2mov(mpcfg.offset.B,mov.total,'beg');
end
% end movement
if any(ismember(phstok,'E'))
    mov.E = offset2mov(mpcfg.offset.E,mov.total,'end');
end
% arbitrary window of interest
if any(ismember(phstok,'W'))
    mov.W = offset2mov(mpcfg.offset.W.offset,mov.total,mpcfg.offset.W.ref);
end
% arbitrary reference
if any(ismember(phstok,'R'))
    mov.R = offset2mov(mpcfg.offset.R.offset,mov.total,mpcfg.offset.R.ref);
end
% reference begin movement
if any(ismember(phstok,'Rb'))
    mov.Rb = offset2mov(mpcfg.offset.Rb,mov.total,'beg');
end
% reference end movement
if any(ismember(phstok,'Re'))
    mov.Re = offset2mov(mpcfg.offset.Re,mov.total,'end');
end


%% time
%----------------------------------------
% find the earliest and latest time point in any mov
[mintim,maxtim] = findminmaxtim(mov);
mov.minmax = cellfun(@(x,y) [x y],num2cell(mintim),num2cell(maxtim),'UniformOutput',false);

% obtain time window from data
timminmax = cellfun(@(x,y) x(time2logic(y,x)),data.time,mov.minmax,'UniformOutput',false);


%% velocity & acceleration
%----------------------------------------
% create empty data structures
dat = struct;
dat.vel = [];
dat.dist = [];
dat.acc = [];
dat.accvel = [];

% velocity
if any(~cellfun(@isempty,regexp(prmtok,'|(^v)|(mv)|(pv)$','once')))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.vel = mpcfg.vel(1);
    [dat.vel, dims.vel.ax, dims.vel.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'vel',false);
    dat.vel = cellfun(@abs,dat.vel,'UniformOutput',false);
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.vel),2);
    tmpcfg.vel = mpcfg.vel(i);
    [dat.ref.vel, dims.ref.vel.ax, dims.ref.vel.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'vel',false);
    dat.ref.vel = cellfun(@abs,dat.ref.vel,'UniformOutput',false);
end
% velocity to derive trajectory length
if any(ismember(prmtok,'tl'))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.vel = mpcfg.dist(1);
    [dat.dist, dims.dist.ax, dims.dist.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'vel',false);
    dat.dist = cellfun(@abs,dat.dist,'UniformOutput',false);
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.dist),2);
    tmpcfg.vel = mpcfg.dist(i);
    [dat.ref.dist, dims.ref.dist.ax, dims.ref.dist.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'vel',false);
    dat.ref.dist = cellfun(@abs,dat.ref.dist,'UniformOutput',false);
end
% acceleration
if any(~cellfun(@isempty,regexp(prmtok,'(^a)|(nvc)|(pa)|(pd)$','once')))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.acc = mpcfg.acc(1);
    [dat.acc, dims.acc.ax, dims.acc.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'acc',false);
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.acc),2);
    tmpcfg.acc = mpcfg.acc(i);
    [dat.ref.acc, dims.ref.acc.ax, dims.ref.acc.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'acc',false);
end
% velocity to aid in peak ac-/de-celeration determination
if any(~cellfun(@isempty,regexp(prmtok,'(pa)|(pd)$','once')))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.vel = mpcfg.acc(1);
    [dat.accvel, dims.accvel.ax, dims.accvel.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'vel',false);
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.acc),2);
    tmpcfg.vel = mpcfg.acc(i);
    [dat.ref.accvel, dims.ref.accvel.ax, dims.ref.accvel.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'vel',false);
end


%% grip descriptives
%----------------------------------------
% create empty data structures
dat.apt = [];
dat.aptvel = [];
dat.ori = [];

% grip aperture
if any(~cellfun(@isempty,regexp(prmtok,'ga$','once')))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.gripapt = mpcfg.gripapt(1);
    [dat.apt, dims.apt.ax, dims.apt.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'gripapt',false);
    dat.apt = cellfun(@abs,dat.apt,'UniformOutput',false);
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.gripapt),2);
    tmpcfg.gripapt = mpcfg.gripapt(i);
    [dat.ref.apt, dims.ref.apt.ax, dims.ref.apt.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'gripapt',false);
    dat.ref.apt = cellfun(@abs,dat.ref.apt,'UniformOutput',false);
end
% aperture velocity
if any(~cellfun(@isempty,regexp(prmtok,'gv$','once')))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.gripaptvel = mpcfg.gripaptvel(1);
    [dat.aptvel, dims.aptvel.ax, dims.aptvel.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'gripaptvel',false);
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.gripaptvel),2);
    tmpcfg.gripaptvel = mpcfg.gripaptvel(i);
    [dat.ref.aptvel, dims.ref.aptvel.ax, dims.ref.aptvel.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'gripaptvel',false);
end
% grip orientation
if any(~cellfun(@isempty,regexp(prmtok,'go$','once')))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.gripori = mpcfg.gripori(1);
    [dat.ori, dims.ori.ax, dims.ori.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'gripori',false);
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.gripori),2);
    tmpcfg.gripori = mpcfg.gripori(i);
    [dat.ref.ori, dims.ref.ori.ax, dims.ref.ori.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'gripori',false);
end


%% position
%----------------------------------------
% create empty data structures
dat.pos = [];

% position
if any(~cellfun(@isempty,regexp(prmtok,'pos$','once')))
    % timeseries of interest
    tmpcfg = mpcfg;
    tmpcfg.pos = mpcfg.pos(1);
    [dat.pos, dims.pos.ax, dims.pos.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'pos',false);
    dat.pos(cellfun(@isempty,dat.pos)) = {nan(length(dims.pos.sns),length(dims.pos.ax))};
    % timeseries for reference
    tmpcfg = mpcfg; i = min(length(mpcfg.pos),2);
    tmpcfg.pos = mpcfg.pos(i);
    [dat.ref.pos, dims.ref.pos.ax, dims.ref.pos.sns] = km_getmovdat(tmpcfg,data,mov.minmax,'pos',false);
    dat.ref.pos(cellfun(@isempty,dat.ref.pos)) = {nan(length(dims.ref.pos.sns),length(dims.ref.pos.ax))};
end


%% calculate movement parameters
%----------------------------------------
% get number of sensors and axes of reference data
if isfield(dims,'ref')
    flds = fieldnames(dims.ref);
    for f = 1:length(flds)
        ns.ref.(flds{f}) = length(dims.ref.(flds{f}).sns);
        na.ref.(flds{f}) = length(dims.ref.(flds{f}).ax);
    end
end
% get number of sensors and axes
flds = fieldnames(dims); flds = flds(~ismember(flds,'ref'));
for f = 1:length(flds)
    ns.(flds{f}) = length(dims.(flds{f}).sns);
    na.(flds{f}) = length(dims.(flds{f}).ax);
    % check if reference matches the original
    if ns.ref.(flds{f}) ~= ns.(flds{f})
        error('KM:MovParam:RefDims','For the parameter ''%s'' the number of reference sensors [%i] does not match the number of original sensors [%i].',flds{f},ns.ref.(flds{f}),ns.(flds{f}));
    end
    if na.ref.(flds{f}) ~= na.(flds{f})
        error('KM:MovParam:RefDims','For the parameter ''%s'' the number of reference axes [%i] does not match the number of original axes [%i].',flds{f},na.ref.(flds{f}),na.(flds{f}));
    end
end
% add 'none' field for compatibility
ns.none = 1;        na.none = 1;
ns.ref.none = 1;    na.ref.none = 1;

% loop over parameters
pnamelist = {};
for p = 1:length(prmlst)
    
    % get parameter name
    prm = prmlst{p};
    
    % get phase name
    phs = phslst{p};
    tim = mov.(phs);
    [~,timidx] = cellfun(@(x,y) time2logic(x,y),tim,timminmax,'UniformOutput',false);
    
    % get sensor name
    sns = snslst{p};
    % FIXME: do something with the sensors name
    
    % set variable input arguments for get_par
    switch prm
        case 'rmt'
            varargin = {mov.total};
        case {'tl','zgv','zpos','pzpos','nzpos'}
            varargin = {data.fsample};
        otherwise
            varargin = {[]};
    end
    
    % get parameter and timeseries name ofr each movement separately
    par = cell(1,nmov);
    for m = 1:nmov
        tim_mov = cellfun(@(x) x(m,:),tim,'UniformOutput',false);
        timidx_mov = cellfun(@(x) x(m,:),timidx,'UniformOutput',false);
        [par{m},tsname] = get_par(prm,tim_mov,timidx_mov,dat,ns,na,timminmax,varargin{:});
    end
    % concatenate movements in fourth dimension [ntrl nsns nax nmov];
    par = cat(4,par{:});
    
    % get reference time window if this is a difference parameter
    if prmdiff(p)>0
        phsref = phsreflst{p};
        timref = mov.(phsref);
        timrefidx = cellfun(@(x,y) time2logic(y,x),timminmax,timref,'UniformOutput',false);
        
        % if a dedicated reference time window is selected for an
        % instantaneous parameter, the mean is used for referencing
        prmref = prm; datref = dat;
        if ismember(phsref,{'R','Rb','Re'})
            switch prm
                case {'v','mv','pv','ppv','apv','npv'}
                    prmref = 'mv';  datref.vel = dat.ref.vel;
                case {'a','ma','pa','pd'}
                    prmref = 'ma';  datref.vel = dat.ref.acc;
                case {'ga','mga','pga'}
                    prmref = 'mga'; datref.apt = dat.ref.apt;
                case {'gv','mgv','pgv','npgv'}
                    prmref = 'mgv'; datref.aptvel = dat.ref.aptvel;
                case {'go'}
                    prmref = 'mgo';	datref.ori = dat.ref.ori;
                case {'pos','maxpos','amaxpos','nmaxpos'}
                    prmref = 'mpos'; datref.pos = dat.ref.pos;
                case{'tl','zgv','zpos','pzpos','nzpos'}
                    warning('KM:MovParam:Reference','Proper referencing is currently not implemented for parameter [%s]',prm);
            end
        end
        
        % get reference parameter
        ref = get_par(prmref,timref,timrefidx,datref,ns,na,timminmax,varargin{:});
        
        % make parameter and reference absolute before subtraction if requested
        if prmabs(p) == 2
            ref = abs(ref);
        	par = abs(par);
        end
        
        % subtract reference from parameter
        par = par - ref;
    end
    
    % make parameter absolute if requested
    if prmabs(p) == 1
        par = abs(par);
    end
    
    % construct 'absolute' and 'difference' prefix
    if prmabs(p)==1 && prmdiff(p)==0
        prefix = 'a';
    elseif prmabs(p)==0 && prmdiff(p)==1
        prefix = 'd';
    elseif prmabs(p)==1 && prmdiff(p)==2
        prefix = 'ad';
    elseif prmabs(p)==2 && prmdiff(p)==1
        prefix = 'da';
    else
        prefix = '';
    end
    
    % construct parameter name for storage
    if prmdiff(p)>0
        if strcmp(phs,'total')
            pname = sprintf('%s%s_%s',prefix,prm,phsref);
        else
            pname = sprintf('%s%s_%s_%s',prefix,prm,phs,phsref);
        end
    else
        if strcmp(phs,'total')
            pname = sprintf('%s%s',prefix,prm);
        else
            pname = sprintf('%s%s_%s',prefix,prm,phs);
        end
    end
    
    % convert pname to cell structure
    pname = repmat({pname},[ns.(tsname) na.(tsname) nmov]);
    
    % add sensor indicator
    if ns.(tsname) > 1
        tmp = repmat(dims.(tsname).sns',[1 na.(tsname) nmov]);
        pname = cellfun(@(x,y) [x '_' y],pname,tmp,'UniformOutput',false);
    end
    
    % add axis indicator
    if na.(tsname) > 1
        tmp = repmat(dims.(tsname).ax,[ns.(tsname) 1 nmov]);
        pname = cellfun(@(x,y) [x '_' y],pname,tmp,'UniformOutput',false);
    end
    
    % add movement number
    if nmov > 1
        tmp = cell(1,1,nmov);
        tmp(1,1,:) = cellfun(@(x) sprintf('%i',x),num2cell(movnr),'UniformOutput',false);
        tmp = repmat(tmp,[ns.(tsname) na.(tsname) 1]);
        pname = cellfun(@(x,y) [x '_' y],pname,tmp,'UniformOutput',false);
    end
    
    % check pname against existing vars
    pname = reshape(pname,[1 ns.(tsname)*na.(tsname)*nmov]);
    [~, idx_vars] = ismember(pname,vars);
    nvars = length(vars);
    idx_vars(idx_vars==0) = (nvars+1):(nvars+sum(idx_vars==0));
    
    % store variable name (overwrite duplicates)
    vars(idx_vars) = pname;
    pnamelist = [pnamelist pname];
    
    % store parameter values (overwrite duplicates)
    par = reshape(par,[size(par,1) ns.(tsname)*na.(tsname)*nmov]);
    trl(:,idx_vars) = par;
    
end

% check
npar = size(trl,2);
nvars = length(vars);
if size(trl,2) ~= length(vars)
    error('The number of parameters in the vars [%d] does not match the values [%d]',nvars,npar);
end

% store trl and vars in configuration
cfg.trl = trl;
cfg.vars = vars;
cfg.movparam.pname = pnamelist;

% store new movement
if isfield(cfg,'movement')
    cfg.movement = movement;
    if exist('tmovement','var'),    cfg.tmovement = tmovement;  end
end
if isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
    for i = 1:length(cfg.movdef)
        cfg.movdef(i).movement = movement;
        if exist('tmovement','var'),    cfg.movdef(i).tmovement = tmovement;  end
    end
end
if isfield(data,'movement')
    data.movement = movement;
    if exist('tmovement','var'),    data.tmovement = tmovement;  end
end

% execute study specific function if requested
if ~isfalse(cfg.movparam.expfun)
    cfg = feval(cfg.movparam.expfun,cfg,data,'end');
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


% function list_tokens
%----------------------------------------
function [tok,lst_first,lst_second] = list_tokens(str,tokstr)
% capture tokens
tok = regexp(str,tokstr,'tokens');
% add empty strings where no tokens were found
lst_tok = tok; lst_tok(cellfun(@isempty,tok)) = {{{''}}};
% take first token for each entry
lst_first = cellfun(@(x) x{1}{1},lst_tok,'UniformOutput',false);
% find entries where two tokens were captured
idx_second = cellfun(@(x) size(x,2)>1,lst_tok);
% take second token for each entry (if available)
lst_tmp = cellfun(@(x) x{2}{1},lst_tok(idx_second),'UniformOutput',false);
% replace entries without a second token by empty strings
lst_second = repmat({''},size(lst_tok)); lst_second(idx_second) = lst_tmp(:);
% concatenate tokens (and leave out empty ones)
tok = [tok{:}];
% find unique tokens
tok = unique([tok{:}]);
if isempty(tok), tok = cell(0); end;


% function get_mov_TA
%----------------------------------------
function [mov_T,mov_A] = get_mov_TA(mpcfg,data,mov,trl,vars)
% movement number
% FIXME: mov_T should be created for all selected movnr.
movnr = km_getmovidx(mpcfg,'time');
% transport velocity
vel_T = km_getmovdat(mpcfg,data,mov,'vel');
vel_T = cellfun(@abs,vel_T,'UniformOutput',false);
% apply transport velocity criterium
if ischar(mpcfg.vel.T.crit)
    %mov_T = eval(mpcfg.vel.T.crit);
    mov_T = cellfun(@(dat) eval(mpcfg.vel.T.crit),vel_T,'UniformOutput',false);
elseif length(mpcfg.vel.T.crit) == 2
    mov_T = cellfun(@(x) ...
        x > min(mpcfg.vel.T.crit) & x < max(mpcfg.vel.T.crit),...
        vel_T,'UniformOutput',false);
else
    mov_T = cellfun(@(x) x>mpcfg.vel.T.crit,vel_T,'UniformOutput',false);
end
% slide borders to local minimum or criterium
if ~isfalse(mpcfg.vel.T.slide)
    %mov_T = km_slidesel(vel_T,tmov,mpcfg.vel.T.slide,mpcfg.vel.T.consecflg);
    mov_T = cellfun(@(x,y) km_slidesel(x,y,mpcfg.vel.T.slide,mpcfg.vel.T.consecflg),vel_T,mov_T,'UniformOutput',false);
end
% go from the sample to the index relative to movement onset
idx_T = cellfun(@logic2idx,mov_T,'UniformOutput',false);
idx_T(cellfun(@isempty,idx_T)) = {NaN};
idx_T = cellfun(@(x) max(x(:)),idx_T);
% if mov_T includes the last sample, go one sample back
nsamp = cellfun(@length,mov_T);
idx_T(idx_T==nsamp & nsamp>1) = idx_T(idx_T==nsamp & nsamp>1) - 1;
% let mov_A start one sample after mov_T ends
idx_A = idx_T+1;
idx_A(idx_A>nsamp) = nsamp(idx_A>nsamp);
% go from the index to the time domain relative to trial start
mov_T = km_gettim(data.time,trl,vars,idx_T)';
mov_T = mov_T + cellfun(@(x) x(movnr,1),mov);
mov_A = km_gettim(data.time,trl,vars,idx_A)';
mov_A = mov_A + cellfun(@(x) x(movnr,1),mov);
% add movement start and transform to cell array for consistency in handling
%tmov = num2cell([mov_T;mov_T]',2)';
mov_T = cellfun(@(x,y) [x(movnr,1) y],mov,num2cell(mov_T),'UniformOutput',false);
mov_A = cellfun(@(x,y) [y x(movnr,2)],mov,num2cell(mov_A),'UniformOutput',false);
% make sure that all trials have the same number of movements
mov_T = km_setequalnmov(mov_T);
mov_A = km_setequalnmov(mov_A);

% function offset2mov
%----------------------------------------
function mov = offset2mov(offset,mov,ref)
if nargin<3,	ref = 'on';	end
offset = repmat({offset},size(mov));
switch lower(ref)
	case {'onset','on','begin','beg','b'}
        mov = cellfun(@(x,y) [x(:,1)+y(1) x(:,1)+y(end)],mov,offset,'UniformOutput',false);
	case {'offset','off','end','e'}
        mov = cellfun(@(x,y) [x(:,2)+y(1) x(:,2)+y(end)],mov,offset,'UniformOutput',false);
    otherwise
        % relative to t = 0
        error('not implemented yet');
        mov = cellfun(@(y) [y(1) y(end)],offset,'UniformOutput',false);
end

% function findminmaxtim
%----------------------------------------
function [mintim,maxtim] = findminmaxtim(mov)
flds = fieldnames(mov);
for f = 1:length(flds)
    if isstruct(mov.(flds{f}))
        [mintim_tmp,maxtim_tmp] = findminmaxtim(mov.(flds{f}));
    else
        mintim_tmp = cellfun(@(x) min(min(x)),mov.(flds{f}));
        maxtim_tmp = cellfun(@(x) max(max(x)),mov.(flds{f}));
    end
    if f == 1
        mintim = mintim_tmp;
        maxtim = maxtim_tmp;
    else
        mintim = min(mintim,mintim_tmp);
        maxtim = max(maxtim,maxtim_tmp);
    end
end


% function get_par
%----------------------------------------
function [par,ts,varargout] = get_par(prm,tim,idx,dat,ns,na,timminmax,varargin)

% switch over parameters
switch prm
    
    %% INDEPENDENT operations
    %----------------------------------------
    case 'firstval'
        ts = 'tmp'; dims = [ns.(ts) na.(ts) length(idx)];
        % get the first sample of the window
        tmpidx = cellfun(@logic2idx,idx,'UniformOutput',false);
        % replace empty indices by NaN
        tmpidx(cellfun(@isempty,tmpidx)) = {nan(1,2)};
        % get first index
        tmpidx = cell2mat(tmpidx'); tmpidx = tmpidx(:,1);
        idx_nan = isnan(tmpidx); tmpidx(idx_nan) = 1;
        % get the value of that sample
        par = cellfun(@(x,i) x(:,:,i(1)),dat.(ts),num2cell(tmpidx'),'UniformOutput',false);
        par = permute(reshape(cell2mat(par),dims),[3 1 2]);
        % if no proper index could be found, return NaN
        par(idx_nan,:) = NaN;
        
    case 'mean'
        ts = 'tmp'; dims = [ns.(ts) na.(ts) length(idx)];
        % calculate the mean value
        par = cellfun(@(x,i) nanmean(x(:,:,i),3),dat.(ts),idx,'UniformOutput',false);
        par = permute(reshape(cell2mat(par),dims),[3 1 2]);
        % if no proper index could be found, return NaN
        par(cellfun(@isempty,idx),:) = NaN;
        
    case 'integral'
        ts = 'tmp'; dims = [ns.(ts) na.(ts) length(idx)];
        % sum the values over all samples
        par = cellfun(@(x,i) nansum(x(:,:,i),3),dat.(ts),idx,'UniformOutput',false);
        par = permute(reshape(cell2mat(par),dims),[3 1 2]);
        % divide by the sampling frequency
        par = bsxfun(@rdivide,par,varargin{1}');
        % if no proper index could be found, return NaN
        par(cellfun(@isempty,idx),:) = NaN;
        
    case 'peak'
        ts = 'tmp'; dims = [ns.(ts) na.(ts) length(idx)];
        % find the maximum value (and its index)
        [par, I] = cellfun(@(x) max(x,[],3),idx2dat(dat.(ts),idx),'UniformOutput',false);
        par = permute(reshape(cell2mat(par),dims),[3 1 2]);
        % if no proper index could be found, return NaN
        par(cellfun(@isempty,idx),:) = NaN;
        % return index to output if requested
        if nargout > 2
            I = permute(reshape(cell2mat(I),dims),[3 1 2]);
            % if no maximum is present, return also NaN for the time
            I(isnan(par)) = NaN;
            varargout = {I};
        end
        
    case 'timetopeak'
        tmpprm = get_tmpprm(prm,varargin{:});
        % get the index of the peak value
        [~,ts,I] = get_par(tmpprm,tim,idx,dat,ns,na,timminmax,varargin{:});
        % convert the index of the peak to a time point
        dims = size(I);
        if length(dims) > 2
            tpv = cell2mat(cellfun(@(x,y) idx2time(x,y),mat2cell(I,ones(dims(1),1),dims(2),dims(3)),timminmax','UniformOutput',false));
        else
            tpv = cell2mat(cellfun(@(x,y) idx2time(x,y),mat2cell(I,ones(dims(1),1),dims(2)),timminmax','UniformOutput',false));
        end
        % get the onset of the movement
        rt = get_par('rt',tim,idx,dat,ns,na,timminmax,varargin{:});
        par = bsxfun(@minus,tpv,rt);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'reltimetopeak'
        tmpprm = get_tmpprm(prm,varargin{:});
        % get the time of the peak value
        [tpv,ts] = get_par(tmpprm,tim,idx,dat,ns,na,timminmax,varargin{:});
        % get the movement duration
        mt = get_par('mt',tim,idx,dat,ns,na,timminmax,varargin{:});
        % calculate the relative time to peak velocity
        par = bsxfun(@rdivide,tpv,mt);
        
    case 'abspeak'
        ts = 'vel'; ns.tmp = ns.(ts); na.tmp = na.(ts);
        % make the data absolute
        dat.tmp = cellfun(@abs,dat.(ts),'UniformOutput',false);
        % get the absolute peak velocity (using 'peak')
        [par,~,I] = get_par('peak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'negpeak'
        ts = 'vel'; ns.tmp = ns.(ts); na.tmp = na.(ts);
        % flip the data around zero
        dat.tmp = cellfun(@(x) -x,dat.(ts),'UniformOutput',false);
        % get the negative peak velocity (using 'peak')
        [par,~,I] = get_par('peak',tim,idx,dat,ns,na,timminmax,varargin{:});
        par = -par;
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    
    %% TIME related
    %----------------------------------------
    case 'rt'       % reaction time (s)
        ts = 'none';
        % get the reaction time
        par = cellfun(@(x) x(:,1),tim,'UniformOutput',false);
        par = cell2mat(par)';
        % if no proper time could be found, return NaN
        par(cellfun(@isempty,tim),:) = NaN;
        par(cell2mat(cellfun(@(x) any(isnan(x),2)',tim,'UniformOutput',false)')) = NaN;
        
    case 'mt'       % movement time (s)
        ts = 'none';
        % get the movement duration
        par = cellfun(@(x) x(:,2)-x(:,1),tim,'UniformOutput',false);
        par = cell2mat(par)';
        % if no proper time could be found, return NaN
        par(cellfun(@isempty,tim),:) = NaN;
        par(cell2mat(cellfun(@(x) any(isnan(x),2)',tim,'UniformOutput',false)')) = NaN;
        
    case 'rmt'    	% relative movement time (%mt)
        warning('please check code');
        ts = 'none';
        % get the duration of the whole movement
        mt = cellfun(@(x) x(:,2)-x(:,1),varargin{1},'UniformOutput',false);
        mt = cell2mat(mt)';
        % get the durationof the phase of interest
        mtx = cellfun(@(x) x(:,2)-x(:,1),tim,'UniformOutput',false);
        mtx = cell2mat(mtx)';
        % calculate relative movement time
        par = mtx./mt;
        % if no proper time could be found, return NaN
        par(cellfun(@isempty,tim),:) = NaN;
        par(cell2mat(cellfun(@(x) any(isnan(x),2)',tim,'UniformOutput',false)')) = NaN;
        
        
    %% VELOCITY related
    %----------------------------------------
    case 'v'        % instantaneous velocity (m/s)
        ts = 'vel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the instantaneous velocity (using 'firstval')
        par = get_par('firstval',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'mv'       % mean velocity (m/s)
        ts = 'vel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the mean velocity (using 'mean')
        par = get_par('mean',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case {'tl','zv'}        % trajectory length (integral of velocity) (m*s)
        ts = 'dist'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the integral of position (using 'integral')
        par = get_par('integral',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case {'pv','ppv'}       % peak velocity (m/s)
        ts = 'vel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the peak velocity (using 'peak')
        [par,~,I] = get_par('peak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case {'tpv','tppv'}     % time to peak velocity from movement onset (s)
        tmpprm = 'pv';
        % get the index of the peak value (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case {'rtpv','rtppv'}   % relative time to peak velocity (%mt)
        tmpprm = 'tpv';
        % get  the relative time to the peak value (using 'reltimetopeak')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        
    case 'apv'      % absolute peak velocity (m/s)
        ts = 'vel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the absolute peak velocity (using 'abspeak')
        [par,~,I] = get_par('abspeak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'tapv'     % time to absolute peak velocity from movement onset (s)
        tmpprm = 'apv';
        % get the index of the absolute peak velocity (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'rtapv'    % relative time to absolute peak velocity (%mt)
        tmpprm = 'tapv';
        % get  the relative time to the absolute peak velocity (using 'reltimetopeak')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        
    case 'npv'       % negative peak velocity (m/s)
        ts = 'vel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the negative peak velocity (using 'abspeak')
        [par,~,I] = get_par('negpeak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'tnpv'     % time to negative peak velocity from movement onset (s)
        tmpprm = 'npv';
        % get the index of the negative peak velocity (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'rtnpv'    % relative time to negative peak velocity (%mt)
        tmpprm = 'tnpv';
        % get  the relative time to the negative peak velocity (using 'reltimetopeak')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
    
        
    %% ACCELERATION related
    %----------------------------------------
    case 'a'        % instantaneous acceleration (m/s²)
        ts = 'acc'; dat.tmp = dat.(ts); ns.tmpl = ns.(ts); na.tmp = na.(ts);
        % get the instantaneous acceleration (using 'firstval')
        par = get_par('firstval',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'ma'       % mean acceleration (m/s²)
        ts = 'acc'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the mean acceleration (using 'mean')
        par = get_par('mean',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'pa'       % peak acceleration [<tpv] (m/s²)
        ts = 'accvel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the index of the peak velocity (using 'timetopeak')
        [~,~,Ipv] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax);
        Ipv = num2cell(Ipv); dims = size(Ipv);
        if length(dims)<3, dims(3) = 1; end
        % copy idx for all sensors and axes
        tmpidx = repmat(idx',[1 dims(2:3)]);
        % cut window of interest off after the time of peak velocity
        tmpidx = cellfun(@(i,j) i&idx2logic(1:j,length(i)),tmpidx,Ipv,'UniformOutput',false);
        % split the acceleration data for each sensor and axis (like idx)
        ts = 'acc'; dims = [ns.(ts) na.(ts) length(idx)];
        tmpdat = cellfun(@(x) mat2cell(x,ones(1,dims(1)),ones(1,dims(2)),size(x,3)),dat.(ts),'UniformOutput',false);
        tmpdat = permute(reshape([tmpdat{:}],dims),[3 1 2]);
        % now tmpidx and tmpdat have the same structure
        % get the peak acceleration
        [par, I] = cellfun(@(x) max(x,[],3),idx2dat(tmpdat,tmpidx));
        % if no proper index could be found, return NaN
        par(cellfun(@isempty,idx),:) = NaN;
        % return index to output if requested
        if nargout > 2
            % if no peak acceleration is found, return NaN for the index
            I(isnan(par)) = NaN;
            varargout = {I};
        end
        
    case 'tpa'      % time to peak acceleration [<tpv] (s)
        tmpprm = get_tmpprm(prm,varargin{:});
        % get the time to peak acceleration (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'rtpa'     % relative time to peak acceleration [<tpv] (%mt)
        tmpprm = get_tmpprm(prm,varargin{:});
        % get the relative time to peak acceleration (using 'reltimetopeak')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        
    case 'pd'       % peak deceleration [>tpv] (m/s²)
        ts = 'accvel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the index of the maximum velocity (using 'timetopeak')
        [~,~,Ipv] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax);
        Ipv = num2cell(Ipv); dims = size(Ipv);
        if length(dims)<3, dims(3) = 1; end
        % copy idx for all sensors and axes
        tmpidx = repmat(idx',[1 dims(2:3)]);
        % limit the window of interest to only after the time of peak velocity
        tmpidx = cellfun(@(i,j) i&idx2logic(j:length(i),length(i)),tmpidx,Ipv,'UniformOutput',false);
        % split the acceleration data for each sensor and axis (like idx)
        ts = 'acc'; dims = [ns.(ts) na.(ts) length(idx)];
        tmpdat = cellfun(@(x) mat2cell(x,ones(1,dims(1)),ones(1,dims(2)),size(x,3)),dat.(ts),'UniformOutput',false);
        tmpdat = permute(reshape([tmpdat{:}],dims),[3 1 2]);
        % now tmpidx and tmpdat have the same structure
        % get the peak deceleration
        [par, I] = cellfun(@(x) min(x,[],3),idx2dat(tmpdat,tmpidx));
        % if no proper index could be found, return NaN
        par(cellfun(@isempty,idx),:) = NaN;
        % return index to output if requested
        if nargout > 2
            % if no peak acceleration is found, return NaN for the index
            I(isnan(par)) = NaN;
            varargout = {I};
        end
        
    case 'tpd'      % time to peak deceleration [>tpv] (s)
        tmpprm = get_tmpprm(prm,varargin{:});
        % get the time to peak deceleration (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'rtpd'     % relative time to peak deceleration [>tpv] (%mt)
        tmpprm = get_tmpprm(prm,varargin{:});
        % get the relative time to peak deceleration (using 'reltimetopeak')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        
    case 'nvc'      % number of velocity changes (ordinal)
        ts = 'acc'; dims = [ns.(ts) na.(ts) length(idx)];
        % find the number of sign changes of acc
        par = cellfun(@(x,i) sum(abs(diff(sign(x(:,:,i)),1,3)),3),dat.(ts),idx,'UniformOutput',false);
        par = permute(reshape(cell2mat(par),dims),[3 1 2]);
        % if no proper index could be found, return NaN
        par(cellfun(@isempty,idx),:) = NaN;
        
        
    %% GRIP related
    %----------------------------------------
    case 'ga'       % instantaneous grip aperture (m)
        ts = 'apt'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the instantaneous grip aperture (using 'firstval')
        par = get_par('firstval',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'mga'      % mean grip aperture (m)
        ts = 'apt'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the mean grip aperture veloctiy (using 'mean')
        par = get_par('mean',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'pga'      % peak grip aperture (m)
        % get the peak grip aperture after the minimum (using 'pd')
        ts = 'apt'; ns.accvel = ns.(ts); na.accvel = na.(ts);
        dat.accvel = cellfun(@(x) -x,dat.(ts),'UniformOutput',false);
        ns.acc = ns.(ts); na.acc = na.(ts);
        dat.acc = cellfun(@(x) -x,dat.(ts),'UniformOutput',false);
        [par,~,I] = get_par('pd',tim,idx,dat,ns,na,timminmax,varargin{:});
        par = -par;
        
        % % get the peak of the grip aperture (using 'peak')
        % ts = 'apt'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % [par,~,I] = get_par('peak',tim,idx,dat,ns,na,timminmax,varargin{:});
        
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'tpga'     % time to maximum grip aperture (s)
        tmpprm = 'pga';
        % get the time to peak grip aperture (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'rtpga'	% relative time to maximum grip aperture (%mt)
        tmpprm = 'tpga';
        % get the relative time to peak grip aperture (using 'reltimetopeak')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        
    case 'gv'       % instantaneous grip aperture velocity (m/s)
        ts = 'aptvel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the instantaneous grip aperture veloctiy (using 'firstval')
        par = get_par('firstval',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'mgv'      % mean grip aperture velocity (m/s)
        ts = 'aptvel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the mean grip aperture veloctiy (using 'mean')
        par = get_par('mean',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'zgv'      % integral of grip aperture velocity (m*s)
        ts = 'aptvel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the integral of the grip aperture veloctiy (using 'integral')
        par = get_par('integral',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case {'pgv','ppgv'} % peak grip aperture velocity (m/s)
        ts = 'aptvel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the peak of the grip aperture veloctiy (using 'peak')
        [par,~,I] = get_par('peak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case {'tpgv','tppgv'}	% time to peak grip aperture velocity (s)
        tmpprm = 'pgv';
        % get the time to peak grip aperture veloctiy (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case {'rtpgv','rtppgv'}	% relatice time to peak grip aperture velocity (%mt)
        tmpprm = 'tpgv';
        % get the relative time to peak grip aperture veloctiy (using 'rtpv')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        
    case 'npgv'     % negative peak grip aperture velocity (m/s)
        ts = 'aptvel'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the negative peak of the grip aperture veloctiy (using 'negpeak')
        [par,~,I] = get_par('negpeak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'tnpgv'     % time to negative peak grip velocity from movement onset (s)
        tmpprm = 'npgv';
        % get the index of the maximally negative velocity (using 'timetopeak')
        [par,ts,I] = get_par('timetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'rtnpgv'    % relative time to negative peak grip velocity (%mt)
        tmpprm = 'tnpgv';
        % get  the relative time to the maximally negative velocity (using 'reltimetopeak')
        [par,ts] = get_par('reltimetopeak',tim,idx,dat,ns,na,timminmax,varargin{:},tmpprm);
        
    case 'go'       % instantaneous grip orientation (deg)
        ts = 'ori'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the instantaneous grip orientation (using 'firstval')
        par = get_par('firstval',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'mgo'      % mean grip orientation (deg)
        ts = 'ori'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the mean grip aperture veloctiy (using 'mean')
        par = get_par('mean',tim,idx,dat,ns,na,timminmax,varargin{:});
        
        
    %% POSITION related
    %----------------------------------------
    case 'pos'      % position (m)
        ts = 'pos'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the instantaneous position (using 'firstval')
        par = get_par('firstval',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'mpos'      % mean position (m)
        ts = 'pos'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the mean grip aperture veloctiy (using 'mean')
        par = get_par('mean',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case {'maxpos','pmaxpos','ppos','pppos'}	% maximally positive position (m)
        ts = 'pos'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the maximal position (using 'peak')
        [par,~,I] = get_par('peak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case {'amaxpos','appos'}	% maximal absolute position (m)
        ts = 'pos'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the maximal absolute position (using 'abspeak')
        [par,~,I] = get_par('abspeak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case {'nmaxpos','nppos'}	% maximally negative position (m)
        ts = 'pos'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the maximal position (using 'negpeak')
        [par,~,I] = get_par('negpeak',tim,idx,dat,ns,na,timminmax,varargin{:});
        % return index to output if requested
        if nargout > 2, varargout = {I}; end
        
    case 'zpos'     % integral of position, area under curve (m*s)
        ts = 'pos'; dat.tmp = dat.(ts); ns.tmp = ns.(ts); na.tmp = na.(ts);
        % get the integral of position (using 'integral')
        par = get_par('integral',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'pzpos'	% integral of positive position (m*s)
        ts = 'pos'; ns.tmp = ns.(ts); na.tmp = na.(ts);
        % use only positive positions
        dat.tmp = cellfun(@(x) x.*(x>0),dat.(ts),'UniformOutput',false);
        % get the integral of the remaining positions (using 'integral')
        par = get_par('integral',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    case 'nzpos'	% integral of negative position (m*s)
        ts = 'pos'; ns.tmp = ns.(ts); na.tmp = na.(ts);
        % use only negative positions
        dat.tmp = cellfun(@(x) x.*(x<0),dat.(ts),'UniformOutput',false);
        % get the integral of the remaining positions (using 'integral')
        par = get_par('integral',tim,idx,dat,ns,na,timminmax,varargin{:});
        
    otherwise
        error('parameter [%s] is not recognized',prm);
end


% function get_tmpprm
%----------------------------------------
function tmpprm = get_tmpprm(prm,varargin)
if length(varargin) > 1
    tmpprm = varargin{end};
    if strcmpi(tmpprm,prm)
        tmpprm = prm(2:end);
    end
else
    tmpprm = prm(2:end);
end
% allow independent operations
if strcmpi(tmpprm,'imetopeak')
    tmpprm = 'peak';
elseif strcmpi(tmpprm,'eltimetopeak')
    tmpprm = 'timetopeak';
end


% function idx2dat
%----------------------------------------
function tmpdat = idx2dat(dat,idx)
% assignment by '=' is not allowed in cellfun, this can be circumvented in
% two ways: by division by zero (requires complex reshaping etc), or by
% calling a dedicated subfunction (i2x)
switch 'i2x'
    case 'zerodivision'
        % replace non-indexed samples by NaN (create NaN by 0/0)
        dims = size(dat{1});
        tmpdat = cellfun(@(x,i) x.*(0./reshape(repmat(i,[prod(dims(1:2)) 1]),size(x))+1),dat,idx,'UniformOutput',false);
    case 'i2x'
        % use eval to circumvent cellfun limitations
        tmpdat = cellfun(@i2x,dat,idx,'UniformOutput',false);
end

% function i2x
%----------------------------------------
function x = i2x(x,i)
x(:,:,~i) = NaN;