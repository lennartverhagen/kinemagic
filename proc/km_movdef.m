function varargout = km_movdef(cfg,data)
%--------------------------------------------------------------------------
% Define movements
%
% Recently two new features were added to this function, (1) the
% possibility to define on- and off-sets separately and (2) the possibility
% to combine multiple source of information to define movements. This opens
% up many combinations which means that this function is not thorougly
% checked for all of them.
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

% FIXME: legacy reasons for app > apt
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

% determine if movements are defined as a whole or on the basis of separate
% on- and off-sets.
defonoff = length(cfg.movdef)>1;

% determine if Multiple Sources of Information approach is used for
% movement definition
msi = ~isfalse(cfg.movdef(1).msi);

% desired feedback
fdbk = cfg.movdef(1).feedback;

% loop over runs
movwin = cell(length(data.time),length(cfg.movdef));
for r = 1:length(data.time)
    
    % loop over movement definitions (in one go, or on- and off-sets
    % seperately)
    for i = 1:length(cfg.movdef)
        
        % loop over parameters
        counter = struct;
        movwin{r,i} = nan(length(data.time{r}),length(cfg.movdef(i).param));
        for p = 1:length(cfg.movdef(i).param)
            pname = cfg.movdef(i).param{p};
            
            % counter
            if ~isfield(counter,pname)
                counter.(pname) = 1;
            else
                counter.(pname) = counter.(pname) + 1;
            end
            
            % parameter configuration
            tcfg = cfg.movdef(i).(pname);
            if iscell(tcfg)
                tcfg = tcfg{counter.(pname)};
            end
            
            % switch between parameters
            switch lower(pname)
                
                case 'pos'
                    mov = movdef_pos(tcfg,data,r);
                    
                case {'vel','acc','gripvel'}
                    mov = movdef_deriv(tcfg,data,r,pname);
                    
                case {'gripapt','gripacc'}
                    mov = movdef_grip(tcfg,data,r,pname);
                    
                case 'logfile'
                    if ~isfield(cfg,'log')
                        error('to use the logfile to determine the movement, the cfg structure must contain a ''log'' field');
                    end
                    mov = movdef_log(tcfg,data,r,cfg.log);
                    
                case 'probwin'
                    mov = movdef_probwin(tcfg,data.time{r},movwin{r,i});
                    
                case 'compwin'
                    mov = movdef_compwin(tcfg,data,r);
                    
                otherwise
                    error('parameter ''%s'' not supported',pname);
                    
            end
            
            % store movement windows
            movwin{r,i}(:,p) = mov;
            
        end
        
    end
end

% FIXME: give feedback
if ~isfalse(fdbk)
    
    % loop over movement definitions (in one go, or on- and off-sets
    % seperately)
    for i = 1:length(cfg.movdef)
        
        switch fdbk
            
            case 'movplot'
                % visually inspect the movements
                % setup tcfg structure
                tcfg        = [];
                tcfg.save   = 'no';
                tcfg.exp    = 'movdef';
                tcfg.dataset= 'movdef';
                % fill in movplot parameters
                if isfield(cfg.movdef(i),'movplot')
                    tcfg.movplot = cfg.movdef(i).movplot;
                elseif isfield(cfg,'movplot')
                    tcfg.movplot = cfg.movplot;
                else
                    tcfg.movplot.param           = {'vel','pos'};
                    tcfg.movplot.pos.marker      = 'mti';
                    tcfg.movplot.pos.axis        = {'x'};
                    tcfg.movplot.vel.marker      = 'mti';
                    tcfg.movplot.vel.axis        = 'yz';
                    tcfg.movplot.xlim            = [-0.5 3];
                end
                % fill in movplot fields specific for movdef
                tcfg.movplot.continuous = true;
                tcfg.movplot.movdef     = cfg.movdef(i);
                % setup tdata structure
                tdata                   = data;
                tdata.movwin            = movwin(:,i)';
                % run movplot
                km_movplot(tcfg,tdata);
                
            otherwise
                
        end
    end
end

% loop over runs
movement = cell(1,length(data.time));
for r = 1:length(data.time)
    
    % loop over movement definitions (in one go, or on- and off-sets
    % seperately)
    for i = 1:length(cfg.movdef)
        
        % msi: combine movement definition parameters
        if msi
            movwin{r,i} = prod(movwin{r,i},2);
        end
        
        % movement definition as a whole or with on- and off-sets
        if defonoff
            % find local maxima
            if i == 1
                [pks,idx] = findlocmax(movwin{r,i},'onset');
            else
                [pks,idx] = findlocmax(movwin{r,i},'offset');
            end
            
            % select highest maxima within a time window
            twin = cfg.movdef(i).twin * data.fsample(r);
            cl = clusterdiff(idx,twin);
            n = max(cl);
            tmp = zeros(1,n);
            for c = 1:n
                clpks = pks;
                clpks(cl~=c) = -Inf;
                [~,j] = max(clpks);
                tmp(c) = idx(j);
            end
            
            % store the on-/offsets
            movwin{r,i} = tmp;
        end
        
    end
    
    % match on- and offsets
    if defonoff
        % get indeces
        idxon = movwin{r,1};
        idxoff = movwin{r,2};
        mov = false(size(data.time{r}));
        
        % loop over on- and off-sets until all are used
        while ~isempty(idxon) && ~isempty(idxoff)
            
            % discard outliers
            idxon = idxon(idxon<max(idxoff));
            idxoff = idxoff(idxoff>min(idxon));
            
            % onset
            ion = max(idxon(idxon<min(idxoff)));
            idxon = idxon(idxon>ion);
            
            % offset
            ioff = min(idxoff(idxoff>ion));
            idxoff = idxoff(idxoff>ioff);
            
            mov(ion:ioff) = true;
        end
        
        % go from the sample to the time domain
        movement{r} = logic2time(mov,data.time{r});
    else
        % make the movement window a logical array
        movwin{r,1} = movwin{r,1}>0;
        
        % go from the sample to the time domain
        movement{r} = logic2time(movwin{r,1},data.time{r});
    end
    
end


% update dataset
if isempty(regexpi(cfg.dataset,'_mov','once'))
    if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
    cfg.dataset = [cfg.dataset '_MOV'];
end

% update processing directory
cfg = km_setcfg(cfg,'dirana');
cfg.dir.proc = cfg.dir.ana;

% return movement definition or store in data and configuration structure
if nargout == 1
    varargout = {movement};
else
    if length(cfg.movdef) == 1
        cfg.movdef.movement = movement;
    end
    cfg.movement	= movement;
    data.movement	= movement;
    varargout       = {cfg,data};
end
%--------------------------------------------------------------------------


% function movdef_pos
%----------------------------------------
function mov = movdef_pos(cfg,data,r)

% select markers
markersel	= km_labelselection(cfg.marker,data.label);
markersel	= ismember(data.label,markersel);
if ~isfalse(cfg.refmarker)
    refsel 	= km_labelselection(cfg.refmarker,data.label);
    refsel	= ismember(data.label,refsel);
else
    refsel = [];
end

% get axis
if isfield(cfg,'crit') && ischar(cfg.crit)
    ax = regexp(cfg.crit,'(?<![a-zA-Z])[a-zA-Z](?![a-zA-Z])','match');
    ax = unique(ax);
elseif isfield(cfg,'axis')
    ax = cfg.axis;
else
    ax = data.axis;
end

% get data
dat = data.pos{r}(markersel,:,:);
if ~isempty(refsel)
    refpos = data.pos{r}(refsel,:,:);
    dat = dat - repmat(refpos,[size(dat,1) 1 1]);
end

% get axis data
for a = 1:length(ax)
    axissel = km_labelselection(ax{a},data.axis);
    axissel = ismember(data.axis,axissel);
    axpos.(ax{a}) = reshape(dat(:,axissel,:),[size(dat,1) 1 size(dat,3)]);
end

% apply criterium
if ischar(cfg.crit)
    for a = 1:length(ax)
        eval(sprintf('%s = axpos.(ax{a});',ax{a}));
    end
    mov = eval(cfg.crit);
else
    mov = true(size(dat,1),1,size(dat,3));
    for a = 1:length(ax)
        if length(cfg.crit) == 2
            mov = mov & axpos.(ax{a}) > min(cfg.crit) & axpos.(ax{a}) < max(cfg.crit);
        else
            mov = mov & axpos.(ax{a}) > cfg.crit;
        end
    end
end

% the data may be used in cfg.fun, so keep only selected axis
axissel = km_labelselection(ax,data.axis);
axissel = ismember(data.axis,axissel);
dat = dat(:,axissel,:);

% lookaround in a specified time window
if isfield(cfg,'lookaround') && ~isfalse(cfg.lookaround)
    mov = lookaround(cfg.lookaround,data,r,mov);
end

% evaluate objective function or not
if isfield(cfg,'fun') && ~isfalse(cfg.fun)
    mov = evalfun(dat,mov,cfg.fun);
end

% combine markers and axes
mov = combmov(mov,cfg.comb);

% function movdef_deriv
%----------------------------------------
function mov = movdef_deriv(cfg,data,r,pname)

% select data
if strcmpi(pname,'vel')
    markersel	= km_labelselection(cfg.marker,data.label);
    markersel	= ismember(data.label,markersel);
    axissel     = km_labelselection(cfg.axis,data.derivaxis);
    axissel     = ismember(data.derivaxis,axissel);
    dat = data.vel{r}(markersel,axissel,:);
elseif strcmpi(pname,'acc')
    markersel	= km_labelselection(cfg.marker,data.label);
    markersel	= ismember(data.label,markersel);
    axissel     = km_labelselection(cfg.axis,data.derivaxis);
    axissel     = ismember(data.derivaxis,axissel);
    dat = data.acc{r}(markersel,axissel,:);
elseif strcmpi(pname,'gripvel')
    markersel	= km_labelselection(cfg.marker,data.grip.label);
    markersel	= ismember(data.grip.label,markersel);
    axissel     = km_labelselection(cfg.axis,data.grip.aptaxis);
    axissel     = ismember(data.grip.aptaxis,axissel);
    dat = data.grip.aptvel{r}(markersel,axissel,:);
end

% take absolute value or not
origdat = dat;
if ~isfalse(cfg.abs),   dat = abs(dat); end

% apply criterium
if ischar(cfg.crit)
    mov = eval(cfg.crit);
elseif length(cfg.crit) == 2
    mov = dat > min(cfg.crit) & dat < min(cfg.crit);
else
    mov = dat > cfg.crit;
end

% slide borders to local minimum or criterium
if isequal(cfg.slide,0) || ~isfalse(cfg.slide)
    mov = km_slidesel(origdat,mov,cfg.slide,cfg.consecflg);
end

% lookaround in a specified time window
if isfield(cfg,'lookaround') && ~isfalse(cfg.lookaround)
    mov = lookaround(cfg.lookaround,data,r,mov);
end

% evaluate objective function or not
if isfield(cfg,'fun') && ~isfalse(cfg.fun)
    mov = evalfun(dat,mov,cfg.fun);
end

% combine markers and axes
mov = combmov(mov,cfg.comb);

% function movdef_grip
%----------------------------------------
function mov = movdef_grip(cfg,data,r,pname)

% select data
markersel	= km_labelselection(cfg.marker,data.grip.label);
markersel	= ismember(data.grip.label,markersel);
axissel     = km_labelselection(cfg.axis,data.grip.aptaxis);
axissel     = ismember(data.grip.aptaxis,axissel);
if strcmpi(pname,'gripapt')
    dat = data.grip.apt{r}(markersel,axissel,:);
elseif strcmpi(pname,'gripacc')
    dat = data.grip.aptacc{r}(markersel,axissel,:);
end

% apply criterium
if ischar(cfg.crit)
    mov = eval(cfg.crit);
elseif length(cfg.crit) == 2
    mov = dat > min(cfg.crit) & dat < max(cfg.crit);
else
    mov = dat > cfg.crit;
end

% lookaround in a specified time window
if isfield(cfg,'lookaround') && ~isfalse(cfg.lookaround)
    mov = lookaround(cfg.lookaround,data,r,mov);
end

% evaluate objective function or not
if isfield(cfg,'fun') && ~isfalse(cfg.fun)
    mov = evalfun(dat,mov,cfg.fun);
end

% combine markers and axes
mov = combmov(mov,cfg.comb);

% function movdef_log
%----------------------------------------
function mov = movdef_log(cfg,data,r,log)

% get logfile variable
idx = ~cellfun(@isempty,regexp(log.vars,cfg.var));
logdat = cellfun(@(x,i) x(:,i),log.data,repmat({idx},size(log.data)),'UniformOutput',false);

% find log indeces
lognr = find(log.runnr == data.run(r));
tim = data.time{r};
if ~isempty(lognr)
    tmp = logdat{lognr}(logdat{lognr}>=min(tim) & logdat{lognr}<=max(tim));
    idx = nearest(tim,tmp);
else
    error('KM:MissingLogRun','the data contains a run [%s] that is not present in the logfile',num2str(lognr));
    %idx = [];
end

% create mov vector
mov = zeros(size(tim));
mov(idx) = 1;

% convolve with specified function
switch cfg.conv
    case 'betapdfmaxcentered'
        % maximum centered beta probability density function
        b = betapdfmaxcentered(cfg.betapdf.a,cfg.betapdf.b,data.fsample(r),cfg.twin);
    case 'betapdf'
        % get x values
        bx = 0:(1/(data.fsample(r)*cfg.twin)):1;
        % beta probability density function
        b = betapdf(bx,cfg.betapdf.a,cfg.betapdf.b);
    case 'twin'
        error('FIXME:I am too lazy to create a vector of ones the size of cfg.twin');
        b = [1 1 1];
    case 'pre'
        error('FIXME:I am too lazy to create a vector of ones and zeros the size of cfg.twin');
        b = [1 1 0];
    case 'post'
        error('FIXME:I am too lazy to create a vector of zeros and ones the size of cfg.twin');
        b = [0 1 1];
    case 'no'
        b = 1;
    otherwise
        b = eval(cfg.conv);
        %warning('KM:MovDef:UnsupportedLogConv','This convolution setting [%s] is not supported',cfg.conv);
        %b = 1;
end
mov = conv(mov,b,'same');

% evaluate objective function or not
if isfield(cfg,'fun') && ~isfalse(cfg.fun)
    mov = evalfun([],mov,cfg.fun);
end
mov = mov(:);

% function movdef_probwin
%----------------------------------------
function mov = movdef_probwin(cfg,tim,movwin)

% find on and offsets of timewindow
win = prod(movwin(:,cfg.paramidx),2);
win = win>0;
mov = single(win(:))';
win = logic2idx(win);

% select or combine windows
switch lower(cfg.winsel)
    case {'minmax','comb'}        
        win = [min(win(:)) max(win(:))];
    otherwise
        if ~isfalse(cfg.winsel) && isnumeric(cfg.winsel)
            % limit selection to number of available windows
            winsel = cfg.winsel(cfg.winsel <= size(win,1));
            % select windows
            win = win(winsel,:);
        end
end

% loop over time window
for i = 1:size(win,1)
    
    % should the probfun operate relative to the whole window, or only to
    % the onset of the window?
    switch lower(cfg.funref)
        case {'on','onset','constant'}
            x = tim(win(i,1):win(i,2)) - tim(win(i,1));
        otherwise
            % samples of current window
            n = win(i,2)-win(i,1);
            x = linspace(0,1,n+1);
    end
    
    % data window
    y = mov(win(i,1):win(i,2));
    
    % convolve with specified function
    if isfield(cfg,'conv') && ~isfalse(cfg.conv)
        switch cfg.conv
            case 'betapdfmaxcentered'
                % maximum centered beta probability density function
                b = betapdfmaxcentered(cfg.betapdf.a,cfg.betapdf.b,n,1);
            case 'betapdf'
                % beta probability density function
                b = betapdf(x,cfg.betapdf.a,cfg.betapdf.b);
            otherwise
                warning('KM:MovDef:UnsupportedProbWinConv','This convolution setting [%s] is not supported',cfg.conv);
        end
        y = conv(y,b,'same');
    end
    
    % evaluate a specified function
    if isfield(cfg,'fun') && ~isfalse(cfg.fun)
        y = eval(cfg.fun);
    end
    
    % restore data window
    mov(win(i,1):win(i,2)) = y;
end

% function movdef_compwin
%----------------------------------------
function mov = movdef_compwin(cfg,data,r)

% select data
if strcmpi(cfg.param,'pos')
    % select markers
    markersel	= km_labelselection(cfg.marker,data.label);
    markersel	= ismember(data.label,markersel);
    if ~isfalse(cfg.refmarker)
        refsel 	= km_labelselection(cfg.refmarker,data.label);
        refsel	= ismember(data.label,refsel);
    else
        refsel = [];
    end
    % select axes
    axissel     = km_labelselection(cfg.axis,data.axis);
    axissel     = ismember(data.axis,axissel);
    % get data
    dat = data.pos{r}(markersel,axissel,:);
    if ~isempty(refsel)
        refpos = data.pos{r}(refsel,axissel,:);
        dat = dat - repmat(refpos,[size(dat,1) 1 1]);
    end
elseif ~isempty(regexp(cfg.param,'vel|acc','once'))
    markersel	= km_labelselection(cfg.marker,data.label);
    markersel	= ismember(data.label,markersel);
    axissel     = km_labelselection(cfg.axis,data.derivaxis);
    axissel     = ismember(data.derivaxis,axissel);
    if strcmpi(cfg.param,'vel')
        dat = data.vel{r}(markersel,axissel,:);
    elseif strcmpi(cfg.param,'acc')
        dat = data.acc{r}(markersel,axissel,:);
    end
elseif ~isempty(regexp(cfg.param,'grip','once'))
    markersel	= km_labelselection(cfg.marker,data.grip.label);
    markersel	= ismember(data.grip.label,markersel);
    axissel     = km_labelselection(cfg.axis,data.grip.aptaxis);
    axissel     = ismember(data.grip.aptaxis,axissel);
    if strcmpi(cfg.param,'gripvel')
        dat = data.grip.aptvel{r}(markersel,axissel,:);
    elseif strcmpi(pname,'gripapt')
        dat = data.grip.apt{r}(markersel,axissel,:);
    elseif strcmpi(pname,'gripacc')
        dat = data.grip.aptacc{r}(markersel,axissel,:);
    end
else
    error('parameter ''%s'' is not supported (yet)',cfg.param);
end

% transform time windows to sample windows
fs = data.fsample(r);
win1 = sort(round(cfg.twin1*fs));
win2 = sort(round(cfg.twin2*fs));
n1 = win1(2)-win1(1) + 1;
n2 = win2(2)-win2(1) + 1;

% loop over every dimension
mov = zeros(size(dat));
for m = 1:size(dat,1)
    for a = 1:size(dat,2)
        
        % add nan-padding to prevent edge artifacts
        tmp = [nan(n1+n2,1); squeeze(dat(m,a,:)); nan(n1+n2,1)];
        
        % check function to perform on windows
        if ~strcmpi(cfg.compfun,'median')
            error('for now only comparison of medians is allowed');
        end
        
        % use medfilt to obtain the median of the windows
        val1 = medfilt2(tmp,[n1,1]);
        if n2 == n1,
            val2 = val1;
        else
            val2 = medfilt2(tmp,[n2,1]);
        end
        
        % calculate relative window shift
        shift1 = floor(n1/2 - 0.5) + win1(1);
        shift2 = floor(n2/2 - 0.5) + win2(1);
        
        % pad values with nan to shift the time windows
        if shift1 < 0
            val1 = [nan(abs(shift1),1); val1(1:end+shift1)];
        else
            val1 = [val1(shift1:end); nan(abs(shift1),1)];
        end
        if shift2 < 0
            val2 = [nan(abs(shift2),1); val2(1:end+shift2)];
        else
            val2 = [val2(1+shift2:end); nan(abs(shift2),1)];
        end
        
        % cut nan-padding (previously added to prevent edge artifacts)
        idx = (n1+n2+1):(length(tmp)-n1-n2);
        %tmp = tmp(idx);
        val1 = val1(idx);
        val2 = val2(idx);
        
        % FIXME: place this option within the cfg.fun!!!
        % take absolute value or not
        if ~isfalse(cfg.abs)
            val1 = abs(val1);
            val2 = abs(val2);
        end
        
        % FIXME: place this option within the cfg.fun!!!
        % apply criterium
        mov(m,a,:) = (val2-val1) > cfg.crit;
    end
end

% replace nans with zeros
mov(isnan(mov)) = 0;

% lookaround in a specified time window
if isfield(cfg,'lookaround') && ~isfalse(cfg.lookaround)
    mov = lookaround(cfg.lookaround,data,r,mov);
end

% evaluate objective function or not
if isfield(cfg,'fun') && ~isfalse(cfg.fun)
    mov = evalfun(dat,mov,cfg.fun);
end

% combine markers and axes
mov = combmov(mov,cfg.comb);

% function evalfun
%----------------------------------------
function mov = evalfun(dat,mov,fun)
if strcmpi(fun,'mov'),  return;	end

% prepare data descriptives
if ~isempty(regexp(fun,'meandat','once'))
    meandat = repmat(nanmean(dat,3),[1 1 size(dat,3)]);
end
if ~isempty(regexp(fun,'meddat','once'))
    meddat = repmat(nanmedian(dat,3),[1 1 size(dat,3)]);
end
if ~isempty(regexp(fun,'maxabsdat','once'))
    maxabsdat = repmat(max(abs(dat),[],3),[1 1 size(dat,3)]);
end
if ~isempty(regexp(fun,'maxdat','once'))
    maxdat = repmat(max(dat,[],3),[1 1 size(dat,3)]);
end
if ~isempty(regexp(fun,'mindat','once'))
    mindat = repmat(min(dat,[],3),[1 1 size(dat,3)]);
end
if ~isempty(regexp(fun,'maxmov','once'))
    maxmov = max(mov);
end
if ~isempty(regexp(fun,'minmov','once'))
    minmov = min(mov);
end
mov = eval(fun);

% function combmov
%----------------------------------------
function mov = combmov(mov,comb)

% quick return
if isfalse(comb)
    return;
end

% force mov double
if ~any(strcmpi(comb,{'any','all'}))
    mov = double(mov);
end

% see if zeros need to be replaced
nonzero = false;
if ~isempty(regexp(comb,'nonzero','once'))
    nonzero = true;
    comb = regexprep(comb,'nonzero','');
    mov(mov==0) = nan;
end

% combine markers and axes
switch comb
    case {'min','max'}
        mov = squeeze(feval(comb,feval(comb,mov,[],1),[],2));
    otherwise
        mov = squeeze(feval(comb,feval(comb,mov,1),2));
end

% if zeros need to be replaced
if nonzero
    mov(isnan(mov)) = 0;
end

% limit mov to the [0 1] range
mov(mov>1) = 1;
mov(mov<0) = 0;

% function lookaround
%----------------------------------------
function mov = lookaround(cfg,data,r,mov)

mov = single(mov);
fs = data.fsample(r);
idx = -cfg.twin*fs;
maxidx = max(abs(idx));
tmp = -maxidx:maxidx;
twin = zeros(length(tmp),1);
twin(tmp>=min(idx) & tmp<=max(idx)) = 1;
for m = 1:size(mov,1)
    for a = 1:size(mov,2)
        mov(m,a,:) = conv(squeeze(mov(m,a,:)),twin,'same')/sum(twin);
    end
end
switch lower(cfg.find)
    case 'all', mov(mov<1) = 0;
    case 'any', mov(mov>0) = 1;
    otherwise,  error('lookaround.find: use only ''all'' or ''any''');
end