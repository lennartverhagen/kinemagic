function cfg = km_setcfg(cfg,task,method)
%--------------------------------------------------------------------------
%
% See also KM_SETTASK
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% check input
if nargin < 1
    cfg = []; return;
end
if nargin < 2
    if ischar(cfg) || iscell(cfg)
        task = cfg;
        cfg = [];
    else
        task = {'basic'};
    end
end
task = km_settask(task);
if nargin < 3
    method = 'free';
end

% (re)set the cfg.save field with every task (except when km_setcfg is
% called by km_save)
[ST,~] = dbstack(1);
if ~strcmpi(ST(1).name,'km_save')
    cfg = setcfg_save(cfg);
end

% list of tasks that operate on the first level of the configuration structure
listbasic   = {'basic','dataset','save'};
listio      = {'subj','sess','dir','dirdata','dirlog','dirpreproc','dirana','dirreport','dirproc'};
listproc    = {'read','readlog','alignlog','trlrej'};
listplot    = {'plot'};
cfg1levlist = [listbasic listio listproc listplot];

% evalute appropriate subfunction
if ~iscell(task),   task = {task};  end
for t = 1:length(task)
    
    % check if the task specific function exists
    funname = ['setcfg_' task{t}];
    if exist(funname,'file')==2
        
        % organize input
        % see if the function operates on the first level of cfg or not
        if any(strcmpi(task{t},cfg1levlist))
            varargin = {cfg};
        elseif ~isfield(cfg,task{t})
            varargin = {[]};
        else
            varargin = {cfg.(task{t})};
        end
        if nargin(funname) > 1
            varargin{2} = method;
        end       
        
        % evaluate subfunction and assign output
        if any(strcmpi(task{t},cfg1levlist))
            cfg = feval(funname,varargin{:});
        else
            cfg.(task{t}) = feval(funname,varargin{:});
        end
        
    % not supported: error
    elseif strcmpi(method,'strict')
        error('task ''%s'' not recognized',task{t});
        
    % not supported: warning
    elseif strcmpi(method,'loose')
        warning('KM:NotSupported','%s: This task (%s) is not supported',mfilename,task{t});
        
    end
end
%--------------------------------------------------------------------------


%% SUBFUNCTIONS
%--------------------------------------------------------------------------

% function setcfg_basic
%----------------------------------------
function cfg = setcfg_basic(cfg) %#ok<*DEFNU>
cfg = setcfg_dataset(cfg);
cfg = setcfg_save(cfg);
if ~isfield(cfg,'feedback'),        cfg.feedback    = 'no';                 end     % 'no', 'yes'

% function setcfg_dataset
%----------------------------------------
function cfg = setcfg_dataset(cfg)
if ~isfield(cfg,'dataset'),         cfg.dataset     = '';                   end     % 'MASK', 'FILT'
if ~isempty(cfg.dataset) && ~strcmpi(cfg.dataset(1),'_')
    cfg.dataset = ['_' cfg.dataset];
end

% function setcfg_save
%----------------------------------------
function cfg = setcfg_save(cfg)
if ~isfield(cfg,'save'),            cfg.save        = 'no';                 end     % 'no', 'yes'
if ~isfalse(cfg.save)
    cfg.save = 'yes';
    %cfg = setcfg_dirproc(cfg);
    %% FIXME: the line above should be commented out: km_plot calls
    %% setcfg_dirana (via a very indirect route) erroneously, before
    %% cfg.dir.ana could be prepend with the root and the subject name
end

% function setcfg_exp
%----------------------------------------
function cfg = setcfg_exp(cfg)
if ~isfield(cfg,'exp') || isfalse(cfg.exp),         error('no experiment name: cfg.exp');	end

% function setcfg_idxr
%----------------------------------------
function cfg = setcfg_idxr(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'dataset'),         cfg.dataset     = 'IDXR';               end
cfg = setcfg_dataset(cfg);
if ~isfield(cfg,'listwise')
    if isfield(cfg,'param')
        cfg.listwise = cfg.param;
    elseif isfield(cfg,'mode')
        cfg.listwise = cfg.mode;
    else
        cfg.listwise = {'apriori','rt','mt','td','mv','pv','mga','ega','dego'};
    end
end
if ~iscell(cfg.listwise),           cfg.listwise    = {cfg.listwise};       end
% cfg.mode is only supported for legacy reasons
cfg.mode = cfg.listwise;
% cfg.param is only supported for legacy reasons
cfg.param = cfg.listwise;
if ~isfield(cfg,'pairwise'),        cfg.pairwise    = {'none'};             end
if ~iscell(cfg.pairwise),           cfg.pairwise    = {cfg.pairwise};       end


%% PREPROCESSING
%--------------------------------------------------------------------------

% function setcfg_read
%----------------------------------------
function cfg = setcfg_read(cfg)
cfg = setcfg_subj(cfg);
cfg = setcfg_sess(cfg);
cfg = setcfg_dirdata(cfg);
if ~isfield(cfg,'read'),            cfg.read        = struct;               end
if isfalse(cfg.read),               return;                                 end
if ~isfield(cfg.read,'type'),       cfg.read.type   = 'polhemus';           end
if ~isfield(cfg.read,'expfun'),     cfg.read.expfun = 'no';                 end
if ~isfalse(cfg.read.expfun)
    if strcmpi(cfg.read.expfun,'exp')
        cfg = setcfg_exp(cfg);
        cfg.read.expfun = ['km_readfun_' cfg.exp];
    elseif exist(cfg.read.expfun,'file') ~= 2 && isempty(regexp(cfg.read.expfun,'^km_readfun_','once'))
        cfg.read.expfun = ['km_readfun_' cfg.read.expfun];
    end
end

% function setcfg_cutbreaks
%----------------------------------------
function cfg = setcfg_cutbreaks(cfg)
if isempty(cfg),                    cfg             = struct;               end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'maxbreak'),        cfg.maxbreak	= '3*samp';             end     % cut a run in chunks if samples are missing for more than maxbreak: time in seconds or '2*samp', '3*samp', etc. Strings are evaluated, the most common sample time is stored in the variable 'samp'
if ~isfield(cfg,'minchunk'),        cfg.minchunk	= 4;                    end     % throw a chunk away if it is shorter than minchunk: time in seconds or '500*samp', '1000*samp', etc. Strings are evaluated, the most common sample time is stored in the variable 'samp'

% function setcfg_readlog
%----------------------------------------
function cfg = setcfg_readlog(cfg)
if ~isfield(cfg,'readlog'),         cfg.readlog     = struct;               end
if isfalse(cfg.readlog),            return;                                 end
if istrue(cfg.readlog,'free'),      cfg.readlog     = struct;               end
if ~isfield(cfg.readlog,'expfun')
    if isfield(cfg.readlog,'fun')
        cfg.readlog.expfun = cfg.readlog.fun;
        cfg.readlog = rmfield(cfg.readlog,'fun');
    else
        cfg.readlog.expfun = 'exp';
    end
end
if strcmpi(cfg.readlog.expfun,'exp')
    cfg = setcfg_exp(cfg);
    cfg.readlog.expfun = ['km_readlogfun_' cfg.exp];
elseif exist(cfg.readlog.expfun,'file') ~= 2 && isempty(regexp(cfg.readlog.expfun,'^km_readlogfun_','once'))
    cfg.readlog.expfun = ['km_readlogfun_' cfg.readlog.expfun];
end
cfg = setcfg_subj(cfg);
cfg = setcfg_sess(cfg);
cfg = setcfg_dirlog(cfg);

% function setcfg_alignlog
%----------------------------------------
function cfg = setcfg_alignlog(cfg,method)
if ~isfield(cfg,'alignlog'),        cfg.alignlog     = struct;              end
if isfalse(cfg.alignlog),           return;                                 end
if strcmpi(method,'strict') && ~isfield(cfg,'log'),
    error('No log provided in cfg.log.');
end
if ~isfield(cfg.alignlog,'param'),  cfg.alignlog.param	= 'event';          end
if strcmpi(cfg.alignlog.param,'event')
    if ~isfield(cfg.alignlog,'evtcode'),    cfg.alignlog.evtcode	= 1;	end
    if ~isfield(cfg.alignlog,'logvar'),     error('No logvar provided.');   end
    if ~isfield(cfg.alignlog,'initidx'),    cfg.alignlog.initidx	= 1:3;	end
end
if ~isfield(cfg.alignlog,'err'),            cfg.alignlog.err        = 0.01;	end
if ~isfield(cfg.alignlog,'adjust'),  cfg.alignlog.adjust	= 'log';     	end
if strcmpi(cfg.alignlog.adjust,'log')
    if ~isfield(cfg.alignlog,'timevar'),    cfg.alignlog.timevar	= '^t_';end
end
% function setcfg_hemifieldcross
%----------------------------------------
function cfg = setcfg_hemifieldcross(cfg)
if isempty(cfg),                    cfg             = 'no';                 end     % 'no', 'loose', 'strict'

% function setcfg_axes
%----------------------------------------
function cfg = setcfg_axes(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'order'),           cfg.order       = 'no';                 end     % 'no', 'xyz', 'xzy', etc
if ~isfield(cfg,'mirror'),          cfg.mirror      = 'no';                 end     % 'no', 'x', 'yz', etc;
if ~isfield(cfg,'rotate'),          cfg.rotate      = 'no';                 end     % 'no' or rotation in degrees over [x y z] axes

% function setcfg_calcmarker
%----------------------------------------
function cfg = setcfg_calcmarker(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
for c = 1:length(cfg)
    if ~isfield(cfg(c),'type'),     cfg(c).type     = 'no';                 end     % 'no', 'mean'
    if ~isfield(cfg(c),'marker'),   cfg(c).marker	= {'all'};              end     % 'no', 'all', {'A','B'}, {'all','-A'}
    if ~iscell(cfg(c).marker),      cfg(c).marker	= {cfg(c).marker};      end     % marker needs to be a cell to streamline the marker selection
    if ~isfield(cfg(c),'label'),    cfg(c).label    = strcat(cfg(c).type,'_',cfg(c).marker{:});	end     % 'mean_AB'
end

% function setcfg_interp
%----------------------------------------
function cfg = setcfg_interp(cfg)
if isempty(cfg),                    cfg             = struct;               end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'method'),          cfg.method      = 'pchip';              end     % 'no', 'pchip', 'spline', type <help interp1> for more
if ~isfield(cfg,'freq'),            cfg.freq        = 'max';                end     % 'max' or any frequency in Hz, e.g. 200

% function setcfg_artfctdef
%----------------------------------------
function cfg = setcfg_artfctdef(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isstruct(cfg),
    tmp = cfg;  clear cfg;
    cfg.type	= tmp;
end
if ~isfield(cfg,'type'),            cfg.type        = 'no';                 end     % 'no', 'tms'

% function setcfg_maskart
%----------------------------------------
function cfg = setcfg_maskart(cfg)
if isempty(cfg),                    cfg             = 'no';                 end     % 'no', 'nan' , 'interp', 'nearest', 'linear', 'spline', 'pchip', 'cubic'
if strcmpi(cfg,'interp'),           cfg             = 'pchip';              end     % 

% function setcfg_filter
%----------------------------------------
function cfg = setcfg_filter(cfg)
if isempty(cfg),                    cfg             = struct;               end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'tseries') && isfield(cfg,'param'), cfg.tseries = cfg.param;end     % for legacy only
if ~isfield(cfg,'tseries'),         cfg.tseries     = {'pos','ori'};        end     % 'no', 'pos', 'vel', 'acc', 'ori', 'orivel', 'oriacc', 'gripapt', 'gripvel', 'gripacc', 'gripori', 'griporivel', 'griporiacc'
if ~iscell(cfg.tseries),            cfg.tseries     = {cfg.tseries};        end
if ~isfield(cfg,'mode'),            cfg.mode        = 'lowpass';            end     % 'lowpass', 'bandpass', 'highpass', 'stoppass'
if ~isfield(cfg,'type'),            cfg.type        = 'but';                end     % 'but', 'dft', 'cheby1', 'cheby2'
if ~isfield(cfg,'freq'),            cfg.freq        = 20;                   end     % frequency cutoff in hertz. bandpass and stoppass require two frequencies [low high]
if ~isfield(cfg,'dir'),             cfg.dir         = 'twopass';            end     % 'twopass', 'onepass', 'onepass-reverse'
if ~isfield(cfg,'ord')
    if strcmpi(cfg.type,'but')
        cfg.ord	= 6;                % the order of the filter, 6 for butterworth
    elseif strcmpi(cfg.type,'fir')
        cfg.ord	= 26;               % the order of the filter, 26 for finite impulse response
    end
end

% function setcfg_grip
%----------------------------------------
function cfg = setcfg_grip(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
% marker and label selection
if ~isfield(cfg,'marker'),          cfg.marker      = {'all'};              end     % 'all', {'A','B'}, {'all','-A'}
if ~iscell(cfg.marker),             cfg.marker      = {cfg.marker};         end     % marker needs to be a cell to streamline the marker selection
if ~isfield(cfg,'label'),           cfg.label       = {strcat('grip_',cfg.marker{:})};	end     % 'grip_AB' FIXME: this only works for single grip specs
if length(cfg.marker) > length(cfg.label),          cfg.marker = {cfg.marker};          end     % organize the marker selection per grip label
% aperture
if ~isfield(cfg,'apt'),             cfg.apt         = struct;               end
cfg.apt = setcfg_grip_apt(cfg.apt);
% orientation
if ~isfield(cfg,'ori'),             cfg.ori         = struct;               end
cfg.ori = setcfg_grip_ori(cfg.ori);

% function setcfg_grip_apt
%----------------------------------------
function cfg = setcfg_grip_apt(cfg)
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'axis'),            cfg.axis        = {'all'};              end     % {'all'}, {'xyz'}, {'xz', 'x'}, etc
if ~iscell(cfg.axis),               cfg.axis        = {cfg.axis};           end
if ~isfield(cfg,'abs'),             cfg.abs         = 'no';                 end     % take absolute value of single dimensional distances 'yes', or not 'no', please note Euclidian multi-dimension distances are always absolute

% function setcfg_grip_ori
%----------------------------------------
function cfg = setcfg_grip_ori(cfg)
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'axis'),            cfg.axis        = {'all'};              end     % {'all'}, {'xy', 'x'}, etc (but not 'xyz')
if ~iscell(cfg.axis),               cfg.axis        = {cfg.axis};           end
if ~isfield(cfg,'dir'),             cfg.dir         = {'ccw'};              end     % 'ccw' or 'cw' (degrees increase counter-clockwise, but you can flip them)
if ~iscell(cfg.dir),                cfg.dir         = {cfg.dir};            end
if ~isfield(cfg,'offset'),          cfg.offset      = 0;                    end     % offset in degrees (change the zero degrees orientation)

% function setcfg_deriv
%----------------------------------------
function cfg = setcfg_deriv(cfg)
if isempty(cfg),                    cfg             = struct;               end     
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'vel'),             cfg.vel         = 'yes';                end
if ~isfield(cfg,'acc'),             cfg.acc         = 'yes';                end
if ~isfield(cfg,'method'),          cfg.method      = 'gradient';           end     % 'gradient', 'diff'
if ~isfield(cfg,'axis'),            cfg.axis        = {'all'};              end     % {'all'}, {'xyz'}, {'xz', 'x'}, etc
if ~iscell(cfg.axis),               cfg.axis        = {cfg.axis};           end
if ~isfield(cfg,'abs'),             cfg.abs         = 'no';                 end     % take absolute value of single dimensional derivatives 'yes', or not 'no', please note Euclidian multi-dimension derivatives are always absolute


%% PROCESSING
%--------------------------------------------------------------------------

% function setcfg_movdef
%----------------------------------------
function cfg = setcfg_movdef(cfg)
if isempty(cfg),                    cfg             = struct;               end
if isfalse(cfg),                    return;                                 end
% loop over on- and offsets
for c = 1:length(cfg)
    if ~isfield(cfg(c),'feedback') || isempty(cfg(c).feedback)
        cfg(c).feedback	= 'no';
    end
    if ~isfalse(cfg(c).feedback) && strcmpi(cfg(c).feedback,'movplot')
        if isfield(cfg(c),'movplot')
            cfg(c) = setcfg_movplot(cfg(c));
        end
    end
    if ~isfield(cfg(c),'param') || isempty(cfg(c).param)
        cfg(c).param	= {'vel'};                                                  % 'vel', 'gripvel'
    end
    if ~iscell(cfg(c).param),       cfg(c).param    = {cfg(c).param};       end
    if length(cfg(c).param) > 1,	cfg(c).msi   	= 'yes';
    elseif ~isfield(cfg(c),'msi'),	cfg(c).msi    	= 'no';
    end
    if length(cfg)>1 && ~isfield(cfg(c),'twin')
        cfg(c).twin = 0.5;
    end
    
    ip = struct;
    for p = 1:length(cfg(c).param)
        pname = cfg(c).param{p};
        
        % counter
        if ~isfield(ip,pname)
            ip.(pname) = 1;
        else
            ip.(pname) = ip.(pname) + 1;
        end
        if ~isfield(cfg(c),pname),  cfg(c).(pname)	= struct;               end
        tcfg = cfg(c).(pname);
        if length(tcfg) < ip.(pname)
            if ~iscell(tcfg)
                tcfg = {tcfg(:)};
            end
            tcfg{ip.(pname)} = tcfg{1};
        end
        if ~iscell(tcfg)
            if ~isfalse(cfg(c).msi) && ~isfield(tcfg,'fun')
                error('no objective function specified to enter MSI algorithm for %s',pname);
            end
            tcfg = feval(['setcfg_movdef_' pname],tcfg);
            tcfg = setcfg_movdef_lookaround(tcfg,pname);
        else
            for i = 1:length(tcfg)
                if ~isfalse(cfg(c).msi) && ~isfield(tcfg{i},'fun')
                    error('no objective function specified to enter MSI algorithm for %s',pname);
                end
                tcfg{i} = feval(['setcfg_movdef_' pname],tcfg{i});
                tcfg{i} = setcfg_movdef_lookaround(tcfg{i},pname);
            end
        end
        cfg(c).(pname) = tcfg;
    end
end

% function setcfg_movdef_lookaround
%----------------------------------------
function cfg = setcfg_movdef_lookaround(cfg,pname)
if isfield(cfg,'twin') && length(cfg.twin) > 1
    warning('KM:twin2lookaround','rename cfg.%s.twin to cfg.%s.lookaround.twin',pname,pname);
    cfg.lookaround.twin = cfg.twin;
    cfg = rmfield(cfg,'twin');
end
if ~isfield(cfg,'lookaround')
    cfg.lookaround = 'no';
elseif ~isfalse(cfg.lookaround)
    if ~isfield(cfg.lookaround,'twin')
        error('please give a time window for the cfg.%s.lookaround.twin',pname);
    end
    if ~isfield(cfg.lookaround,'find')
        cfg.lookaround.find = 'all';
    end
    if ~any(strcmpi(cfg.lookaround.find,{'all','any'}))
        error('cfg.%s.lookaround.find can only be ''all'' or ''any''',pname);
    end
end

% function setcfg_movdef_vel
%----------------------------------------
function cfg = setcfg_movdef_vel(cfg)
if ~isfield(cfg,'marker'),      error('No marker specified for movdef');	end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end     % axis of the movdef parameter
if ~isfield(cfg,'crit'),            cfg.crit        = 0.05;                 end     % critical value
if ~isfield(cfg,'abs'),             cfg.abs         = 'yes';                end     % take the absolute values of the velocity parameter
if ~isfield(cfg,'slide'),           cfg.slide       = cfg.crit/5;           end     % let the on- and off-sets slide away till a local minimum or the slide criterium is reached
if ~isfield(cfg,'consecflg'),       cfg.consecflg   = 'onset';              end     % if a consecutive off- and on-set come together, move the 'onset', 'offset', of 'both' one sample away from the border
if ~isnumeric(cfg.slide) && istrue(cfg.slide,'free')
    cfg.slide	= 0;
end
cfg = setcfg_movdef_comb(cfg);

% function setcfg_movdef_acc
%----------------------------------------
function cfg = setcfg_movdef_acc(cfg)
cfg = setcfg_movdef_vel(cfg);

% function setcfg_movdef_pos
%----------------------------------------
function cfg = setcfg_movdef_pos(cfg)
if ~isfield(cfg,'marker'),      error('No marker specified for movdef');	end
if ~isfield(cfg,'refmarker'),       cfg.refmarker   = 'no';                 end     % reference marker, serving as a zero point
if ~isfield(cfg,'crit'),            cfg.crit        = -Inf;                 end     % critical value
if isfield(cfg,'axis') && ~iscell(cfg.axis), cfg.axis = {cfg.axis};         end
cfg = setcfg_movdef_comb(cfg);

% function setcfg_movdef_gripapt
%----------------------------------------
function cfg = setcfg_movdef_gripapt(cfg)
if ~isfield(cfg,'marker'),      error('No marker specified for movdef');	end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end     % axis of the movdef parameter
if ~isfield(cfg,'crit'),            cfg.crit        = -Inf;                 end     % critical value
cfg = setcfg_movdef_comb(cfg);

% function setcfg_movdef_gripvel
%----------------------------------------
function cfg = setcfg_movdef_gripvel(cfg)
cfg = setcfg_movdef_vel(cfg);

% function setcfg_movdef_gripacc
%----------------------------------------
function cfg = setcfg_movdef_gripacc(cfg)
if ~isfield(cfg,'marker'),      error('No marker specified for movdef');	end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end     % axis of the movdef parameter
if ~isfield(cfg,'abs'),             cfg.abs         = 'yes';                end     % take the absolute values of the velocity parameter
cfg = setcfg_movdef_comb(cfg);

% function setcfg_movdef_logfile
%----------------------------------------
function cfg = setcfg_movdef_logfile(cfg)
if ~isfield(cfg,'var'),  	error('No variable name specified for movdef'); end

% function setcfg_movdef_probwin
%----------------------------------------
function cfg = setcfg_movdef_probwin(cfg)
if ~isfield(cfg,'paramidx'),error('no parameter window indices provided');	end
if ~isfield(cfg,'winsel'),          cfg.winsel      = 'all';                end     % 'minmax' or an integer index
if ~isfield(cfg,'funref'),      	cfg.funref  	= 'win';                end     % 'win', 'onset'

% function setcfg_movdef_compwin
%----------------------------------------
function cfg = setcfg_movdef_compwin(cfg)
if ~isfield(cfg,'param'),   error('No parameter specified for movdef_compwin');	end
if ~isfield(cfg,'marker'),  error('No marker specified for movdef_compwin');    end
if ~isfield(cfg,'refmarker'),       cfg.refmarker = 'no';                   end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end     % axis of the movdef parameter
if ~isfield(cfg,'twin1'),	error('No time-window-1 specified for movdef_compwin');	end
if ~isfield(cfg,'twin2'),	error('No time-window-2 specified for movdef_compwin');	end
if ~isfield(cfg,'compfun'),         cfg.compfun     = 'median';             end     % change
if ~isfield(cfg,'crit'),            cfg.crit        = 0.05;                 end     % incorporate in fun
if ~isfield(cfg,'abs'),             cfg.abs         = 'no';                 end     % incorporate in fun
if ~isfield(cfg,'fun'),             cfg.fun         = 'mov';                end     % change
cfg = setcfg_movdef_comb(cfg);

% function setcfg_movdef_comb
%----------------------------------------
function cfg = setcfg_movdef_comb(cfg)
if ~isfield(cfg,'comb')
    if isfield(cfg,'fun') && ~isfalse(cfg.fun)
        if strcmpi(cfg.fun,'mov')
            cfg.comb = 'prod';
        else
            cfg.comb = 'mean';
        end
    else
        cfg.comb = 'any';
    end
end

% function setcfg_movcheck
%----------------------------------------
function cfg = setcfg_movcheck(cfg)
if isempty(cfg),                    cfg             = struct;               end
if isfalse(cfg),                    return;                                 end
% loop over multiple options
for c = 1:length(cfg)
    % force parameters of next options to match the first
    if c > 1
        cfg(c).param = cfg(1).param;
    end
    if ~isfield(cfg(c),'refmarker') || isempty(cfg(c).refmarker)
        cfg(c).refmarker = 'no';
    end
    if ~isfield(cfg(c),'param') || isempty(cfg(c).param)
        cfg(c).param = {'dur'};
    end
    if ~iscell(cfg(c).param),          	cfg(c).param        = {cfg.param};	end
    for p = 1:length(cfg(c).param)
        pname = lower(cfg(c).param{p});
        if ~isfield(cfg(c),pname),
            cfg(c).(pname)	= struct;
        end
        switch lower(pname)
            case 'dur'
                % durations
                cfg(c).(pname) = setcfg_check_dur(cfg(c).(pname));
            case 'bvel'
                % begin velocity
                cfg(c).(pname) = setcfg_check_bvel(cfg(c).(pname));
            case 'peakvel'
                % peak velocity
                cfg(c).(pname) = setcfg_check_peakvel(cfg(c).(pname));
            case {'bga','ega'}
                % end grip aperture
                cfg(c).(pname) = setcfg_check_apt(cfg(c).(pname));
            case {'pos','bpos','epos'}
                % position
                cfg(c).(pname) = setcfg_check_pos(cfg(c).(pname));
        end
    end
    idx = ismember(cfg(c).param,'dur');
    cfg(c).param = [cfg(c).param(find(idx,1,'first')) cfg(c).param(~idx)];
end
for c = 1:length(cfg)
    flds = fieldnames(cfg(c));
    for f = length(flds)
        if isempty(cfg(c).(flds{f}))
            error('KM:SetCFG:MovCheck','All fields of cfg.movcheck must be specified in all options.');
        end
    end
end

% function setcfg_check_dur
%----------------------------------------
function cfg = setcfg_check_dur(cfg)
if ~isfield(cfg,'off'),             cfg.off         = 0.01;                 end     % minimum duration of no movement (rest/break) in seconds
if ~isfield(cfg,'on'),              cfg.on          = [0.3 3];              end     % minimum (and maximum) duration of movement in seconds

% function setcfg_check_bvel
%----------------------------------------
function cfg = setcfg_check_bvel(cfg)
if istrue(cfg,'free'),              cfg             = struct;               end
if ~isfield(cfg,'abs'),             cfg.abs         = 'yes';                end     % take the absolute values of the velocity parameter
if ~isfield(cfg,'marker'),      error('No marker specified for bvel');      end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end     % axis of the velocity parameter
if ~isfield(cfg,'crit'),            cfg.crit        = 0.05;                 end     % minimum (and maximum) peak velocity of a movement in m/s

% function setcfg_check_peakvel
%----------------------------------------
function cfg = setcfg_check_peakvel(cfg)
if istrue(cfg,'free'),              cfg             = struct;               end
if ~isfield(cfg,'abs'),             cfg.abs         = 'yes';                end     % take the absolute values of the velocity parameter
if ~isfield(cfg,'marker'),      error('No marker specified for peakvel');	end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end     % axis of the velocity parameter
if ~isfield(cfg,'crit'),            cfg.crit        = 0.5;                  end     % minimum (and maximum) peak velocity of a movement in m/s

% function setcfg_check_ega
%----------------------------------------
function cfg = setcfg_check_apt(cfg)
if istrue(cfg,'free'),              cfg             = struct;               end
if ~isfield(cfg,'marker'),      error('No marker specified for grip apt');	end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end     % axis of the end grip aperture parameter
if ~isfield(cfg,'crit'),            cfg.crit        = [-Inf 0.15];          end     % minimum (and maximum) end grip aperture in meters

% function setcfg_check_pos
%----------------------------------------
function cfg = setcfg_check_pos(cfg)
if istrue(cfg,'free'),              cfg             = struct;               end
if ~isfield(cfg,'marker'),      error('No marker specified for position');	end
if ~isfield(cfg,'refmarker'),       cfg.refmarker   = 'no';                 end
if ~isfield(cfg,'type'),            cfg.type        = 'incl';               end
if ~isfield(cfg,'crit'),        error('No criterium specified for position');	end
if ~ischar(cfg.crit)
    if ~isfield(cfg,'axis'),    error('No axis specified for position');	end
    if length(cfg.axis) > 1,    error('Use evaluable strings as criteria when using more than one axis');   end
end

% function setcfg_movplot
%----------------------------------------
function cfg = setcfg_movplot(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if istrue(cfg,'free'),              cfg             = struct;               end
if ~isfield(cfg,'param'),           cfg.param       = {'vel'};              end
if ~iscell(cfg.param),              cfg.param       = {cfg.param};          end
if isfield(cfg,'marker') && iscell(cfg.marker),	cfg.marker	= cfg.marker{1};end
if isfield(cfg,'axis') && iscell(cfg.axis),     cfg.axis	= cfg.axis{1};  end
for p = 1:length(cfg.param)
    pname = cfg.param{p};
    if ~isfield(cfg,pname), continue;   end
    if isfield(cfg.(pname),'marker') && iscell(cfg.(pname).marker)
        cfg.(pname).marker	= cfg.(pname).marker{1};
    end
    if isfield(cfg.(pname),'axis') && iscell(cfg.(pname).axis)
        cfg.(pname).axis	= cfg.(pname).axis{1};
    end
end
if ~isfield(cfg,'xlim'),            cfg.xlim        = [-0.5 5];             end
if ~isfield(cfg,'ylim'),            cfg.ylim        = 'auto';               end

% function setcfg_trldef
%----------------------------------------
function cfg = setcfg_trldef(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'param'),       error('No trial def parameter specified');  end
if strcmpi(cfg.param,'log'),        cfg         = setcfg_trldef_log(cfg);   end
if ~isfield(cfg,'catruns'),         cfg.catruns     = 'no';                 end

% function setcfg_trldef_log
%----------------------------------------
function cfg = setcfg_trldef_log(cfg)
if ~isfield(cfg,'beg'),             cfg.beg         = 't_beg';              end
if ~isfield(cfg,'stim'),            cfg.stim        = cfg.beg;              end
if ~isfield(cfg,'end'),             cfg.end         = 't_end';              end
if ~isfield(cfg,'prepad'),          cfg.prepad      = 0;                    end     % prepad in seconds relative to begin of trial
if ~isfield(cfg,'postpad'),         cfg.postpad     = 0;                    end     % postpad in seconds relative to end of trial
if ~isfield(cfg,'offset'),          cfg.offset      = 0;                    end     % offset in seconds relative to stimulus

% function setcfg_cuttrials
%----------------------------------------
function cfg = setcfg_cuttrials(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'align'),           cfg.align       = 'no';                 end
if ~isfield(cfg,'timevar'),         cfg.timevar     = 'no';                 end

% function setcfg_adjusttime
%----------------------------------------
function cfg = setcfg_adjusttime(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'timevar'),         cfg.timevar     = 'no';                 end

% function setcfg_trlcheck
%----------------------------------------
function cfg = setcfg_trlcheck(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'refmarker'),       cfg.refmarker   = 'no';                 end
if ~isfield(cfg,'param'),           cfg.param       = {};                   end
if ~iscell(cfg.param),              cfg.param       = {cfg.param};          end
for p = 1:length(cfg.param)
    pname = lower(cfg.param{p});
    if ~isfield(cfg,pname),
        cfg.(pname)	= struct;
    end
    switch lower(pname)
        case 'bvel'
            % begin velocity
            cfg.(pname) = setcfg_check_bvel(cfg.(pname));
        case {'bga','ega'}
            % end grip aperture
            cfg.(pname) = setcfg_check_apt(cfg.(pname));
        case {'pos','bpos','epos'}
            % position
            cfg.(pname) = setcfg_check_pos(cfg.(pname));
    end
end

% function setcfg_movselect
%----------------------------------------
function cfg = setcfg_movselect(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'movidx'),          cfg.movidx      = [];                   end
if ~isfield(cfg,'nmov') || isfalse(cfg.nmov)
    if ~isfalse(cfg.movidx)
        cfg.nmov = length(cfg.movidx);
    else
        cfg.nmov = 1;
    end
end
cfg.nmov = [min(cfg.nmov) max(cfg.nmov)];
if ~isfield(cfg,'param'),       	cfg.param       = 'index';              end
if ~iscell(cfg.param),              cfg.param       = {cfg.param};          end
if length(cfg.param) > 1,           error('only a single ordering parameter is allowed (the combination is multiple parameters is not implemented yet).');   end
if strcmpi(cfg.param{1},'index')
    if ~isfield(cfg,'sort'),      	cfg.sort        = 'ascend';             end
else
    if ~isfield(cfg,'sort'),      	cfg.sort        = 'descend';            end
    cfg = setcfg_movparam(cfg);
end

% function setcfg_movparam
%----------------------------------------
function cfg = setcfg_movparam(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'param'),       error('no movement parameters defined');    end
if ~iscell(cfg.param),              cfg.param       = {cfg.param};          end
if ~isfield(cfg,'expfun'),          cfg.expfun      = 'no';                 end
if ~isfalse(cfg.expfun) && ~strcmpi(cfg.expfun,'exp') && exist(cfg.expfun,'file') ~= 2 && isempty(regexp(cfg.expfun,'^km_movparamfun_','once'))
    cfg.expfun = ['km_movparamfun_' cfg.expfun];
end
if ~isfield(cfg,'clean'),       	cfg.clean       = 'yes';                end

% FIXME: this setcfg function is completely outdated. Please adapt on the
% basis of km_movparam. Good luck...

% FIXME: implement prmabs and prmdiff, see km_movparam

% isolate parameter names
prm = regexp(cfg.param,'^d?([a-zA-Z]+)','tokens');
prm = [prm{:}];         % concatenate parameters (and leave out empty ones)
prm = unique([prm{:}]); % find unique parameters
if isempty(prm), prm = cell(0); end;

% isolate phase indicators in parameter names
phs = regexp(cfg.param,'_([A-Z][a-zA-Z]*)','tokens');
phs = [phs{:}];         % concatenate parameters (and leave out empty ones)
phs = unique([phs{:}]); % find unique phases
if isempty(phs), phs = cell(0); end;

% isolate sensor indicators in parameter names
sns = regexp(cfg.param,'_([a-z]+)','tokens');
sns = [sns{:}];         % concatenate parameters (and leave out empty ones)
sns = unique([sns{:}]); % find unique sensors
if isempty(sns), sns = cell(0); end;

% isolate axis indicators in parameter names
% NOT IMPLEMENTED YET

% check cfg based on parameters
% TIME related
if any(ismember(prm,{'rt','mt','rmt'}))
	if ~isfield(cfg,'time'),        cfg.time        = struct;               end
	if ~isfield(cfg.time,'movnr'),  cfg.time.movnr  = 1;                    end
end
% VELOCITY related
if any(ismember(prm,{'v','mv','pv','tpv','rtpv','ppv','tppv','rtppv','apv','tapv','rtapv','npv','tnpv','rtnpv'}))
	if ~isfield(cfg,'vel'),     error('no velocity param specs provided');  end
end
if any(ismember(prm,{'tl','zv'}))
	if ~isfield(cfg,'dist'),     error('no distance param specs provided');  end
end
% ACCELERATION related
if any(ismember(prm,{'a','ma','pa','tpa','rtpa','pd','tpd','rtpd','nvc'}))
	if ~isfield(cfg,'acc'),     error('no acceleration param specs provided');  end
end
% GRIP related
if any(ismember(prm,{'mga'}))
    warning('KM:SetCFG:mga','Are you sure you intend to calculate the MGA [MEAN grip aperture]? Maybe you are looking for the PGA [PEAK grip aperture].');
end
if any(ismember(prm,{'tmga','rtmga'}))
    error('KM:SetCFG:mga','The TMGA and RTMGA parameters do not exist [(relative) time to MEAN grip aperture]. Maybe you are looking for the TPGA or RTPGA instead [(relative) time to PEAK grip aperture]?');
end
if any(ismember(prm,{'ga','mga','pga','tpga','rtpga'}))
    if ~isfield(cfg,'gripapt')
        if isfield(cfg,'gripapp')       % legacy app > apt
            cfg.gripapt = cfg.gripapp; cfg = rmfield(cfg,'gripapp');
        else
            error('no grip apt param specs provided')
        end
    end
end
if any(ismember(prm,{'gv','pgv','tpgv','rtpgv','npgv','tnpgv','rtnpgv'}))
    if ~isfield(cfg,'gripaptvel')
        if isfield(cfg,'gripappvel')    % legacy app > apt
            cfg.gripaptvel = cfg.gripappvel; cfg = rmfield(cfg,'gripappvel');
        else
            error('no grip apt vel param specs provided')
        end
    end
end
if any(ismember(prm,{'go'}))
	if ~isfield(cfg,'gripori'), error('no grip ori param specs provided');  end
end
% POSITION related
if any(ismember(prm,{'pos','maxpos','pmaxpos','amaxpos','nmaxpos','zpos','pzpos','nzpos'}))
	if ~isfield(cfg,'pos'),     error('no position param specs provided');  end
end

% check cfg based on phases
if any(ismember(phs,{'T','A'})) || any(ismember(prm,{'pa','tpa','rtpa','pd','tpd','rtpd'}))
	if ~isfield(cfg.vel,'T'),    error('no transport velocity param specs provided');	end
    if ~isfield(cfg.vel.T,'crit'),          cfg.vel.T.crit          = 0.3;               	end     % critical value
    if ~isfield(cfg.vel.T,'slide'),         cfg.vel.T.slide         = cfg.vel.T.crit/5;   	end     % let the on- and off-sets slide away till a local minimum or the slide criterium is reached
    if ~isfield(cfg.vel.T,'consecflg'),     cfg.vel.T.consecflg     = 'onset';          	end     % if a consecutive off- and on-set come together, move the 'onset', 'offset', of 'both' one sample away from the border
end
if any(ismember(phs,{'Tt','At'}))
	if ~isfield(cfg.vel,'Tt'),    error('no thumb transport velocity param specs provided');	end
    if ~isfield(cfg.vel.Tt,'crit'),     	cfg.vel.Tt.crit         = 0.3;               	end     % critical value
    if ~isfield(cfg.vel.Tt,'slide'),    	cfg.vel.Tt.slide        = cfg.vel.Tt.crit/5;   	end     % let the on- and off-sets slide away till a local minimum or the slide criterium is reached
    if ~isfield(cfg.vel.Tt,'consecflg'),	cfg.vel.Tt.consecflg	= 'onset';          	end     % if a consecutive off- and on-set come together, move the 'onset', 'offset', of 'both' one sample away from the border
end
if any(ismember(phs,{'Ti','Ai'}))
	if ~isfield(cfg.vel,'Ti'),    error('no index transport velocity param specs provided');	end
    if ~isfield(cfg.vel.Ti,'crit'),     	cfg.vel.Ti.crit         = 0.3;               	end     % critical value
    if ~isfield(cfg.vel.Ti,'slide'),    	cfg.vel.Ti.slide        = cfg.vel.Ti.crit/5;   	end     % let the on- and off-sets slide away till a local minimum or the slide criterium is reached
    if ~isfield(cfg.vel.Ti,'consecflg'),	cfg.vel.Ti.consecflg	= 'onset';          	end     % if a consecutive off- and on-set come together, move the 'onset', 'offset', of 'both' one sample away from the border
end

if any(ismember(phs,{'B','E','W','R','Rb','Re'}))
    if ~isfield(cfg,'offset'),      cfg.offset      = struct;               end
end
if any(ismember(phs,'B'))
	if ~isfield(cfg.offset,'B'), 	cfg.offset.B    = 0;                    end
end
if any(ismember(phs,'E'))
	if ~isfield(cfg.offset,'E'), 	cfg.offset.E    = 0;                    end
end
if any(ismember(phs,'Rb'))
	if ~isfield(cfg.offset,'Rb'), 	cfg.offset.Rb   = 0;                    end
end
if any(ismember(phs,'Re'))
	if ~isfield(cfg.offset,'Re'), 	cfg.offset.Re   = 0;                    end
end
% arbitrary reference
if any(ismember(phs,'R'))
    if ~isfield(cfg.offset,'R'),        error('a reference field must be supplied in cfg.movparam.offset.R');       end
    if ~isfield(cfg.offset.R,'offset'),	error('a window offset must be supplied in cfg.movparam.offset.R.offset');	end
    if ~isfield(cfg.offset.R,'ref'),    cfg.offset.R.ref  	= 'end';    	end
end
% arbitrary window
if any(ismember(phs,'W'))
    if ~isfield(cfg.offset,'win'),	error('a window field must be supplied in cfg.movparam.offset.W');              end
    if ~isfield(cfg.offset.W,'offset'),	error('a window offset must be supplied in cfg.movparam.offset.W.offset');	end
    if ~isfield(cfg.offset.W,'ref'),   	cfg.offset.W.ref	= 'beg';       	end
end

% legacy app > apt
if isfield(cfg,'gripapp')
    cfg.gripapt = cfg.gripapp;	cfg = rmfield(cfg,'gripapp');
end
if isfield(cfg,'gripappvel')
    cfg.gripaptvel = cfg.gripappvel;	cfg = rmfield(cfg,'gripappvel');
end
cfg.list = fieldnames(cfg);
cfg = setcfg_tseries(cfg);  % should this be: cfg.tseries = setcfg_tseries(cfg.tseries); ???

% function setcfg_tseries
%----------------------------------------
function cfg = setcfg_tseries(cfg)
if ~isfield(cfg,'proc'),            cfg.proc        = 'no';                 end
if ~isfield(cfg,'list'),            cfg.list        = {};                   end
if ~iscell(cfg.list),               cfg.list        = {cfg.list};           end
if ~isempty(cfg.list) && isfalse(cfg.list{1}),      cfg.list = {};          end
if isfalse(cfg.list),               return;                                 end
tseries = {'pos','vel','tvel','ttvel','itvel','acc','gripapt','gripori','gripaptvel','griporivel'};
cfg.list = tseries(ismember(tseries,cfg.list));
tseries = cfg.list;
for ts = 1:length(tseries)
    if ~isfield(cfg,tseries{ts}),	cfg.(tseries{ts})	= struct;         	end
    for i = 1:length(cfg.(tseries{ts}))
        if ~isfield(cfg.(tseries{ts}),'movname') || isempty(cfg.(tseries{ts})(i).movname)
            cfg.(tseries{ts})(i).movname = 'movement';
        end
        if ~isfield(cfg.(tseries{ts}),'movnr') || isempty(cfg.(tseries{ts})(i).movnr)
            cfg.(tseries{ts})(i).movnr = 1;
        end
        if ~isfield(cfg.(tseries{ts}),'marker') || isempty(cfg.(tseries{ts})(i).marker)
            error('plotting time-series %s: no marker',tseries{ts});
        end
        if ~isfield(cfg.(tseries{ts}),'axis') || isempty(cfg.(tseries{ts})(i).axis)
            error('plotting time-series %s: no axis',tseries{ts});
        end
        if ~iscell(cfg.(tseries{ts})(i).marker),	cfg.(tseries{ts})(i).marker	= {cfg.(tseries{ts})(i).marker};	end
        if ~iscell(cfg.(tseries{ts})(i).axis),      cfg.(tseries{ts})(i).axis	= {cfg.(tseries{ts})(i).axis};      end
        if isfield(cfg.(tseries{ts}),'trlrej') && ~iscell(cfg.(tseries{ts})(i).trlrej)
            cfg.(tseries{ts})(i).trlrej = {cfg.(tseries{ts})(i).trlrej};
        end
        if any(strcmpi(tseries{ts},{'pos','gripapt'}))
            if ~isfield(cfg.(tseries{ts}),'basecorr') || isempty(cfg.(tseries{ts})(i).basecorr)
                cfg.(tseries{ts})(i).basecorr = 'no';
            end
            if ~isfalse(cfg.(tseries{ts})(i).basecorr)
                if ~isstruct(cfg.(tseries{ts})(i).basecorr)
                    cfg.(tseries{ts})(i).basecorr = struct;
                end
                if ~isfield(cfg.(tseries{ts}).basecorr,'win') || isempty(cfg.(tseries{ts})(i).basecorr.win)
                    cfg.(tseries{ts})(i).basewin = 'no';
                end
                if ~isfield(cfg.(tseries{ts}).basecorr,'val') || isempty(cfg.(tseries{ts})(i).basecorr.val)
                    cfg.(tseries{ts})(i).basecorr.val = 0;
                end
                if ~isfield(cfg.(tseries{ts}).basecorr,'descrip') || isempty(cfg.(tseries{ts})(i).basecorr.descrip)
                    cfg.(tseries{ts})(i).basecorr.descrip = 'nanmedian';
                end
                if ~isfield(cfg.(tseries{ts}).basecorr,'corrflg') || isempty(cfg.(tseries{ts})(i).basecorr.corrflg)
                    cfg.(tseries{ts})(i).basecorr.corrflg = 'sess';
                end
            end
        end
    end
end

% function setcfg_trlrej
%----------------------------------------
function cfg = setcfg_trlrej(cfg)
if ~isfield(cfg,'trlrej') || istrue(cfg.trlrej,'free')
    cfg.trlrej = struct;
end
if isfalse(cfg.trlrej),             return;                                 end
cfg = setcfg_dirana(cfg);
if ~isfield(cfg.trlrej,'fun'),      cfg.trlrej.fun	= 'auto';               end
if ~strcmpi(cfg.trlrej.fun,'auto') && ~isfalse(cfg.trlrej.fun)
    if strcmpi(cfg.trlrej.fun,'exp')
        cfg = setcfg_exp(cfg);
        cfg.trlrej.fun = ['km_trlrejfun_' cfg.exp];
    elseif exist(cfg.trlrej.fun,'file') ~= 2 && isempty(regexp(cfg.trlrej.fun,'^km_trlrejfun_','once'))
        cfg.trlrej.fun = ['km_trlrejfun_' cfg.trlrej.fun];
    end
end
if ~isfield(cfg.trlrej,'apriori'),      cfg.trlrej.apriori      = struct;	end
if isfield(cfg.trlrej,'apriorifun')	% legacy
    if ~isfield(cfg.trlrej.apriori,'fun')
        cfg.trlrej.apriori.fun = cfg.trlrej.apriorifun;
    end
    cfg.trlrej = rmfield(cfg.trlrej,'apriorifun');
end
if ~isfield(cfg.trlrej.apriori,'fun'),	cfg.trlrej.apriori.fun	= 'no';     end
if ~isfalse(cfg.trlrej.apriori.fun)
    if istrue(cfg.trlrej.apriori.fun,'free'),cfg.trlrej.apriori.fun	= 'exp';end
    if strcmpi(cfg.trlrej.apriori.fun,'exp')
        cfg = setcfg_exp(cfg);
        cfg.trlrej.apriori.fun = ['km_aprioritrlrejfun_' cfg.exp];
    elseif exist(cfg.trlrej.apriori.fun,'file') ~= 2 && isempty(regexp(cfg.trlrej.apriori.fun,'^km_aprioritrlrejfun_','once'))
        cfg.trlrej.apriori.fun = ['km_aprioritrlrejfun_' cfg.trlrej.apriori.fun];
    end
end
if isfield(cfg.trlrej,'mode')   % legacy
    cfg.trlrej.param	= cfg.trlrej.mode;
    cfg.trlrej          = rmfield(cfg.trlrej,'mode');
end
if ~isfield(cfg.trlrej,'param')
    if strcmpi(cfg.trlrej.fun,'auto')
        cfg.trlrej.param  	= {'apriori','rt','mt','td','mv','pv','mga','ega','dego'};
    else
        cfg.trlrej.param    = {};
    end
end
if ~iscell(cfg.trlrej.param),     	cfg.trlrej.param    = {cfg.trlrej.param};           end
if ~isfield(cfg.trlrej.apriori,'param'),cfg.trlrej.apriori.param    = cfg.trlrej.param;	end
if ~iscell(cfg.trlrej.apriori.param),   cfg.trlrej.apriori.param    = {cfg.trlrej.apriori.param};	end
cfg.trlrej.apriori.param = cfg.trlrej.apriori.param(~strcmpi(cfg.trlrej.apriori.param,'apriori'));
for p = 1:length(cfg.trlrej.apriori.param)
    pname = cfg.trlrej.apriori.param{p};
    if ~isfield(cfg.trlrej,pname)
        warning('KM:SetCFG:TrlRejAprioriParam','No cutoffs are specified for the apriori parameter [%s]',pname);
        cfg.trlrej.(pname) = [-Inf Inf];
    end
end
if ~isfield(cfg.trlrej.apriori,'cutmode'),	cfg.trlrej.apriori.cutmode	= 'higher';  	end
if ~isfield(cfg.trlrej,'check'),	cfg.trlrej.check	= struct;       	end
if isfield(cfg.trlrej,'reportmode')     % legacy
    cfg.trlrej.reportparam = cfg.trlrej.reportmode;
    cfg.trlrej = rmfield(cfg.trlrej,'reportmode');
end
if ~isfield(cfg.trlrej,'reportparam')
    if strcmpi(cfg.trlrej.fun,'auto')
        cfg.trlrej.reportparam   = cfg.trlrej.param;
    else
        cfg.trlrej.reportparam   = {'auto'};
    end
end
if ~iscell(cfg.trlrej.reportparam),	cfg.trlrej.reportparam = {cfg.trlrej.reportparam};	end

% function setcfg_movparamCI
%----------------------------------------
function cfg = setcfg_movparamCI(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'trlrej'),          cfg.trlrej      = 'no';                 end
cfg.trlrej = setcfg_idxr(cfg.trlrej);
cfg = setcfg_factlvl(cfg);

% function setcfg_eposerr
%----------------------------------------
function cfg = setcfg_eposerr(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'trlrej'),          cfg.trlrej      = 'no';                 end
cfg.trlrej = setcfg_idxr(cfg.trlrej);
if ~isfield(cfg,'movname'),         cfg.movnr       = 'movement';           end
if ~isfield(cfg,'movnr'),           cfg.movnr       = 1;                    end
if ~isfield(cfg,'marker'),          error('no marker specified');           end
if ~iscell(cfg.marker),             cfg.marker      = {cfg.marker};         end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end
if ~isfield(cfg,'endoffset'),       cfg.endoffset   = 0;                    end
cfg = setcfg_factlvl(cfg);
if ~isfield(cfg,'expfun'),          cfg.expfun      = 'no';                 end
if ~isfalse(cfg.expfun) && exist(cfg.expfun,'file') ~= 2 && isempty(regexp(cfg.expfun,'^km_eposerrfun_','once'))
	cfg.expfun = ['km_eposerrfun_' cfg.expfun];
end

% function setcfg_traject
%----------------------------------------
function cfg = setcfg_traject(cfg)
if isempty(cfg),                    cfg             = 'no';                 end
if isfalse(cfg),                    return;                                 end
if ~isfield(cfg,'type'),            cfg.type        = {'ellipsoid','project','intersect'};	end
if ~iscell(cfg.type),               cfg.type        = {cfg.type};           end
cfg.type = km_labelselection(cfg.type,{'ellipsoid','project','intersect'});
if isempty(cfg.type),   error('No (supported) type of trajectory specification is selected, so there is nothing to do.');   end
if ~isfield(cfg,'npnts_rept'),      cfg.npnts_rept  = 1000;                 end     % number of points for each trial trajectory
if ~isfield(cfg,'npnts_avg'),       cfg.npnts_avg   = 300;                  end     % number of points for the average trajectory
if ~isfield(cfg,'npnts_circle'),    cfg.npnts_circle= 51;                   end     % number of points around a circle (used for plotting)
if ~isfield(cfg,'calc_pnts'),       cfg.calc_pnts    = 'yes';               end
if ~isfield(cfg,'avgcalc_iter'),    cfg.avgcalc_iter= 1;                    end
if ~isfield(cfg,'avgcalc_outlier'), cfg.avgcalc_outlier = 'best';           end
if ~isfield(cfg,'nsesscomb'),       cfg.nsesscomb   = 0;                    end
if ~isfield(cfg,'sesscounter'),     cfg.sesscounter = 1;                    end
if ~isfield(cfg,'trlrej'),          cfg.trlrej      = 'no';                 end
cfg.trlrej = setcfg_idxr(cfg.trlrej);
if ~isfield(cfg,'movname'),         cfg.movnr       = 'movement';           end
if ~isfield(cfg,'movnr'),           cfg.movnr       = 1;                    end
if ~isfield(cfg,'marker'),          error('no marker specified');           end
if ~iscell(cfg.marker),             cfg.marker      = {cfg.marker};         end
if ~isfield(cfg,'axis'),            cfg.axis        = 'xyz';                end
cfg = setcfg_factlvl(cfg);
if ~isfield(cfg,'ntrl'),            cfg.ntrl.default= [4 Inf];              end
for f = 1:length(cfg.fname)
    if ~isfield(cfg.ntrl,cfg.fname{f}), cfg.ntrl.(cfg.fname{f}) = cfg.ntrl.default; end
end
if ~isfield(cfg,'expfun'),          cfg.expfun      = 'no';                 end
if ~isfalse(cfg.expfun) && exist(cfg.expfun,'file') ~= 2 && isempty(regexp(cfg.expfun,'^km_trajectfun_','once'))
	cfg.expfun = ['km_trajectfun_' cfg.expfun];
end

% function setcfg_addvar
%----------------------------------------
function cfg = setcfg_addvar(cfg)
if ~isfield(cfg,'addvar')
    cfg.addvar = struct('newvar',{},'oldvar',{},'fun',{});
end
if isempty(cfg.addvar),             return;                                 end
for i = 1:length(cfg.addvar)
    if ~isfield(cfg.addvar(i),'newvar') || isempty(cfg.addvar(i).newvar) || ~ischar(cfg.addvar(i).newvar)
        error('no (correct) name for the new variable give in xxx.addvar(%i).newvar',i);
    end
    if ~isfield(cfg.addvar(i),'oldvar') || isempty(cfg.addvar(i).oldvar)
        cfg.addvar(i).oldvar = {''};
    end
    if ~iscell(cfg.addvar(i).oldvar)
        cfg.addvar(i).oldvar = {cfg.addvar(i).oldvar};
    end
    if ~all(cellfun(@ischar,cfg.addvar(i).oldvar))
        error('all cell fields in xxx.addvar(%i).oldvar should be character arrays',i);
    end
    if ~isfield(cfg.addvar(i),'fun') || isempty(cfg.addvar(i).fun) || ~ischar(cfg.addvar(i).fun)
        error('an evaluable function should be entered in xxx.addvar(%i).fun',i);
    end
end
    
% function setcfg_factlvl
%----------------------------------------
function cfg = setcfg_factlvl(cfg)
if ~isfield(cfg,'fact') || isempty(cfg.fact)
    error('you should specify factors in xxx.fact');
end
if ~iscell(cfg.fact),               cfg.fact = {cfg.fact};                  end
if ~isfield(cfg,'lvl') || isempty(cfg.lvl)
    error('you should specify levels in xxx.lvl');
end
% select factor names
cfg.fact	= km_parsefactstr(cfg.fact);
cfg.nfact	= length(cfg.fact);
cfg.fname   = cellfun(@(x) sprintf('%sx',x{:}),cfg.fact,'UniformOutput',false);
cfg.fname   = cellfun(@(x) x(1:end-1),cfg.fname,'UniformOutput',false);
cfg.bfact	= unique([cfg.fact{:}]);
cfg.nbfact	= length(cfg.bfact);
flds = fieldnames(cfg.lvl);
if ~all(ismember(cfg.bfact,flds))
    error('not all factors have denoted levels in xxx.lvl');
end

% function setcfg_contrfactrem
%----------------------------------------
function cfg = setcfg_contrfactrem(cfg,flg)
if ~isfield(cfg,flg) || isfalse(cfg.(flg))
    cfg.(flg) = {};
end
if ~iscell(cfg.(flg)),              cfg.(flg)       = {cfg.(flg)};          end
if isempty(cfg.(flg)) || ~iscell(cfg.(flg){1}),	cfg.(flg)	= {cfg.(flg)};  end
if length(cfg.(flg)) == 1 && length(cfg.fact) > 1
    cfg.(flg) = repmat(cfg.(flg),1,length(cfg.fact));
elseif length(cfg.fact) > 1
    error('the number of factors in cfg.%s does not match those in cfg.fact',flg);
end
for f = 1:length(cfg.fact)
    if strcmpi(flg,'contr')
        for c = 1:length(cfg.(flg){f})
            if length(cfg.lvl.(cfg.(flg){f}{c})) ~= 2
                warning('KMContr:badnrlvls','the factor [%s] has %d levels while contrasts can only be calculated for 2 levels',cfg.(flg){f}{c},length(cfg.lvl.(cfg.(flg){f}{c})));
            end
        end
        cfg.(flg){f} = cfg.(flg){f}(cellfun(@(str) length(cfg.lvl.(str)),cfg.(flg){f})==2);
    end
    if ~all(ismember(cfg.(flg){f},cfg.fact{f}))
        warning('KMContr:nomatch','some elements in cfg.%s do not match those in cfg.fact',flg);
        cfg.(flg){f} = cfg.(flg){f}(ismember(cfg.(flg){f},cfg.fact{f}));
    end
end

%% PLOTTING OF RESULTS
%--------------------------------------------------------------------------

% function setcfg_plot
%----------------------------------------
function cfg = setcfg_plot(cfg)
cfg = setcfg_basic(cfg);
cfg = setcfg_subj(cfg);
cfg = setcfg_sess(cfg);
if ~iscell(cfg.sess),               cfg.sess        = {cfg.sess};           end
cfg = setcfg_dirroot(cfg);
if ~isfield(cfg,'load'),            cfg.load        = 'no';                 end
if ~isfield(cfg,'plot')
    error('you must specify configuration settings in cfg.plot');
end
cfg.plot = setcfg_plot_plot(cfg.plot);
cfg.tseries.list = cfg.plot.tseries;
cfg.tseries = setcfg_tseries(cfg.tseries);
if length(cfg.subj) > 1
    cfg.plot.subjstr = 'Group';
else
    cfg.plot.subjstr = cfg.subj{1};
end
cfg.plot.sessstr = strcat(cfg.sess{:});
% remove low dashes
cfg.plot.sesstr = strrep(cfg.plot.sessstr,'_','');

% function setcfg_plot_plot
%----------------------------------------
function cfg = setcfg_plot_plot(cfg)
if ~isfield(cfg,'fetchdatafromcfg'),cfg.fetchdatafromcfg	= 'no';         end
if ~isfield(cfg,'reptcontrib'),     cfg.reptcontrib	= 6;                    end
if ~isfield(cfg,'subjcontrib'),     cfg.subjcontrib	= 20;                   end
if ~isfield(cfg,'pimp'),            cfg.pimp        = 'yes';                end
if ~isfield(cfg,'export'),          cfg.export      = {'-pdf'};             end
if ~iscell(cfg.export),             cfg.export      = {cfg.export};         end
for i = 1:length(cfg.export)
    if ~strcmpi(cfg.export{i}(1),'-')
        cfg.export{i} = ['-' cfg.export{i}];
    end
end
if ~isfield(cfg,'trlrej'),          cfg.trlrej      = 'no';                 end
cfg.trlrej = setcfg_idxr(cfg.trlrej);
if ~isfield(cfg,'param'),           cfg.param       = {};                   end
if isfalse(cfg.param),              cfg.param       = {};                   end
if ~iscell(cfg.param),              cfg.param       = {cfg.param};          end
if ~isfield(cfg.trlrej,'pairwise'), cfg.trlrej.pairwise = cfg.param;        end
% set tseries configuration
if ~isfield(cfg,'tseries'),         cfg.tseries     = {};                   end
if isfalse(cfg.tseries),            cfg.tseries     = {};                   end
% if isfalse(cfg.param) && isfalse(cfg.tseries)
%     error('either a kinematic parameter or a time series should be selected.');
% else
if ~isfalse(cfg.param) && ~isfalse(cfg.tseries)
    error('for now parameters and time series can not be plotted simultaneously.');
end
% confidence interval
idx_CI = ~cellfun(@isempty,regexp(cfg.param,'CI$'));
if any(idx_CI)
    if sum(idx_CI) == 1 && strcmpi(cfg.param{idx_CI},'ci')
        if ~isfield(cfg.CI,'param') || ~iscell(cfg.CI.param) || isempty(cfg.CI.param)
            error('confidence interval parameters are requested to be plotted but not specified.');
        end
    else
        paramCI = regexprep(cfg.param(idx_CI),'CI$','');
        if ~isstruct(cfg.CI), cfg.CI = struct;  end
        if isfield(cfg.CI,'param') && iscell(cfg.CI.param) && ~isempty(cfg.CI.param)
            cfg.CI.param = unique([cfg.CI.param(:); paramCI(:)]);
        else
            cfg.CI.param = paramCI;
        end
    end
    cfg.param = cfg.param(~idx_CI);
else
%elseif isfalse(cfg.CI) || ~isfield(cfg.CI,'param') || isfalse(cfg.CI.param)
    cfg.CI = false;
end

% end-point-error
idx_epe = ismember(cfg.param,'epe');
cfg.param = cfg.param(~idx_epe);
if any(idx_epe)
    if ~isstruct(cfg.epe),          cfg.epe         = struct();             end
    if ~isfield(cfg.epe,'marker'),  cfg.epe.marker  = {'index','thumb'};    end
    if ~iscell(cfg.epe.marker),     cfg.epe.marker  = {cfg.epe.marker};     end
    if ~isfield(cfg.epe,'axis'),    cfg.epe.axis    = 'xyz';                end
    if ~isfield(cfg.epe,'param'),   cfg.epe.param   = 'a';                  end
    if ~isfield(cfg.epe,'scale')
        if strcmpi(cfg.epe.param,'volume')
            cfg.epe.scale = 100^3;  % m³ to cm³
        elseif strcmpi(cfg.epe.param,'area')
            cfg.epe.scale = 100^2;  % m² to cm²
        else
            cfg.epe.scale = 100;    % m to cm
        end
    end
    % % select factors and levels or conditions
    % cfg = setcfg_factlvl(cfg);
else
    cfg.epe = false;
end

% factors and levels to exclude
if ~isfield(cfg,'factexcl'),        cfg.factexcl    = {};                   end
if ~iscell(cfg.factexcl),           cfg.factexcl    = {cfg.factexcl};       end
if ~isfield(cfg,'lvlexcl'),         cfg.lvlexcl     = struct;               end
flds = fieldnames(cfg.lvlexcl);
if ~all(ismember(cfg.factexcl,flds))
    error('not all excluding factors have denoted levels in xxx.lvlexcl');
end

% subject and group descriptives
if ~isfield(cfg,'descrip'),         cfg.descrip     = struct;               end
if ~isstruct(cfg.descrip)
    tmpdescrip = cfg.descrip;
    cfg.descrip = struct;
    cfg.descrip.group	= tmpdescrip;
    cfg.descrip.subj    = tmpdescrip;
end
if ~isfield(cfg.descrip,'group'),   cfg.descrip.group   = 'mean';           end
if ~isfield(cfg.descrip,'subj'),    cfg.descrip.subj    = 'mean';           end
if ~isfield(cfg.descrip,'err'),     cfg.descrip.err     = 'sem';            end
if ~isfield(cfg.descrip,'contrerr'),cfg.descrip.contrerr= 'no';             end
if ~ismember(cfg.descrip.group,{'mean','median'})
    error('unknown group descriptive (%s) for parameter plotting',cfg.descrip.group);
end
if ~ismember(cfg.descrip.subj,{'mean','median','logmean','logmedian','meanlog','medianlog','var'})
    error('unknown descriptive (%s) for parameter plotting',cfg.descrip.subj);
end
if ~ismember(cfg.descrip.err,{'var','sem'})
    error('unknown error descriptive (%s) for parameter plotting',cfg.descrip.err);
end

% add variables to trl matrix if requested
cfg = setcfg_addvar(cfg);

% select factors and levels or conditions
cfg = setcfg_factlvl(cfg);

% select contrasts based on factors
cfg = setcfg_contrfactrem(cfg,'contr');
cfg = setcfg_contrfactrem(cfg,'factrem');

% set plot specifications (fill/line colors, style, width, markers)
cfg = setcfg_plot_plotspec(cfg);

% The code below is used when time-series are defined using conditions and
% contrasts
if ~isfalse(cfg.tseries)
    cfg = setcfg_factlvl(cfg);
    if ~isfield(cfg,'cond')||isempty(cfg.cond),	cfg.cond	= {'all'};   	end
    % select condition names, combinations and contrasts
    cfg.ncond   = length(cfg.cond);
    cfg.cond	= parsecondstr(cfg.cond);
    % determine the basis (parts) of all conditions
    cfg.bcond	= getbasiscond(cfg.cond);
    cfg.nbcond  = length(cfg.bcond);
    cfg.bcname	= {cfg.bcond(:).name};
    % restrict length of condition names
    cfg.cname   = spellout({cfg.cond(:).name});
    for c = 1:length(cfg.cname)
        if length(cfg.cname{c}) > 40
            cfg.cname{c} = sprintf('cond%d',c);
        end
    end
%     % color and line specifications
%     if length(cfg.cond)==1
%         flcl	= 'k';
%         lnstyle	= {'-'};
%     elseif length(cfg.cond)==2
%         flcl	= 'rg';
%         lnstyle	= {'-','-'};
%     elseif length(cfg.cond)==3
%         flcl	= 'rgb';
%         lnstyle	= {'-','-','-'};
%     elseif length(cfg.cond)==4
%         flcl	= 'rrgg';
%         lnstyle	= {'-','--','-','--'};
%         % flcl	= 'bmrg';
%         % lnstyle	= {'-','-','-','-'};
%     elseif length(cfg.cond)<8
%         flcl	= 'rgbcmyk';
%         lnstyle	= {'-','-','-','-','-','-','-'};
%     else
%         flcl	= 'rgbcmykrgbcmykrgbcmyk';          % color specifications of line plots
%         lnstyle	= repmat({'-'},1,length(flcl));   % line style specifications of line plots
%     end
%     if ~isfield(cfg,'flcl')
%         cfg.flcl	= flcl;
%     end
%     if ~isfield(cfg,'lnstyle')
%         cfg.lnstyle	= lnstyle;
%     end

% draw sensors (and connections) at defined intervals or positions
if ~isfield(cfg,'drawsens'),        cfg.drawsens = [];                      end
if ~isfalse(cfg.drawsens)
    if ~isstruct(cfg.drawsens),     cfg.drawsens = struct;                  end
    for ts = 1:length(cfg.tseries)
        if length(cfg.drawsens) < ts,   cfg.drawsens(ts).marker = [];       end
        if ~isfield(cfg.drawsens(ts),'marker') || isempty(cfg.drawsens(ts).marker)
            cfg.drawsens(ts).marker = {'all'};
        end
        if ~iscell(cfg.drawsens(ts).marker),	cfg.drawsens(ts).marker = {cfg.drawsens(ts).marker};	end
        if ~isfield(cfg.drawsens(ts),'tpoint')
            cfg.drawsens(ts).tpoint = [];
        end
        if ~isfield(cfg.drawsens(ts),'nmark')
            cfg.drawsens(ts).tpoint = [];
        end
        if isfalse(cfg.drawsens(ts).tpoint)
            cfg.drawsens(ts).tpoint = [];
            if isfalse(cfg.drawsens(ts).nmark)
                cfg.drawsens(ts).nmark = 4;
            end
        end
        if ~isfield(cfg.drawsens(ts),'connect') || isempty(cfg.drawsens(ts).connect)
            cfg.drawsens(ts).connect = 'yes';
        end
    end
end

cfg.plot.drawsens.nmark     = 3;
cfg.plot.drawsens.tpoint    = 66;
cfg.plot.drawsens.connect	= 'yes';

end

% function setcfg_plot_plotspec
%----------------------------------------
function cfg = setcfg_plot_plotspec(cfg)

% default line and marker specifications
fflnstyle = {'-','--',':','-.','-','--',':','-.'};
ffmarker = {'o','*','x ','+','.','d','s','p','h'};

% temporary store color, line and marker specifications before
% constructing defaults
tmp = [];
if isfield(cfg,'flcl')
    tmp.flcl	= cfg.flcl;
    cfg         = rmfield(cfg,'flcl');
end
if isfield(cfg,'lncl')
    tmp.lncl	= cfg.lncl;
    cfg         = rmfield(cfg,'lncl');
    if ~isfield(tmp,'flcl')
        tmp.flcl= tmp.lncl;
    end
end
if isfield(cfg,'lnstyle')
    if ~iscell(cfg.lnstyle), cfg.lnstyle = num2cell(cfg.lnstyle); end
    if ~iscell(cfg.lnstyle{1}), cfg.lnstyle = {cfg.lnstyle}; end
    tmp.lnstyle	= cfg.lnstyle;
    cfg         = rmfield(cfg,'lnstyle');
end
if isfield(cfg,'lnwidth')
    tmp.lnwidth	= cfg.lnwidth;
    cfg         = rmfield(cfg,'lnwidth');
end
if isfield(cfg,'marker')
    tmp.marker	= cfg.marker;
    cfg         = rmfield(cfg,'marker');
end

% loop over factors
for f = 1:cfg.nfact
    tmpfact = cfg.fact{f}(~ismember(cfg.fact{f},cfg.contr{f}));
    nff = length(tmpfact);
    if nff > 4
        error('You requested %d factors, but maximally 4 factors are supported in the parameter plotting. You could use cfg.plot.factexcl to select specific factors.',nff);
    end
    
    % get levels
    lvl = cell(1,nff);
    for ff = 1:nff
        [~, lvl{ff}] = unique(cfg.lvl.(tmpfact{ff})(:)');
    end
    tmplvl = arraycomb(lvl{:});
    lvl = ones(size(tmplvl,1),4);
    lvl(:,1:size(tmplvl,2)) = tmplvl;
    
    % base flcl on the number of levels of the first factor
    if ~isempty(tmpfact)
        nflcl = length(cfg.lvl.(tmpfact{1}));
    else
        nflcl = 1;
    end
    switch nflcl
        case 1,     flcl = 'h';
        case {2,3}, flcl = 'rgb';
        case 4,     flcl = 'borg';
        otherwise,  flcl = 'rgbcmyh'; flcl = repmat(flcl,1,3);
    end
    flcl = str2rgb(flcl);
    
    % base lncl on the number of levels of the second factor
    if length(tmpfact) > 1
        nlncl = length(cfg.lvl.(tmpfact{2}));
    else
        nlncl = 1;
    end
    switch nlncl
        case 1,     lncl = 'k';
        case {2,3}, lncl = 'bom';
        case 4,     lncl = 'cmyg';
        otherwise,  lncl = 'bomrgck'; lncl = repmat(lncl,1,3);
    end
    lncl = str2rgb(lncl);
    
    % assign ffclspec, fflnstyle, ffmarker to all level combinations
    if ~isempty(lvl)
        flcl	= flcl(lvl(:,1),:);
        lncl	= lncl(lvl(:,2),:);
        lnstyle	= fflnstyle(lvl(:,2));
        marker	= ffmarker(lvl(:,3));
        lnwidth	= repmat(2,1,size(flcl,1));
    else
        flcl	= flcl(1,:);
        lncl	= lncl(1,:);
        lnstyle	= fflnstyle(1);
        marker	= ffmarker(1);
        lnwidth	= repmat(2,1,size(flcl,1));
    end
    
    % use defaults or user specified specs
    specnames = {'flcl','lncl','lnstyle','lnwidth','marker'};
    for sn = 1:length(specnames)
        sname = specnames{sn};
        if isfield(tmp,sname)
            if iscell(tmp.(sname))
                cfg.(sname){f}	= tmp.(sname){f};
            else
                cfg.(sname){f}	= tmp.(sname);
            end
        else
            cfg.(sname){f}	= eval(sname);
        end
    end
end

% function setcfg_paramreport
%----------------------------------------
function cfg = setcfg_paramreport(cfg)
if isfalse(cfg),                    cfg             = 'no';              	end
if isfalse(cfg),                    return;                                 end
if ~isstruct(cfg),                  cfg             = struct;              	end
if ~isfield(cfg,'write'),           cfg.write       = {'grouprepmeas'};     end
if ~iscell(cfg.write),              cfg.write       = {cfg.write};          end
if any(ismember(cfg.write,'group'))
    cfg.write{ismember(cfg.write,'group')} = 'grouprepmeas';
end
if any(ismember(cfg.write,{'subj','groupuni'}))
    if ~isfield(cfg,'param'),       cfg.param       = {'all','-trlbeg','-trlend','trloff','-t_*'};              end
    if ~iscell(cfg.param),          cfg.param       = {cfg.param};          end
                                    cfg.param       = unique(cfg.param);
    if ~isfield(cfg,'reject'),      cfg.reject      = {'all'};              end
    if ~iscell(cfg.reject),         cfg.reject      = {cfg.reject};         end
                                    cfg.reject      = unique(cfg.reject);
    if any(ismember(cfg.reject,'comb'))
        cfg.reject{ismember(cfg.reject,'comb')}     = 'reject';
    end
    if ~isfield(cfg,'exclude'),     cfg.exclude     = 'no';                 end
    if istrue(cfg.exclude,'free'),  cfg.exclude     = 'cases';              end
else
    %if isfield(cfg,'reject'),       cfg             = rmfield(cfg,'reject');end
    %if isfield(cfg,'exclude'),      cfg             = rmfield(cfg,'exclude');end
end

% function setcfg_paramplot
%----------------------------------------
function cfg = setcfg_paramplot(cfg)
if isfalse(cfg),                    cfg             = 'no';              	end
if isfalse(cfg),                    return;                                 end
if istrue(cfg,'free'),              cfg             = struct;              	end
if ~isfield(cfg,'type'),            cfg.type            = {'errorbar'};     end
if ~iscell(cfg.type),               cfg.type            = {cfg.type};       end
%if isfalse(cfg.descrip.err)
%    for f = 1:length(cfg.type)
%        if strcmpi(cfg.type{f},'errorbar')
%            cfg.type{f} = 'line';
%        end
%    end
%end

% function setcfg_tseriesplot
%----------------------------------------
function cfg = setcfg_tseriesplot(cfg)
if isfalse(cfg),                    cfg             = 'no';              	end
if isfalse(cfg),                    return;                                 end
if istrue(cfg,'free'),              cfg             = struct;              	end


%% FILE INPUT AND OUTPUT
%--------------------------------------------------------------------------

% function setcfg_subj
%----------------------------------------
function cfg = setcfg_subj(cfg)
if ~isfield(cfg,'subj'),        error('no subject specified');              end

% function setcfg_sess
%----------------------------------------
function cfg = setcfg_sess(cfg)
%if ~isfield(cfg,'sess'),       error('no session specified');              end
if ~isfield(cfg,'sess'),        cfg.sess = '';                              end
if iscell(cfg.sess)
    for ss = 1:length(cfg.sess)
        if isempty(cfg.sess{ss}),   cfg.sess{ss} = '';  end
        if ~isempty(cfg.sess{ss}) && ~strcmpi(cfg.sess{ss}(1),'_')
            cfg.sess{ss} = ['_' cfg.sess{ss}];
        end
    end
elseif ~isempty(cfg.sess) && ~strcmpi(cfg.sess(1),'_')
    cfg.sess = ['_' cfg.sess];
end
%% FIXME: does cfg.sess need to be a cell array? for km_plot it does, for
%% km_save it needs to be a string...
    
% function setcfg_dir
%----------------------------------------
function cfg = setcfg_dir(cfg)
if ~isfield(cfg,'dir'),         error('no directory specified');            end

% function setcfg_dirroot
%----------------------------------------
function cfg = setcfg_dirroot(cfg)
cfg = setcfg_dir(cfg);
if ~isfield(cfg.dir,'root'),	error('no root directory: cfg.dir.root');   end
if ~exist(cfg.dir.root,'dir'),  error('root directory does not exist');     end

% function setcfg_dirdata
%----------------------------------------
function cfg = setcfg_dirdata(cfg)
cfg = setcfg_dir(cfg);
if ~isfield(cfg.dir,'data'),	error('no data directory specified');       end
if ~exist(cfg.dir.data,'dir'),  error('data directory does not exist');     end

% function setcfg_dirlog
%----------------------------------------
function cfg = setcfg_dirlog(cfg)
cfg = setcfg_dir(cfg);
if ~isfield(cfg.dir,'log'),     error('no logfile directory specified');    end
if ~exist(cfg.dir.log,'dir'),   error('logfile directory does not exist');  end

% function setcfg_dirpreproc
%----------------------------------------
function cfg = setcfg_dirpreproc(cfg)
cfg = setcfg_dir(cfg);
if ~isfield(cfg.dir,'preproc'),	error('no preproc directory specified');	end
%if ~exist(cfg.dir.preproc,'dir'),   mkdir(cfg.dir.preproc);                 end
%if ~exist(cfg.dir.preproc,'dir'),error('preproc directory does not exist');end

% function setcfg_dirana
%----------------------------------------
function cfg = setcfg_dirana(cfg)
cfg = setcfg_dir(cfg);
if ~isfield(cfg.dir,'ana'),     error('no analysis dir: cfg.dir.ana');      end
%if ~exist(cfg.dir.ana,'dir'),       mkdir(cfg.dir.ana);                     end

% function setcfg_dirproc
%----------------------------------------
function cfg = setcfg_dirproc(cfg)
try 
    if ~isempty(regexpi(cfg.dataset,'_ANA(_IDX[RC])?$','once'))
        cfg = setcfg_dirana(cfg);
        cfg.dir.proc = cfg.dir.ana;
    else
        cfg = setcfg_dirpreproc(cfg);
        cfg.dir.proc = cfg.dir.preproc;
    end
end

% function setcfg_dirreport
%----------------------------------------
function cfg = setcfg_dirreport(cfg)
cfg = setcfg_dir(cfg);
if ~isfield(cfg.dir,'report'),      error('no report directory: cfg.dir.report');   end
%if ~exist(cfg.dir.report,'dir'),	mkdir(cfg.dir.report);                          end
