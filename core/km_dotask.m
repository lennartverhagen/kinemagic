function varargout = km_dotask(cfg,data)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------


%% Set input and configuration
%----------------------------------------
% set input
if nargin < 1,  cfg  = [];	end
if nargin < 2,  data = [];	end
    
% set basic configuration
cfg = km_setcfg(cfg,'basic');

% let's see what there is to do
if isfield(cfg,'task')
    task = cfg.task;
else
    task = {};
end

% get the correct naming conventions for the different tasks
task = km_settask(task,'strict');

% split tasks in chunks if repetitions are found
[cfg,data,rtrnflg] = km_chunktask(cfg,data);
if rtrnflg,   return; end

% check the configuration and try to fill in the defaults where needed
cfg = km_setcfg(cfg,task);

% set output arguments ready if the function returns early
varargout = {};
if nargout == 1
    varargout = {cfg};
elseif nargout == 2
    varargout = {cfg,data};
end

% select only the tasks where there is actually something to do
if ~iscell(task),   task = {task};  end
check = true(1,length(task));
for t = 1:length(task)
    if isfield(cfg,task{t}) && isfalse(cfg.(task{t}))
        check(t) = false;
    end
end
task = task(check);

% return if there is nothing to do
if isempty(task)
    warning('KM:Return','Nothing to do...');
    return
end


%% Get data in shape for processing
%----------------------------------------
% read in data fresh if requested
if strcmpi(task{1},'read')
    [cfg,data] = km_read(cfg);
    task = task(2:end);
end

% check to see if CFG and DATA need to be loaded
if isempty(data)
    flg_load_cfg = true;
    flg_load_data = true;
else
    flg_load_cfg = false;
    flg_load_data = false;
end

% check to see if only logfile can be read without loading data
if isequal(task,{'trlrej'})
    flg_load_data = false;
    if isfield(cfg,'trl')
        flg_load_cfg = false;
    else
        flg_load_cfg = true;
    end
elseif all(ismember(task,{'readlog','trldef','trlrej'}))
    flg_load_cfg = false;
    flg_load_data = false;
end

% initialize data structure if not loaded
if isempty(data) && ~flg_load_data
    data.dimord = 'marker_axis_time';
end
    
% check configuration if data needs to be loaded and/or saved
if flg_load_data || flg_load_cfg
    cfg = km_setcfg(cfg,{'subj','sess','dirpreproc','dirproc'});
end

% load configuration if requested
if flg_load_cfg
    fname_cfg = fullfile(cfg.dir.proc,sprintf('%s%s%s_CFG',cfg.subj,cfg.sess,cfg.dataset));
    if exist([fname_cfg '.mat'],'file')
        datacfg = load(fname_cfg);
        datacfg = datacfg.cfg;              % bring the loaded configuration structure to the right level
        cfg = combstruct(datacfg,cfg);      % combine the data and input configuration structures
        cfg.dataset = datacfg.dataset;      % but always overwrite the dataset type of the input cfg with the one from the data
    end
end

% load data if requested
if flg_load_data
    fname_data = fullfile(cfg.dir.proc,sprintf('%s%s%s_DATA',cfg.subj,cfg.sess,cfg.dataset));
    load(fname_data);
    % % update dataset name
    % if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
end

% set desired dimension order
data = km_setdimord(data,data.dimord,'marker_axis_time');


%% Execute task
%----------------------------------------
% perform tasks sequentially
for t = 1:length(task)
    
    % define function name
    funname = ['km_' task{t}];
    
    % organize input
    varargin = {};
    if nargin(funname) == 1
        varargin = {cfg};
    elseif nargin(funname) == 2
        varargin = {cfg,data};
    end
    
    % evaluate appropriate processing function
    if nargout(funname) == 0
        feval(funname,varargin{:});
    elseif nargout(funname) == 1
        cfg = feval(funname,varargin{:});
    elseif nargout(funname) == 2 || nargout(funname) == -1
        [cfg,data] = feval(funname,varargin{:});
    end
    
end


%% Set output
%----------------------------------------
% set dimension order back
data = km_setdimord(data,'marker_axis_time',data.dimord);

% save configuration and data
km_save(cfg,data);   

% return output
varargout = {cfg,data};
if nargout == 0 && ~istrue(cfg.save,'free') && ~isequal(cfg.save,'cfg') && isequal(cfg.save,'data')
    warning('KM:NoOutput','No output was requested to be returned or saved');
end
