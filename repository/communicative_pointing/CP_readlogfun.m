function log = CP_readlogfun(cfg)
%--------------------------------------------------------------------------
% This function reads in the log files.
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------


%% Configuration
%----------------------------------------
% set configuration
cfg.expr.startexp	= 'bof';
cfg.expr.stopexp	= 'eof';
cfg.expr.startrun	= 'bof';
cfg.expr.stoprun	= '\r\n\r\n';
cfg.expr.vars       = '\t';
cfg.expr.timevar    = '^tim_';
cfg.timscale        = 1e-3;

% check configuration
cfg = km_setcfg(cfg,{'subj','sess','dirlog'});

%% Read logfile
%----------------------------------------
% find logfile
fname = dir(fullfile(cfg.dir.log,[cfg.prepad cfg.subj cfg.postpad]));
fname = cellfun(@(x) fullfile(cfg.dir.log,x),{fname.name},'UniformOutput',false);

% read the logfile, based on start and stop expressions
tcfg                = [];
tcfg.logfile        = fname;
tcfg.expr           = cfg.expr;
tcfg.allowstrings   = 'yes';
[hdr,logdata] = km_readlog_text(tcfg);

% capture variable names
tok = regexp(hdr{1},'\r','split');
vars = regexp(tok{1},cfg.expr.vars,'split');

%% Adapt logfile to fit standards
%----------------------------------------
nrun = length(logdata);
for r = 1:nrun
    % recode strings to numbers
    %idx_P1 = strcmp('corr_P1',vars);
    %idx_P2 = strcmp('corr_P2',vars);
    %tmp = double(strcmp('player 1 right',logdata{r}(:,idx_P1)));
    %tmp(strcmp('reset',logdata{r}(:,idx_P1))) = -1;
    %logdata{r}(:,idx_P1) = num2cell(tmp);
    %tmp = double(strcmp('player 2 right',logdata{r}(:,idx_P2)));
    %tmp(strcmp('reset',logdata{r}(:,idx_P2))) = -1;
    %logdata{r}(:,idx_P2) = num2cell(tmp);
    
    % convert cell to matrix
    logdata{r} = cell2mat(logdata{r});
        
    % get begin of exp time
    idx_tim_exp_begin = strcmp(vars,'tim_trl_ini');
    tim_exp_beg = logdata{r}(1,idx_tim_exp_begin);
    
    % correct timing based on the triggers
    idx_timvars = ~cellfun(@isempty,regexp(vars,cfg.expr.timevar,'once'));
    logdata{r}(:,idx_timvars) = logdata{r}(:,idx_timvars) - tim_exp_beg;
    
    % and scale times from ms to s
    idx_timvars = regexp(vars,cfg.expr.timevar,'once');
    idx_timvars = cellfun(@(x) ~isempty(x),idx_timvars);
    logdata{r}(:,idx_timvars) = logdata{r}(:,idx_timvars)*cfg.timscale;
end

% rename the variables so that they can be used in matlab as fieldnames
%vars{strcmp(vars,'Trialnumber')} = 'trlnr';
%vars{strcmp(vars,'Action type')} = 'action';
%vars{strcmp(vars,'Target Coding')} = 'coding';
%vars{strcmp(vars,'Total Number configuration')} = 'totcfgnr';
%vars{strcmp(vars,'Number configuration')} = 'cfgnr';
%vars{strcmp(vars,'Difficulty trial')} = 'level';
%vars{strcmp(vars,'P1: Cfg 3 targets')} = 'cfg_P1';
%vars{strcmp(vars,'P2: Cfg 3 targets')} = 'cfg_P2';
%vars{strcmp(vars,'P1: Correct Target')} = 'goal_P1';
%vars{strcmp(vars,'P2: Correct Target')} = 'goal_P2';
%vars{strcmp(vars,'P1: Selected Target')} = 'sel_P1';
%vars{strcmp(vars,'P2: Selected Target')} = 'sel_P2';
%vars{strcmp(vars,'P1: X Coordinate')} = 'x_P1';
%vars{strcmp(vars,'P2: X Coordinate')} = 'x_P2';
%vars{strcmp(vars,'P1: Y Coordinate')} = 'y_P1';
%vars{strcmp(vars,'P2: Y Coordinate')} = 'y_P2';
%vars{strcmp(vars,'P1: Clr_Squares')} = 'color_sq_P1';
%vars{strcmp(vars,'P2: Clr_Circles')} = 'color_ci_P2';
%vars{strcmp(vars,'P2: Clr_Triangles')} = 'color_tr_P2';
%vars{strcmp(vars,'P1: feedback')} = 'corr_P1';
%vars{strcmp(vars,'P2: feedback')} = 'corr_P2';
%vars{strcmp(vars,'P1: Reactiontime')} = 'rt_P1';
%vars{strcmp(vars,'P2: Reactiontime')} = 'rt_P2';
%vars{strcmp(vars,'P1: pointing time')} = 'mt_P1';
%vars{strcmp(vars,'P2: pointing time')} = 'mt_P2';
%vars{strcmp(vars,'P1: touch time')} = 'ht_P1';
%vars{strcmp(vars,'P2: touch time')} = 'ht_P2';
%vars{strcmp(vars,'P1: back time')} = 'bt_P1';
%vars{strcmp(vars,'P2: back time')} = 'bt_P2';
%vars{strcmp(vars,'Trial Begintime')} = 'tim_trl_beg';
%vars{strcmp(vars,'Trial Endtime')} = 'tim_trl_end';
%vars{strcmp(vars,'Trial Duration')} = 'dur_trl';
vars{strcmp(vars,'tim_exp_begin')} = 'tim_exp_beg';

for r = 1:nrun
    % RT from the logfile is the RT combined with MT
    %logdata{r}(:,strcmpi('rt_P1',vars)) = logdata{r}(:,strcmpi('rt_P1',vars)) - logdata{r}(:,strcmpi('mt_P1',vars));
    %logdata{r}(:,strcmpi('rt_P2',vars)) = logdata{r}(:,strcmpi('rt_P2',vars)) - logdata{r}(:,strcmpi('mt_P2',vars));
    % change coding of selected targets: [4 5 6] > [1 2 3]
    %logdata{r}(:,strcmpi('goal_P1',vars)) = logdata{r}(:,strcmpi('goal_P1',vars)) - 3;
end

% code post-error trials
idx_corr_P2 = find(strcmpi('corr_P2',vars));
vars = [vars(1:idx_corr_P2) {'posterr'} vars((idx_corr_P2+1):end)];
for r = 1:nrun
    error_P2 = logdata{r}(:,idx_corr_P2) == 0;
    posterr = [0; error_P2(1:end-1)] + 1;
    logdata{r} = [logdata{r}(:,1:idx_corr_P2) posterr logdata{r}(:,(idx_corr_P2+1):end)];
end

% adapt the time stamps in the logfile
%idx_beg = find(strcmpi('tim_trl_beg',vars));
%for r = 1:nrun
%    new_tim_trl_beg = logdata{r}(1:(end-1),strcmpi('tim_trl_end',vars));
%    new_tim_trl_beg = [logdata{r}(1,idx_beg); new_tim_trl_beg];
%    % add picture presentation duration
%    new_tim_trl_beg(2:end) = new_tim_trl_beg(2:end) + 0.017;
%    % store
%    logdata{r} = [logdata{r}(:,1:(idx_beg-1)) new_tim_trl_beg logdata{r}(:,idx_beg:end)];
%end
%vars = [vars(1:idx_beg) {'tim_stim'} vars((idx_beg+1):end)];

% add extra time stamps in the logfile
%[~,idx] = ismember({'tim_stim','rt_P1','mt_P1','ht_P1','bt_P1','rt_P2','mt_P2','ht_P2','bt_P2'},vars);
%new_vars = {'tim_but_u_P1','tim_sel_d_P1','tim_sel_u_P1','tim_but_d_P1','tim_but_u_P2','tim_sel_d_P2','tim_sel_u_P2','tim_but_d_P2'};
%vars = [vars(1:(idx_beg+1)) new_vars vars((idx_beg+2):end)];
%for r = 1:nrun
%    timdat = logdata{r}(:,idx);
%    % add picture presentation duration
%    timdat(:,2) = timdat(:,2) + 0.017;
%    % store
%    timdat = cumsum(timdat,2);
%    new_logdata = timdat(:,2:end);
%    logdata{r} = [logdata{r}(:,1:(idx_beg+1)) new_logdata logdata{r}(:,(idx_beg+2):end)];
%end

% store information in log structure
log.vars        = vars;
log.data        = logdata;
log.runnr       = 1:nrun;

