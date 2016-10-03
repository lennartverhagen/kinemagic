function [trl vars] = km_trldef_log(cfg)
%--------------------------------------------------------------------------
%
% See also KM_TRLDEF
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% set configuration
if ~isfield(cfg,'log'),             error('no log provided');               end
if ~isfield(cfg,'beg'),             cfg.beg         = 't_beg';              end
if ~isfield(cfg,'stim'),            cfg.stim        = cfg.beg;              end
if ~isfield(cfg,'end'),             cfg.end         = 't_end';              end
if ~isfield(cfg,'prepad'),          cfg.prepad      = 0;                    end     % prepad in seconds relative to begin of trial
if ~isfield(cfg,'postpad'),         cfg.postpad     = 0;                    end     % postpad in seconds relative to end of trial
if ~isfield(cfg,'offset'),          cfg.offset      = 0;                    end     % offset in seconds relative to stimulus

% variable names
vars = cfg.log.vars;

% loop over runs
if isfield(cfg.log,'runnr')
    runnr = cfg.log.runnr;
else
    runnr = 1;
end
trl = cell(1,max(runnr));
for r = 1:length(cfg.log.data)

    % log data
    logdata = cfg.log.data{r};

    % define trial begin, end and offset
    idx_beg     = strcmpi(vars,cfg.beg);
    idx_stim	= strcmpi(vars,cfg.stim);
    idx_end     = strcmpi(vars,cfg.end);
    if ~any(idx_beg), error('KM:trldef','Trial begin variable not found: %s',cfg.beg); end
    if ~any(idx_stim), error('KM:trldef','Trial stimulus variable not found: %s',cfg.stim); end
    if ~any(idx_end), error('KM:trldef','Trial end variable not found: %s',cfg.end); end
    timbeg      = logdata(:,idx_beg)+cfg.prepad;
    timend      = logdata(:,idx_end)+cfg.postpad;
    timstim     = logdata(:,idx_stim)+cfg.offset;
    trl{runnr(r)} = [timbeg timend timstim logdata];
        
end

% update variable names
vars = [{'trlbeg' 'trlend' 'trloff'} vars];