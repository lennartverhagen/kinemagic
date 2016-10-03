function [hdr,logdata] = km_readlog_text(cfg)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% check configuration
if ~isfield(cfg,'logfile'),         error('no logfile specified');          end
if ~iscell(cfg.logfile),            cfg.logfile = {cfg.logfile};            end
nlogfile = length(cfg.logfile);
for f = 1:nlogfile
    if ~exist(cfg.logfile{f},'file'),	error('logfile does not exist');	end
end
if ~isfield(cfg,'expr'),            cfg.expr            = struct;           end
if isfield(cfg.expr,'start')
    cfg.expr.startexp	= cfg.expr.start;
    cfg.expr.startrun   = cfg.expr.start;
    cfg.expr            = rmfield(cfg.expr.start);
end
if isfield(cfg.expr,'stop')
    cfg.expr.stopexp    = cfg.expr.stop;
    cfg.expr.stoprun    = cfg.expr.stop;
    cfg.expr            = rmfield(cfg.expr.stop);
end
if ~isfield(cfg.expr,'startexp'),	cfg.expr.startexp	= 'bof';            end
if ~isfield(cfg.expr,'stopexp'),    cfg.expr.stopexp    = 'eof';           	end
if ~isfield(cfg.expr,'startrun'),   cfg.expr.startrun	= 'bof';           	end
if ~isfield(cfg.expr,'stoprun'),    cfg.expr.stoprun    = 'eof';            end
if ~isfield(cfg,'allowstrings'),    cfg.allowstrings    = 'no';             end

% loop over files
hdr = cell(1,nlogfile);
logdata = cell(1,nlogfile);
for f = 1:nlogfile
    
    % open logfile and read as a string
    fid = fopen(cfg.logfile{f},'r');
    F = fread(fid);
    S = char(F');
    fclose(fid);
    
    % look for last instances of experiment start and stop expressions
    if isempty(cfg.expr.startexp) || strcmpi(cfg.expr.startexp,'bof')
        i_start = 1;
    else
        i_start	= max(regexp(S,cfg.expr.startexp));
    end
    if isempty(cfg.expr.stopexp) || strcmpi(cfg.expr.stopexp,'eof')
        i_stop	= length(S);
    else
        i_stop  = max(regexp(S,cfg.expr.stopexp));
    end
    
    % crop logfile to only last experiment in logfile
    S = S(i_start:i_stop);
    
    % look for all instances of run start expression
    if isempty(cfg.expr.startrun) || strcmpi(cfg.expr.startrun,'bof')
        i_start = 1;
    else
        i_start	= regexp(S,cfg.expr.startrun);
    end
    i_start = [i_start length(S)];
    
    % cut logfile into runs
    R = cell(1,length(i_start)-1);
    for r = 1:length(i_start)-1
        R{r} = S(i_start(r):i_start(r+1)-1);
    end
    
    % crop run logfiles to end of run expression
    for r = 1:length(R)
        if ~isempty(cfg.expr.stoprun) && ~strcmpi(cfg.expr.stoprun,'eof')
            i_stop  = regexp(R{r},cfg.expr.stoprun,'start','once');
            if isempty(i_stop),	i_stop = length(R{r})+1;    end
            R{r} = R{r}(1:i_stop-1);
        end
    end
    
    % separate the header from the log data
    hdr{f} = cell(1,length(R));
    logdata{f} = cell(1,length(R));
    for r = 1:length(R)
        i_hdrstart = 1;
        i_hdrstop = regexp(R{r},'\D\s+[-\d]+\t[-\d]+','start','once')+1;
        hdr{f}{r} = R{r}(i_hdrstart:i_hdrstop);
        
        i_logstart = regexp(R{r},'[-\d]+\t[-\d]+','start','once');
        i_logstop = max(regexp(R{r},'\t[-\d]+','end'));
        tmp = R{r}(i_logstart:i_logstop);
        if istrue(cfg.allowstrings)
            tok = regexp(tmp,'[\r\n]+','split')';
            tok = cellfun(@(x) regexp(x,'\t','split'),tok,'UniformOutput',false);
            tok = vertcat(tok{:});
            toknum = cellfun(@(x) str2num(x),tok,'UniformOutput',false);
            idx = cellfun(@(x) isempty(x),toknum,'UniformOutput',true);
            toknum(idx) = tok(idx);
            logdata{f}{r} = toknum;
        else
            tmpdat = str2num(tmp);
            if ~isempty(tmpdat)
                logdata{f}{r} = tmpdat;
            else
                tok = sprintf('\r\n%s',tmp);
                tok = regexp(tok,'\r\n([\t-\d]+)','tokens');
                tok = cellfun(@(x) str2num(x{1}),tok,'UniformOutput',false);
                if all(cellfun(@(x) isempty(x),tok))
                    logdata{f}{r} = tmp;
                else
                    logdata{f}{r} = tok;
                end
            end
        end
    end
end

% collapse over files
hdr = [hdr{:}];
logdata = [logdata{:}];
    
% remove all empty headers and logdata's
idx = ~cellfun(@isempty,logdata);
hdr = hdr(idx);
logdata = logdata(idx);

