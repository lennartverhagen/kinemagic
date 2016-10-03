function [cfg,data] = km_read_joystick(cfg)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------


%% Configuration
%----------------------------------------
% set configuration
if ~isfield(cfg,'dimord'),      cfg.dimord      = 'marker_axis_time';       end	% the first dimension is the marker, the second the spatial dimension (x,y,z) and the third the time
if ~isfield(cfg,'posscale'),	cfg.posscale	= 2e-5;                     end	% 0.05 mm to 1 m
if ~isfield(cfg,'timscale'),	cfg.timscale	= 1e-3;                     end % 1 ms to 1 s

% check configuration
cfg = km_setcfg(cfg,{'subj','sess','dirdata'});

% check function specific configuration
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
if ~isfield(cfg.expr,'startexp'),	cfg.expr.startexp	= 'Time of start Experiment:';	end
if ~isfield(cfg.expr,'stopexp'),    cfg.expr.stopexp    = 'eof';                        end
if ~isfield(cfg.expr,'startrun'),   cfg.expr.startrun	= '--- Block:';                     end
if ~isfield(cfg.expr,'stoprun'),    cfg.expr.stoprun    = '\r\n\r\n';                   end

% FIXME: calibration string


%% Get data files
%----------------------------------------
% first look for ASCII .txt file
fname = dir(fullfile(cfg.dir.data,sprintf('*%s%s*_pos.txt',cfg.subj,cfg.sess)));

% find the logfile(s) and remove them from the txt-file list
idx = [];
for f = 1:length(fname)
    if isempty(strfind(fname(f).name,'log'))
        idx(end+1) = f;
    end
end
fname = fname(idx);

% give error when no data files are found
if isempty(fname)
    error('No data files found in %s',cfg.dir.data);
end

% store some experimental information
cfg.info.nfiles	= length(fname);
cfg.info.date	= fname(1).date;


%% Read data
%----------------------------------------
data = [];
tr = 1;
for r = 1:cfg.info.nfiles
    
    % open logfile and read as a string
    fname(r).name = fullfile(cfg.dir.data,fname(r).name);
    fid = fopen(fname(r).name, 'r');
    F = fread(fid);
    S = char(F');
    fclose(fid);
    
    % look for last instances of experiment start and stop expressions
    if isempty(cfg.expr.startexp) || strcmpi(cfg.expr.startexp,'bof')
        i_start = 1;
    else
        i_start	= max(regexp(S,cfg.expr.startexp));
    end
    if isempty(i_start),    i_start = 1;    end
    if isempty(cfg.expr.stopexp) || strcmpi(cfg.expr.stopexp,'eof')
        i_stop	= length(S);
    else
        i_stop  = max(regexp(S,cfg.expr.stopexp));
    end
    if isempty(i_stop),    i_stop = 1;    end
    
    % crop logfile to only last experiment in logfile
    S = S(i_start:i_stop);
    
    % get experimental trigger time
    Texp = str2double(regexp(S,'\d*','once','match'));
    
    % look for all instances of run start expression
    if isempty(cfg.expr.startrun) || strcmpi(cfg.expr.startrun,'bof')
        i_start = 1;
    else
        i_start	= regexp(S,cfg.expr.startrun);
    end
    i_start = [i_start length(S)];

    % cut logfile into runs
    R = cell(1,length(i_start)-1);
    for t = 1:length(i_start)-1
        R{t} = S(i_start(t):i_start(t+1)-1);
    end
    
    % crop run logfiles to end of run expression
    for t = 1:length(R)
        if ~isempty(cfg.expr.stoprun) && ~strcmpi(cfg.expr.stoprun,'eof')
            i_stop  = regexp(R{t},cfg.expr.stoprun,'start','once');
            if isempty(i_stop),	i_stop = length(R{t})+1;    end
            R{t} = R{t}(1:i_stop-1);
        end
    end
    
    % separate the header from the log data
    for t = 1:length(R)
        i_hdrstart = 1;
        i_hdrstop = max(regexp(R{t},'---','end'))+1;
        hdr = R{t}(i_hdrstart:i_hdrstop);
        
        % FIXME: this opens up the opportunity to take into account the
        % number of markers, axis, time offset and data chunks
                
        % get fixation and stimulus times
        Tfix = regexp(hdr,'time_fix: (\d*)','tokens');
        Tstim = regexp(hdr,'time_pict: (\d*)','tokens');
        Tfix = str2double(Tfix{1});
        Tstim = str2double(Tstim{1});
        
        % read posdata per line
        i_line = regexp(R{t},'\r\n','start');
        i_line = [i_line length(R{t})];
        posdata = cell(1,length(i_line)-1);
        for il = 1:length(i_line)-1
            posdata{il} = str2num(R{t}(i_line(il):i_line(il+1)));
        end
        
        % calculate time offset of the fixation relative to the stimulus
        % usefull to calculate the time relative to the stimulus
        %fixoffset = Tfix-Tstim;
        
        % time relative to begin of stimulus
        if length(posdata) == 6
            % time relative to the stimulus
            %tim = [posdata{1}+fixoffset posdata{3} posdata{5}];
            % time relative to the begin of experiment
            tim = [posdata{1}+Tfix posdata{3}+Tstim posdata{5}+Tstim] - Texp;
        else
            % time relative to the stimulus
            %tim = [posdata{1}+fixoffset posdata{3}];
            % time relative to the begin of experiment
            tim = [posdata{1}+Tfix posdata{3}+Tstim] - Texp;
        end
        
        % position data (for one axis)
        if length(posdata) == 6
            pos = [posdata{2} posdata{4} posdata{6}];
        else
            pos = [posdata{2} posdata{4}];
        end
        
        % reshape pos into the correct dimensions (1 marker, 1 axis, many
        % time-points)
        dims = [1 1 length(pos)];
        pos = reshape(pos,dims);
        
        % store tim and pos
        data.time{tr} = tim * cfg.timscale;
        data.pos{tr}  = pos * cfg.posscale;
        
        % get sample frequency and experiment duration
        twin = max(data.time{tr}) - min(data.time{tr});
        data.fsample(tr) = (size(data.pos{tr},2)-1)/twin;
    	
        % get run number  
        data.run(tr) = r;
        
        % increment trial counter
        tr = tr + 1;
        
    end
        
end