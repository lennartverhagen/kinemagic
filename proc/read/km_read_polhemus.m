function [cfg,data] = km_read_polhemus(cfg)
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
if ~isfield(cfg,'posscale'),	cfg.posscale	= 1e-2;                     end	% 0.1 mm to 1 m
if ~isfield(cfg,'oriscale'),	cfg.oriscale	= 1;                        end % no scaling
if ~isfield(cfg,'timscale'),	cfg.timscale	= 1e-3;                     end % 1 ms to 1 s

% check configuration
cfg = km_setcfg(cfg,{'subj','sess','dirdata'});


%% Get data files
%----------------------------------------
% first look for binary .dat file
fname = dir(fullfile(cfg.dir.data,sprintf('*%s%s*.dat',cfg.subj,cfg.sess)));

% then go look for the binary .tmot file (header + data)
if isempty(fname)
    fname = dir(fullfile(cfg.dir.data,sprintf('*%s%s*.tmot',cfg.subj,cfg.sess)));
end

% maybe try to find a ASCII .txt file
if isempty(fname)
    fname = dir(fullfile(cfg.dir.data,sprintf('*%s%s*.txt',cfg.subj,cfg.sess)));
end

% finally settle for the comma-seperated ASCII .csv file
if isempty(fname)
    fname = dir(fullfile(cfg.dir.data,sprintf('*%s%s*.csv',cfg.subj,cfg.sess)));
end

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
for r = 1:cfg.info.nfiles
    
    % determine filetype
    [~,~,ext] = fileparts(fname(r).name);
    filetype = lower(ext(2:end));
    fname(r).name = fullfile(cfg.dir.data,fname(r).name);
    
    if ismember(filetype,{'dat','tmot','txt'})
        fid = fopen(fname(r).name, 'r');
    end
    
    switch filetype
        case 'dat'
            
            % file specifics
            hdr = 8;        % Size of frame header in bytes
            mnr = 2;        % Number of bytes to skip before "marker" byte
            cnr = 6;        % Number of bytes to skip before "response body" bytes
            tnr = hdr;      % Number of bytes to skip before "timestamp" bytes
            snr = hdr + 8;  % Number of bytes to skip before "sync" byte
            dnr = hdr + 12; % Number of bytes to skip before "distortion" byte
            xnr = hdr + 12; % Number of bytes to skip before "xpos" bytes
            
            % Get length of response body
            fseek(fid,cnr,0);           % skip frame tag, station number, initation command, error indicator and reserved byte
            c = fread(fid,2,'int8');	% get number of bytes in response body
            
            % check if distortion byte is included
            get_dist = 0;
            if c(1)==42
                get_dist = 1;
                xnr = xnr + 4;
            end
            
            % Get marker
            fseek(fid,mnr,-1);          % fo to first marker byte
            skip = 3+2+c(1)+2;          % skip initiation command, error indicator, reserved byte, response body (c) and frame tag
            m = fread(fid,'*int8',skip);
            
            % Get time and frame stamps
            fseek(fid,tnr,-1);          % go to first time stamp
            skip = hdr + c(1)-(2*4);
            tim = fread(fid,[2,length(m)],'2*uint32',skip)';
            
            % Get external sync byte
            fseek(fid,snr,-1);          % go to first extarnal sync byte
            skip = hdr + c(1)-4;
            sync = fread(fid,length(m),'1*int32',skip);
            
            % Get distortion byte
            if get_dist
                fseek(fid,dnr,-1);      % go to first extarnal dist byte
                skip = hdr + c(1)-4;
                dist = fread(fid,length(m),'1*int32',skip);
            end
            
            % Get pos and ori
            fseek(fid,xnr,-1);          % go to first x position byte
            skip = hdr + c(1)-(6*4);
            M = fread(fid,[6,length(m)],'6*float32',skip)';
            
            % Close data file
            fclose(fid);
            
            % separate markers from data matrix
            dimord = 'marker_time_axis';
            u = unique(m);
            nmark = length(u);
            nsamp = length(m)/nmark;
            naxes = size(M,2);
            N = nan(nmark,nsamp,naxes);
            for i = 1:nmark
                N(u(i),:,:) = M(i:nmark:nmark*nsamp,:);
            end
            
        case 'tmot'
            
            error('The Polhemus standard .tmot file type is not yet supported, but it is just the binary file with an extra header...')
            
        case 'txt'
            
            C = textscan(fid,'%f%s%f%f%f%f%f%f%f%f');
            fclose(fid);
            
            for c = 1:length(C)
                if isnumeric(C{c})
                    M(:,c) = C{c};
                else
                    A = str2mat(C{c});
                    M(:,c) = hex2dec(A(:,3:end));
                end
            end
            
            error('Polhemus txt file reader is under construction, switch to csv if possible');
            
        case 'csv'
            
            M = dlmread(fname(r).name,',',1,0);
            ncol = size(M,2);
            
            % FIXME, read the first line with the variable names
            % option 1:
            %vars = dlmread(fname(r).name,',',[0 0 1 ncol]);
            % option 2:
            fid = fopen(fname(r).name, 'r');
            vars = regexp(fgetl(fid),',','split');
            fclose(fid);
            if length(vars) ~= ncol
                vars = {'sensor','sample','sync','x','y','z'};
            end
            vars = lower(vars);
            
            % get marker, time, sample and sync
            idx_marker  = ismember(vars,'sensor');
            idx_sample  = ismember(vars,'sample');
            idx_sync    = ismember(vars,'sync');
            if ~any(idx_sync)
                idx_sync    = ismember(vars,'marker');
            end
            idx_dat     = ~cellfun(@isempty,regexpi(vars,'^[xyz]|[xyz]$'));
            m = M(:,idx_marker);
            Fs = 240*cfg.timscale;
            tim = [M(:,idx_sample)/Fs M(:,idx_sample)];
            sync = M(:,idx_sync);
            
            % alternatively, use the (less precise) timestamp
            % idx_time    = ismember(vars,'timestamp');
            % tim = [M(:,idx_time) M(:,idx_sample)];
            
            % store axis names
            if ~isfield(data,'axis')
                tmp = regexpi(vars(idx_dat),'^([xyz])|([xyz])$','tokens');
                data.axis = unique(cellfun(@(x) x{1}(1),tmp));
            end
            
            % separate markers from data matrix
            dimord = 'marker_time_axis';
            u = unique(m);
            nmark = length(u);
            nsamp = length(m)/nmark;
            naxes = sum(idx_dat);
            N = nan(nmark,nsamp,naxes);
            for i = 1:nmark
                N(u(i),:,:) = M(i:nmark:nmark*nsamp,idx_dat);
            end
    end
    
    % organize dimensions
    N = permute(N,finddim(dimord,cfg.dimord));
    
    % store position and orientation data
    data.pos{r} = N(:,1:3,:) * cfg.posscale;
    if size(N,2) > 3
        data.ori{r} = N(:,4:6,:) * cfg.oriscale;
    end
    
    % store time, sample and event stamps
    data.time{r} = tim(1:nmark:(nmark*nsamp),1)*cfg.timscale;
    data.samp{r} = tim(1:nmark:(nmark*nsamp),2);
    data.event{r} = sync(1:nmark:(nmark*nsamp));
    
    % store dimensions in configuration
    if r == 1
        data.dimord = cfg.dimord;
        cfg.info.nmark          = size(data.pos{r},1);
        cfg.info.nposaxes       = size(data.pos{r},2);
        if isfield(data,'ori')
            cfg.info.noriaxes	= size(data.ori{r},2);
        else
            cfg.info.noriaxes   = 0;
        end
    end
    cfg.info.nsamp(r) = size(data.pos{r},3);
    
    % get sample frequency and experiment duration
    twin = max(data.time{r}) - min(data.time{r});
    data.fsample(r) = (cfg.info.nsamp(r)-1)/twin;
    
    % get run number
    if ~isempty(cfg.sess)
        expr = [cfg.sess '(\d)'];
    else
        expr = '_run(\d)';
    end
    runname = regexp(fname(r).name,expr,'tokens');
    if ~isempty(runname)
        data.run(r) = str2double(runname{end});
    else
        data.run(r) = r;
    end
    
    % get experiment duration
    dur = round(twin/60);
    cfg.info.expdur{r} = sprintf('%1.0d minutes',dur);
    
end