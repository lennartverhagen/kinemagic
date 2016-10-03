function [cfg,data] = km_filter(cfg,data)
%--------------------------------------------------------------------------
%
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


% check if the sampling frequency is stable
for r = 1:length(data.time)
    % determine time step
    diffx = diff(data.time{r});
    step = 1/data.fsample(r);
    if any(abs(diffx-step)>step/4)
        error('The sampling frequency is not stable enough for filtering, maybe you should first cut-out sections with missing samples (km_cutbreaks) and interpolate (km_interp).')
    end
end

% define artifacts if requested
if isfield(cfg.filter,'artfctdef')
    cfg.artfctdef = cfg.filter.artfctdef;
    cfg.filter.artifact = km_artfctdef(cfg,data);
end

% loop over time series
for ts = 1:length(cfg.filter.tseries)
    
    % create temporary configuration structure
    tcfg = cfg.filter;
    tcfg.tseries = tcfg.tseries{ts};
    
    % filter data
    [tcfg,data] = filterdat(tcfg,data);
    
end

% update configuration structure
tcfg.tseries = cfg.filter.tseries;
cfg.filter = tcfg;

% update dataset
if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
switch lower(cfg.filter.mode)
    case 'lowpass'
        modestr = 'lp';
    case 'highpass'
        modestr = 'hp';
    case 'bandpass'
        modestr = 'bp';
    case 'bandstop'
        modestr = 'bs';
    case 'dft'
        modestr = 'dft';
    otherwise
        modestr = 'other';
end
if numel(cfg.filter.freq) == 2
    filtstr = sprintf('FILT%s%dto%d',modestr,round(cfg.filter.freq(1)),round(cfg.filter.freq(2)));
else
    filtstr = sprintf('FILT%s%d',modestr,round(cfg.filter.freq(1)));
end
cfg.dataset = [cfg.dataset '_' filtstr];


%% function filterdat
%----------------------------------------
function [cfg,data] = filterdat(cfg,data)
% filter data (position, grip, velocity, etc)

% check input
if ~isfield(cfg,'artifact'),	cfg.artifact = [];	end

% return if nothing to do
if isfalse(cfg)
    return
end

% get settings
Fs      = data.fsample;
mode    = cfg.mode;	% 'lowpass'
Ffil    = cfg.freq;	% e.g. 15 for 'lowpass', 0.1 for 'highpass', [0.1 15] for 'bandpass'
type    = cfg.type;	% 'but'
N       = cfg.ord;	% 6 for 'but', 25 for 'fir'
dir     = cfg.dir;	% 'twopass'

% check filter frequency
switch mode
    case 'lowpass'
        Ffil = max(Ffil);
    case 'highpass'
        Ffil = min(Ffil);
    case 'bandpass'
        Ffil = [min(Ffil) max(Ffil)];
    case 'bandstop'
        Ffil = [min(Ffil) max(Ffil)];
end

% get function handles
switch type
    case 'but'
        filterfun = @butter;
    case 'fir'
        filterfun = @fir1;
    case {'cheby','cheby2'}
        filterfun = @cheby2;
    case 'cheby1'
        filterfun = @cheby1;
end

% get data based on cfg.tseries
if isfield(data,cfg.tseries)
    dat = data.(cfg.tseries);
else
    tmptseries = strrep(cfg.tseries,'grip','');
    tmptseries = strrep(tmptseries,'grp','');
    if isfield(data,tmptseries)
        dat = data.grip.(tmptseries);
    else
        error('Parameter ''%s'' not present in the data.',cfg.tseries)
    end
end

% loop over runs
for r = 1:length(dat)
    
    % Nyquist frequency
    Fn = Fs(r)/2;

    % compute filter coefficients
    switch mode
        case 'highpass'
            [B, A] = feval(filterfun, N, Ffil/Fn, 'high');
        case 'bandstop'
            [B, A] = feval(filterfun, N, Ffil/Fn, 'stop');
        otherwise
            [B, A] = feval(filterfun, N, Ffil/Fn);
    end
    
    % original data
    y = permute(dat{r},[3 2 1]);
    
    % apply filter to data
    if ~isempty(cfg.artifact)
        % select artifact free data segments
        art = km_time2logic(cfg.artifact{r},data.time{r});
        idx = km_logic2idx(~art);
        
        % loop over markers
        yi = nan(size(y));
        for m = 1:size(y,3)
            
            % loop over data segments
            for i = 1:size(idx,1)
                yy = y(idx(i,1):idx(i,2),:,m);
                switch dir
                    case 'onepass'
                        yyi = filter(B, A, yy)';
                    case 'onepass-reverse'
                        yy  = flipud(yy);
                        yyi = filter(B, A, yy)';
                        yyi = flipud(yyi);
                    case 'twopass'
                        yyi = filtfilt(B, A, yy);
                end
                yi(idx(i,1):idx(i,2),:,m) = yyi;
            end
        end
        
    elseif any(isnan(y(:)))
        % loop over markers
        yi = nan(size(y));
        for m = 1:size(y,3)
            % select NaN free data segments
            art = any(isnan(y(:,:,m)),2);
            idx = logic2idx(~art);
            
            % loop over data segments
            for i = 1:size(idx,1)
                yy = y(idx(i,1):idx(i,2),:,m);
                switch dir
                    case 'onepass'
                        yyi = filter(B, A, yy)';
                    case 'onepass-reverse'
                        yy  = flipud(yy);
                        yyi = filter(B, A, yy)';
                        yyi = flipud(yyi);
                    case 'twopass'
                        yyi = filtfilt(B, A, yy);
                end
                yi(idx(i,1):idx(i,2),:,m) = yyi;
            end
        end
        
    else
        % no artifact and no NaNs
        % loop over markers
        yi = nan(size(y));
        for m = 1:size(y,3)
            yy = y(:,:,m);
            switch dir
                case 'onepass'
                    yyi = filter(B, A, yy)';
                case 'onepass-reverse'
                    yy  = flipud(yy);
                    yyi = filter(B, A, yy)';
                    yyi = flipud(yyi);
                case 'twopass'
                    yyi = filtfilt(B, A, yy);
            end
            yi(:,:,m) = yyi;
        end
    end
    
    % store filtered data
    dat{r} = permute(yi,[3 2 1]);
    
end

% store data
if ~exist('tmptseries','var')
    data.(cfg.tseries) = dat;
else
    data.grip.(tmptseries) = dat;
end
