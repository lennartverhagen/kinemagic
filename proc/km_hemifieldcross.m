function [cfg,data] = km_hemifieldcross(cfg,data)
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


method = cfg.hemifieldcross;
switch lower(method)
    
    % loose hemifield crossing correction (only y and z axes)
    %----------------------------------------
    case {'loose','soft'}
        
        % check number of axis
        if length(data.axis) < 3
            error('hemifield crossing correction can only be done robustly on 3 axis');
        end
        
        % store the sign of the first axis of the first run
        defaultsign = zeros(length(data.label),length(data.axis));
        
        % loop over runs
        for r = 1:length(data.pos)
            
            fprintf(1,'\n  Hemisphere crossing (%s) in run %1.0f:\n',method,r);
            
            % loop over markers
            for m = 1:length(data.label)
                
                % get data
                pos = squeeze(data.pos{r}(m,:,:));
                
                % fill NaN blocks with last non-NaN value (or the first
                % non-NaN value)
                for a = 1:size(pos,1)
                    tmp = [pos(a,:) 1];
                    idx_nan = logic2idx(isnan(tmp));
                    for i = 1:size(idx_nan)
                        if idx_nan(i,1) == 1
                            pos(a,idx_nan(i,1):idx_nan(i,2)) = tmp(idx_nan(i,2)+1);
                        else
                            pos(a,idx_nan(i,1):idx_nan(i,2)) = tmp(idx_nan(i,1)-1);
                        end
                    end
                end
                
                % get sign (and shifted sign)
                sg = squeeze(sign(pos(:,1:end-1)));
                sg_s = squeeze(sign(pos(:,2:end)));
                
                % determine default sign of the axes
                if r == 1
                    defaultsign(m,:) = mode(sign(pos)');
                end
                
                % Find indeces where the axes flip.
                flip = [nan(size(sg,1),1) single(sg ~= sg_s)];
                
                % determine if the first sample is already flipped or
                % not...
                flip(:,1) = sign(pos(:,1)) ~= defaultsign(m,:)';
                
                % Uniquely identify the axis.
                flip(1,flip(1,:)>0) = 2;
                flip(2,flip(2,:)>0) = 3;
                flip(3,flip(3,:)>0) = 4;
                flip = sum(flip,1);
                
                % find indeces where two axes flip, but one remains, these
                % are most likely the locations where the constant axis
                % bounces of the hemifield border and makes the other two
                % flip signs.
                ic = find(flip>4 & flip<9);
                
                % Identify the bouncing axis
                if isempty(ic)
                    axstr = '';
                else
                    axissel = flip(ic);
                    axissel(axissel==7) = 1;
                    axissel(axissel==6) = 2;
                    axissel(axissel==5) = 3;
                    baxis = unique(axissel);
                    if length(baxis) > 1
                        baxstr = sprintf('%s',data.axis{baxis});
                        baxis = mode(axissel);
                        if sum(axissel==baxis) <= length(axissel)/2
                            error('KM:HemiCrossMultiAxes','More than one axis (%s) seem to be bouncing... but a reliable axis could not be determined',baxstr);
                        else
                            saxstr = data.axis{baxis};
                            warning('KM:HemiCrossMultiAxes','More than one axis (%s) seems to be bouncing... but the %s axis bounces the most',baxstr,saxstr);
                        end
                    end
                    axstr = sprintf(' on the %s axis',data.axis{baxis});
                end
                
                % report
                fprintf(1,'    Marker %s: %1.0f times%s\n',data.label{m},length(ic)/2,axstr);
                
                % continue if no flips are detected
                if isempty(ic), continue;   end
                
                % make flips of the right length
                if isodd(length(ic))
                    ic = [ic size(pos,2)+1];
                end
                
                % correct flips
                idx_on = ic(1:2:end);
                idx_off = ic(2:2:end)-1;
                for i = 1:length(idx_on)
                    idx = idx_on(i):idx_off(i);
                    pos(:,idx) = -pos(:,idx);
                end
                
                % Some kinematic equipment have a low-pass filter
                % implemented. This will lead to a step-response after a
                % sign flip, lasting a few samples. These samples should be
                % interpolated to ensure smooth data. Here I chose to
                % interpolate 10 samples. That is actually quite liberal!
                ic = sort([idx_on idx_off]);
                pos = interp_flip(pos,ic,10);
                
                % set data
                data.pos{r}(m,:,:) = pos;
                
                % update the default sign of the axes
                if r == 1
                    defaultsign(m,:) = mode(sign(pos)');
                end
                
            end
            
        end
        
    % strict hemifield crossing correction (all axes that change sign)
    %----------------------------------------
    case {'strict','hard'}
        
        error('under construction');    % FIXME
        
        cfg.hemifieldcross = 'soft';
        data = km_hemifieldcross(cfg,data);
        cfg.hemifieldcross = 'hard';
        
        for r = 1:length(data.pos)
            fprintf('\n  Hemisphere crossing (%s) in run %1.0f:',method,r);
            
            pos = KMD.data.raw{f}.pos;
            for m = 1:size(pos,1)
                sg = sign(pos(m,:,1:end-1));
                sg_s = sign(pos(m,:,2:end));
                
                ax = 'xyz';
                for d = 2:3
                    mc = find( sg(:,d) ~= sg_s(:,d) & ~isnan(sg(:,d)) &...
                        ~isnan(sg_s(:,d)) ) + 1;
                    disp(sprintf('    Bird %1.0f (%s-axis): %1.0f times',m,ax(d),length(mc)));
                    
                    idx = [];
                    
                    while true
                        i = length(mc);
                        if i >= 2
                            idx = [idx mc(1):mc(2)-1];
                            mc = mc(3:end);
                        elseif i == 1
                            idx = [idx mc(1):size(pos,1)];
                            mc = [];
                        elseif i < 1
                            break;
                        end
                    end
                    pos(idx,d,m) = pos(idx,d,m).*-1;
                end
            end
            KMD.data.raw{f}.pos = pos;
        end
        
    % not recognized
    %----------------------------------------
    otherwise
        error('Hemisphere crossing correction method %s not recognized',method);
    
end
%--------------------------------------------------------------------------


%% function interp_flip
%----------------------------------------
function pos = interp_flip(pos,ic,n)

% check input
if nargin < 3,  n = 10;	end
if n < 1,       return;	end
ic = ic(:);

% get on and offset of flips
art = [ic ic+n-1];

% go to sample domain and back to combine overlapping artifacts
nsamp = size(pos,2);
art = idx2logic(art,nsamp);

% padded artifact
npadsamp = 100;
idx = conv(single(art),ones(1,npadsamp),'same')>0;

% original samples
x = 1:nsamp;
xi = x(idx)';

% time without artifact
x = xi(~art(idx));

% data without artifact
y = pos(:,idx)';
y = y(~art(idx),:);

% interpolate
yi = interp1(x,y,xi,'pchip','extrap');

% reshape data to fit original
pos(:,idx) = yi';