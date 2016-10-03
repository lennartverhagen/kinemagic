function [cfg,param,tseries] = km_plot_descriptives(cfg,Sparam,SCI,Sepe,Stseries)
%--------------------------------------------------------------------------
% calculate descriptives
%
% See also KM_PLOT
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% return if data is already in output format
if ~iscell(Sparam) && ~iscell(Stseries)
    param = Sparam;
    tseries = Stseries;
    return
end

nsubj = length(cfg.subj);
fname  = cfg.plot.fname;
nfact  = length(fname);
if ~isempty(cfg.plot.tseries)
    cname = cfg.plot.cname;
    ncond = length(cname);
end
nparam = length(cfg.plot.param);
ntseries = length(cfg.plot.tseries);

if ~isempty(SCI{1})
    nCI = length(cfg.plot.CI.param);
else
    nCI = 0;
end

if ~isempty(Sepe{1})
    nepe = length(cfg.plot.epe.marker);
else
    nepe = 0;
end

%% calculate descriptives
%----------------------------------------
groupdes = cfg.plot.descrip.group;
subjdes = cfg.plot.descrip.subj;
errdes = cfg.plot.descrip.err;
% parameters
param = struct;

% loop over parameters
for p = 1:nparam+nCI+nepe
    
    % parameter name
    isparam = false;
    isCI = false;
    isepe = false;
    if p <= nparam
        isparam = true;
        pnameload = cfg.plot.param{p};
        pnamesave = pnameload;
    elseif p <= nparam+nCI
        isCI = true;
        pnameload = cfg.plot.CI.param{p-nparam};
        pnamesave = [pnameload 'CI'];
    else
        isepe = true;
        pnameload = sprintf('epe_%s',cfg.plot.epe.marker{p-nparam-nCI});
        pnamesave = pnameload;
    end
    
    % loop over factors
    for f = 1:nfact
        
        % organize dimensions
        if isparam
            dims	= size(Sparam{1}.(fname{f}));
        elseif isCI
            dims	= size(SCI{1}.(fname{f}).(pnameload));
        elseif isepe
            dims	= size(Sepe{1}.(fname{f}).(cfg.plot.epe.marker{p-nparam}));
        end
        [~,cidx] = ismember(cfg.plot.contr{f},cfg.plot.fact{f});
        if cidx < 1,    cidx = [];  end
        [~,ridx] = ismember(cfg.plot.factrem{f},cfg.plot.fact{f});
        if ridx < 1,    ridx = [];  end
        
        % get mean, median or variance descriptives
        subjdat = nan(nsubj,prod(dims));
        for s = 1:nsubj
            % get data
            if isparam
                dat = cellfun(@(x) x(:,p),Sparam{s}.(fname{f})(:),'UniformOutput',false);
            elseif isCI
                dat = num2cell(SCI{s}.(fname{f}).(pnameload)(:));
            elseif isepe
                dat = num2cell(Sepe{s}.(fname{f}).(cfg.plot.epe.marker{p-nparam})(:));
            end
            % log-transform data if requested
            if any(strcmpi(subjdes,{'meanlog','medlog','varlog','semlog'}))
                % the line below is the original one
                dat = cellfun(@(x) log(x),dat,'UniformOutput',false);
                % HACK: line below only for gcd
                %dat = cellfun(@(x) log10((x+0.005)*200),dat,'UniformOutput',false);
                % HACK: line below only for mgv
                %dat = cellfun(@(x) log10((x+0.05)*20),dat,'UniformOutput',false);
                for i = 1:numel(dat)
                    dat{i}(~isfinite(dat{i})) = NaN;
                end
            end
            
            % calculate descriptive
            switch lower(subjdes)
                case {'mean','meanlog','logmean'}
                    subjdat(s,:)	= cellfun(@(x) nanmean(x),dat);
                case {'med','median','medlog','medianlog','logmed','logmedian'}
                    subjdat(s,:)	= cellfun(@(x) nanmedian(x),dat);
                case {'var','variance','varlog','logvar'}
                    subjdat(s,:)	= cellfun(@(x) nanvar(x),dat);
                case {'sem','semlog','logsem'}
                    subjdat(s,:)	= cellfun(@(x) nansem(x),dat);
                otherwise
                    error('subject data descriptive not supported');
            end
        end
        
        % calculate contrasts
        subjdat = calc_contr(subjdat,dims,cidx);
        
        % remove main factor effects
        subjdat = remove_maineff(subjdat,dims,ridx);
        
        % update dimensions
        dims(cidx) = [];
        if isempty(dims),       dims = [1 1];       end
        if length(dims) == 1,   dims = [1 dims];    end
        
        % log-transform descriptives if requested
        if any(strcmpi(subjdes,{'logmean','logmed','logvar','logsem'}))
            subjdat = log(subjdat);
        end
        
        % calculate error over repetitions, or over subjects
        if nsubj == 1
            
            % initiate temporary variable
            tmp = struct;
            
            % store descriptives
            tmp.dat	= reshape(subjdat,dims);
            
            if isempty(cidx)
                % calculate error descriptive
                switch lower(errdes)
                    case 'var'
                        tmp.err = reshape(cellfun(@(x) nanvar(x),dat),dims);
                    case 'sem'
                        tmp.err = reshape(cellfun(@(x) nansem(x),dat),dims);
                    otherwise
                end
                
                % store original data
                tmp.val        = reshape(dat,dims);
            else
                tmp.err = nan(size(tmp.dat));
                tmp.val        = tmp.dat;
            end
            
            % store tmp in correct field
            param.(pnamesave).(fname{f}) = tmp;
            
        else
            
            % initiate temporary variable
            tmp = struct;
            
            % calculate group descriptive
            switch lower(groupdes)
                case 'mean'
                    tmp.dat	= reshape(nanmean(subjdat,1),dims);
                case 'median'
                    tmp.dat	= reshape(nanmedian(subjdat,1),dims);
            end
            
            % calculate group error
            switch lower(errdes)
                case 'var'
                    tmp.err	= reshape(nanvar(subjdat,1),dims);
                case 'sem'
                    tmp.err	= reshape(nansem(subjdat,1),dims);
            end
            
            % store subject data
            tmp.subj = reshape(subjdat,[nsubj dims]);
            
            % store tmp in correct field
            param.(pnamesave).(fname{f}) = tmp;
            
            % calculate the descriptives on the difference
            if (length(dims) > 2 && any(dims(1:end-1) == 2)) || ...
                (length(dims) < 3 && any(dims == 2))
                
                % initiate temporary variable
                tmp = struct;
                
                % update the dimensions to calculate the difference
                subjdims = [nsubj dims];
                didx = find(dims == 2,1,'first')+1;
                diffdims = subjdims;
                diffdims(didx) = 1;
                
                % calculate the difference of the mean and median
                subjdat = reshape(diff(reshape(subjdat,subjdims),1,didx),diffdims);
                
                % update the dimensions
                dims = diffdims(2:end);
                if length(dims) == 1,   dims = [1 dims];    end
                
                % calculate group descriptive
                switch lower(groupdes)
                    case 'mean'
                        tmp.dat	= reshape(nanmean(subjdat,1),dims);
                    case 'median'
                        tmp.dat	= reshape(nanmedian(subjdat,1),dims);
                end
                
                % calculate group error
                switch lower(errdes)
                    case 'var'
                        tmp.err	= reshape(nanvar(subjdat,1),dims);
                    case 'sem'
                        tmp.err	= reshape(nansem(subjdat,1),dims);
                end
                
                % store tmp in correct field
                param.(pnamesave).(fname{f}).diff = tmp;
                
            end
        end

    end
    
    % add epe to param names
    if p > nparam
        cfg.plot.param{p} = pnamesave;
    end
    
end

% time series
tseries = struct;
for ts = 1:ntseries
    
    % get time-series name
    tsname = cfg.plot.tseries{ts};
    
    % store dimord, axis and marker
    tseries.(tsname).dimord = regexprep(Stseries{1}.(tsname).dimord,'^rpt_|^rept_','');
    tseries.(tsname).axis   = Stseries{1}.(tsname).axis;
    tseries.(tsname).marker = Stseries{1}.(tsname).marker;
    
    % loop over factors
    for f = 1:nfact
        
        % organize dimensions
        ddims       = size(Stseries{1}.(tsname).(fname{f}){1});
        ddims       = ddims(2:end);
        fdims       = size(Stseries{1}.(tsname).(fname{f}));
        dims        = [prod(fdims) ddims];
        subjdims	= [nsubj dims];
        
        [~,cidx] = ismember(cfg.plot.contr{f},cfg.plot.fact{f});
        if cidx < 1,    cidx = [];  end
        [~,ridx] = ismember(cfg.plot.factrem{f},cfg.plot.fact{f});
        if ridx < 1,    ridx = [];  end
        
        % loop over subjects
        avg = nan(subjdims);
        for s = 1:nsubj
            % get data & calculate descriptive
            dat         = Stseries{s}.(tsname).(fname{f});
            tmp         = cellfun(@(x) nanmean(x),dat,'UniformOutput',false);
            switch length(ddims)
                case 3,     avg(s,:,:,:,:)	= vertcat(tmp{:});
                case 2,     avg(s,:,:,:)	= vertcat(tmp{:});
                otherwise,  avg(s,:,:)      = vertcat(tmp{:});
            end
        end
        
        % calculate contrasts
        avg = calc_contr(avg,fdims,cidx,ddims);
        
        % remove main factor effects
        if ~isempty(ridx), error('not implemented yet'); end
        avg = remove_maineff(avg,fdims,ridx,ddims);
        
        % update dimensions
        fdims(cidx) = [];
        if isempty(fdims),      fdims = [1 1];      end
        if length(fdims) == 1,  fdims = [1 fdims];  end
        dims        = [prod(fdims) ddims];
        subjdims	= [nsubj dims];
        
        % calculate error over repetitions, or over subjects
        if nsubj == 1
            
            % initiate temporary variable
            tmp = struct;
            
            % calculate average, variance and error
            avg	= cellfun(@(x) reshape(nanmean(x),ddims),dat,'UniformOutput',false);
            var	= cellfun(@(x) reshape(nanvar(x),ddims),dat,'UniformOutput',false);
            err	= cellfun(@(x) reshape(nansem(x),ddims),dat,'UniformOutput',false);
            
            % store descriptives
            tmp.mean	= reshape(avg,fdims);
            tmp.var     = reshape(var,fdims);
            tmp.sem     = reshape(err,fdims);
            tmp.val     = reshape(dat,fdims);
            
            % store tmp in correct field
            tseries.(tsname).(fname{f}) = tmp;
            
        else
            
            % initiate temporary variable
            tmp = struct;
            
            % group descriptives
            % LENVER: this is a good point to run the approach vector
            % orientation script (below)
            % km_approachvector(avg,fdims);
            tmp.mean	= reshape(nanmean(avg,1),dims);
            tmp.var     = reshape(nanvar(avg,1),dims);
            tmp.sem     = reshape(nansem(avg,1),dims);
            tmp.subj	= reshape(avg,subjdims);
            
            % restructure
            tmp.mean    = cellfun(@(x) reshape(x,ddims),num2cell(tmp.mean,2:ndims(tmp.mean)),'UniformOutput',false);
            tmp.var     = cellfun(@(x) reshape(x,ddims),num2cell(tmp.var,2:ndims(tmp.var)),'UniformOutput',false);
            tmp.sem     = cellfun(@(x) reshape(x,ddims),num2cell(tmp.sem,2:ndims(tmp.sem)),'UniformOutput',false);
            tmp.mean    = reshape(tmp.mean,fdims);
            tmp.var     = reshape(tmp.var,fdims);
            tmp.sem     = reshape(tmp.sem,fdims);
            
            % store tmp in correct field
            tseries.(tsname).(fname{f}) = tmp;
            
            % calculate the descriptives on the difference
            if (length(fdims) > 2 && any(fdims(1:end-1) == 2)) || ...
                (length(fdims) < 3 && any(fdims == 2))
                
                % initiate temporary variable
                tmp = struct;
                
                % update the dimensions to calculate the difference
                dim         = find(fdims == 2,1,'first')+1;
                diffdims	= [nsubj fdims ddims];
                diffdims(dim) = 1;
                
                % calculate the difference of the mean and median
                avg     = reshape(diff(reshape(avg,[nsubj fdims ddims]),1,dim),diffdims);
                
                % update the dimensions
                fdims(dim-1)= 1;
                dims        = [prod(fdims) ddims];
                if length(dims) == 1
                    dims = [1 dims];
                end
                
                % group descriptives
                tmp.mean	= reshape(nanmean(avg,1),dims);
                tmp.var     = reshape(nanvar(avg,1),dims);
                tmp.sem     = reshape(nansem(avg,1),dims);
                
                % restructure
                tmp.mean    = cellfun(@(x) reshape(x,ddims),num2cell(tmp.mean,2:ndims(tmp.mean)),'UniformOutput',false);
                tmp.var     = cellfun(@(x) reshape(x,ddims),num2cell(tmp.var,2:ndims(tmp.var)),'UniformOutput',false);
                tmp.sem     = cellfun(@(x) reshape(x,ddims),num2cell(tmp.sem,2:ndims(tmp.sem)),'UniformOutput',false);
                tmp.mean    = reshape(tmp.mean,fdims);
                tmp.var     = reshape(tmp.var,fdims);
                tmp.sem     = reshape(tmp.sem,fdims);
                
                % store tmp in correct field
                tseries.(tsname).(fname{f}).diff = tmp;
                
            end
        end

    end
        
%     tseries{ts} = struct;
%     % loop over conditions
%     for c = 1:ncond
%         % keep single subject or do group analysis
%         if nsubj == 1
%             tseries.(cname{c}) = Stseries{1}.(cname{c});
%         else
%             error('work in progress');            
%             % collect subject data
%             
%         end
%     end
end
%--------------------------------------------------------------------------


% function calc_contr
%----------------------------------------
function x = calc_contr(x,dims,cidx,ddims)
if nargin < 4, ddims = []; end

% first reshape the data to full factorial size
nsubj = size(x,1);
subjdims = [nsubj dims ddims];
cdims = subjdims;
x = reshape(x,subjdims);
cidx = cidx+1;

% loop over contrast factors
for c = 1:length(cidx)
    % update the dimensions
    cdims(cidx(c)) = 1;
    % calculate the difference of the mean and median
    x = reshape(diff(x,1,cidx(c)),cdims);
    %x = reshape(-diff(x,1,cidx(c)),cdims);
end

% reshape data back to nsubj x prod(dims)
x = reshape(x,[nsubj prod(cdims(2:(1+length(dims)))) ddims]);

    
% function remove_maineff
%----------------------------------------
function x = remove_maineff(x,dims,ridx,ddims)
if nargin < 4, ddims = []; end

% return quickly if possible
if isempty(ridx),   return; end;

% first reshape the data to full factorial size
nsubj = size(x,1);
subjdims = [nsubj dims];
rdims = subjdims;
x = reshape(x,subjdims);
ridx = ridx+1;

% loop over effects to be removed
for r = 1:length(ridx)
    dimord = [1 ridx(r) setdiff(2:length(subjdims),ridx(r))];
    x = permute(x,dimord);
    for i = 1:subjdims(ridx(r))
        tmp = x(:,i,:);
        tmp = mean(tmp(:,:),2);
        for s = 1:nsubj
            x(s,i,:) = x(s,i,:)-tmp(s);
        end
    end
    x = ipermute(x,dimord);
end

% remove 2-way interaction effects
ridx = subset(ridx,2);
for r = 1:size(ridx,1)
    dimord = [1 ridx(r,:) setdiff(2:length(subjdims),ridx(r,:))];
    x = permute(x,dimord);
    for i = 1:subjdims(ridx(r,1))
        for j = 1:subjdims(ridx(r,2))
            tmp = x(:,i,j,:);
            tmp = mean(tmp(:,:,:),3);
            for s = 1:nsubj
                x(s,i,j,:) = x(s,i,j,:)-tmp(s);
            end
        end
    end
    x = ipermute(x,dimord);
end

% reshape data back to nsubj x prod(dims)
x = reshape(x,[nsubj prod(rdims(2:end))]);


% function km_approachvector
%----------------------------------------
function km_approachvector(avg,fdims)
dims = size(avg);

a = nan(dims(1:3));
for s = 1:dims(1)
    for c = 1:dims(2)
        for m = 1:dims(3)
            x = squeeze(avg(s,c,m,1,:));
            y = squeeze(avg(s,c,m,2,:));
            z = squeeze(avg(s,c,m,3,:));
            %x = interp1(0:100,x,0:4:100)';
            %y = interp1(0:100,y,0:4:100)';
            %z = interp1(0:100,z,0:4:100)';
            %figure;
            %plot3(x,y,z);
            %plot(y,z);
            %axis equal
            gra = gradient([y z]');
            gra_y = gra(1,end);
            gra_z = gra(2,end);
            a(s,c,m) = 90 - 180*atan(gra_z/gra_y)/pi + (sign(gra_y)==-1)*180;
        end
    end
end

a_avg = reshape(mean(a,1),dims(2:3));
a_sem = reshape(sem(a,1),dims(2:3));
%a_avg(:,2) = a_avg(:,2)-180;
th_avg = a_avg(:,1); th_sem = a_sem(:,1);
id_avg = a_avg(:,2); id_sem = a_sem(:,2);
th_avg = reshape(th_avg,fdims)';
th_sem = reshape(th_sem,fdims)';
id_avg = reshape(id_avg,fdims)';
id_sem = reshape(id_sem,fdims)';
figure;
errorbar(th_avg,th_sem);
figure;
errorbar(id_avg,id_sem);

figure;
errorbar(th_avg-[ 15  15; 75  75],th_sem);
figure;
errorbar(id_avg-[195 195; 255 255],id_sem);

% % test
% gra_y =  0; gra_z =  1;   % b =   0
% gra_y =  1; gra_z =  1;   % b =  45
% gra_y =  1; gra_z =  0;   % b =  90
% gra_y =  1; gra_z = -1;	% b = 135
% gra_y =  0; gra_z = -1;	% b = 180
% gra_y = -1; gra_z = -1;	% b = 225
% gra_y = -1; gra_z =  0;	% b = 270
% gra_y = -1; gra_z =  1;	% b = 315
% b = 90 - 180*atan(gra_z/gra_y)/pi + (sign(gra_y)==-1)*180