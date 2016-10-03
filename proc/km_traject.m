function [cfg,data] = km_traject(cfg,data)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2012-11-01
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end

% overwrite previous trajectory analyses
data.traject = [];

% retrieve options from cfg
flg_ellipsoid = ismember('ellipsoid',cfg.traject.type);
flg_project = ismember('project',cfg.traject.type);
flg_intersect = ismember('intersect',cfg.traject.type);
% number of points for each trial trajectory
npnts_rept = cfg.traject.npnts_rept;
% number of points for the average trajectory
npnts_avg = cfg.traject.npnts_avg;
% calculate points around a circle
if ~isfalse(cfg.traject.calc_pnts)
    q = linspace(0,2*pi,cfg.traject.npnts_circle)';
    pnts = zeros(length(q),3);
    pnts(:,[2 3]) = [cos(q) sin(q)];
end

% get trial variables and variable names
if isfield(cfg,'trl')
    trl = cfg.trl;
else
    error('no trl provided');
end
if isfield(cfg,'vars')
    vars = cfg.vars;
else
    error('no vars provided');
end

% get movement
if isfield(cfg.traject,'movname');
    movname = cfg.traject.movname;
else
    movname = 'movment';
end
if isfield(cfg,movname) && ~isstruct(cfg.(movname)) && ~isempty(cfg.(movname))
    mov = cfg.(movname);
elseif isfield(cfg,'movdef') && isfield(cfg.movdef,movname)
    mov = cfg.movdef.(movname);
elseif isfield(data,movname)
    mov = data.(movname);
else
    error('no %s provided',movname);
end

% make sure that all trials have the same number of movements
mov = km_setequalnmov(mov);

% get position data
tcfg = [];
tcfg.marker = cfg.traject.marker;
[pos,ax,sns] = km_getmovdat(tcfg,data,mov,'pos');
cfg.traject.marker = sns;
nsns = length(cfg.traject.marker);
if nsns == 1, pos = cellfun(@(x) shiftdim(x,-1),pos,'UniformOutput',false); end

% get velocity data
tcfg = [];
tcfg.marker = cfg.traject.marker;
tcfg.axis = 'xyz';
[vel,ax,sns] = km_getmovdat(tcfg,data,mov,'vel');
if nsns == 1, vel = cellfun(@(x) shiftdim(x,-1),vel,'UniformOutput',false); end
dist = cellfun(@(x) cumsum(x,3),vel,'UniformOutput',false);

% reject trials if requested
if ~isfalse(cfg.traject.trlrej)
    idxr = km_getidxr(cfg,cfg.traject);
    
    % remove trials from list
    pos = pos(~idxr.reject);
    dist = dist(~idxr.reject);
    trl = trl(~idxr.reject,:);
end

% interpolate pos to npnts_trl equidistant points
interppos = cellfun(@(x,y) interp2xi(x,y,npnts_rept),dist,pos,'UniformOutput',false);

% cell 2 matrix
interppos = km_datcell2mat(interppos);

% store interppos and trl in a temporary field
if cfg.traject.sesscounter < 2 && isfield(cfg,'tmp')
    cfg = rmfield(cfg,'tmp');
end
cfg.tmp(cfg.traject.sesscounter).interppos = interppos;
cfg.tmp(cfg.traject.sesscounter).trl = trl;

% only continue processing of enough sessions have been collected
if cfg.traject.sesscounter < cfg.traject.nsesscomb
    % update session counter
    cfg.traject.sesscounter = cfg.traject.sesscounter + 1;
    return;
end

% collect interppos and trl from temporary field and remove tmp
interppos = vertcat(cfg.tmp(:).interppos);
trl = vertcat(cfg.tmp(:).trl);
cfg = rmfield(cfg,'tmp');

% loop over factors
for f = 1:length(cfg.traject.fact)
    % get factor(s)
    fact = cfg.traject.fact{f};
    fname = cfg.traject.fname{f};
    [~,iv] = ismember(fact,vars);
    
    % get levels
    factlvl = cell(1,length(fact));
    for ff = 1:length(fact)
        factlvl{ff} = cfg.traject.lvl.(fact{ff});
    end
    
    % combine levels of all factors
    lvl = arraycomb(factlvl{:});
    nlvl = size(lvl,1);
    
    % determine dimensions
    dims = cellfun(@(x) length(x),factlvl);
    if length(dims) == 1
        dims = [dims 1];
    end
    
    % loop over markers
    trlsel = cell(1,nlvl);
    for sr = 1:nsns
        srname = cfg.traject.marker{sr};
        
        % initialize CI structures (store levels in cells)
        tmp_avgpos	= cell(size(lvl,1),1);
        tmp_avggra	= cell(size(lvl,1),1);
        if flg_project,     CIellipsoid	= cell(size(lvl,1),1);  end
        if flg_project,     CIproject   = cell(size(lvl,1),1);  end
        if flg_intersect,   CIintersect	= cell(size(lvl,1),1);  end
        
        % devide end positions over factor x level combinations and calculate
        % multi-dimensional confidence intervals
        for i = 1:nlvl
            idx = all(trl(:,iv) == repmat(lvl(i,:),size(trl,1),1),2);
            
            % if not enough trials are selected, just go on...
            nrpt = sum(idx);
            if nrpt < min(cfg.traject.ntrl.(fname))
                continue
            elseif sr == 1
                if nrpt <= max(cfg.traject.ntrl.(fname))
                    trlsel{i} = idx;
                else
                    tmpidx = find(idx);
                    tmpidx = tmpidx(randperm(length(tmpidx)));
                    n = max(cfg.traject.ntrl.(fname));
                    trlsel{i} = false(size(idx));
                    trlsel{i}(tmpidx(1:n)) = true;
                end
            elseif sr > 1
                idx = trlsel{i};
            end
            nrpt = sum(idx);
            
            % select data and reshape for processing
            tmp = interppos(idx,sr,:,:);
            siz = size(tmp);
            tmp = reshape(tmp,siz([1 3 4]));
            
            % iteratively create spatial average to converge on solution
            % FIXME: make avgpos converge (and test for it
            avgpos = shiftdim(nanmean(tmp,1));
            dconv = nan(1,8);
            for c = 1:cfg.traject.avgcalc_iter
                
                % calculate trajectory distance from pos
                x = cumsum([0 sqrt(sum(diff(avgpos,1,2).^2))]);
                % remove equal points (which might arise by averaging)
                d = [true (diff(x)>0)];
                if ~all(d), x = x(d); avgpos = avgpos(:,d); end
                    
                % resample avgpos to n equidistance positions
                avgpos = interp2xi(x,avgpos,npnts_avg);
                
                % transform tmp to nearest trajectory to average
                tmppos = nan([nrpt size(avgpos)]);
                for p = 1:size(avgpos,2)
                    
                    % calculate distance from all points to the average
                    d = sqrt(sum(bsxfun(@minus,shiftdim(tmp,1),avgpos(:,p)).^2,1));
                    % find closest points
                    [d,idx] = min(shiftdim(d,1),[],1);
                    
                    % extract closest points
                    C = cellfun(@(x) tmp(x,:,idx(x)),num2cell(1:nrpt),'UniformOutput',false);
                    C = vertcat(C{:});
                    if sum(~any(isnan(C),2))<4
                        A = 1;
                    end
                    
                    % replace far away points with NaN
                    if ~ischar(cfg.traject.avgcalc_outlier)
                        C(d>cfg.traject.avgcalc_outlier,:) = NaN;
                    elseif strcmpi(cfg.traject.avgcalc_outlier,'best')
                        tmpCI = km_mdCI(C,'medium');
                        C(tmpCI.xyz.d>3,:) = NaN;
                    elseif strcmpi(cfg.traject.avgcalc_outlier,'dist')
                        C(d>0.05,:) = NaN;
                    end
                    
                    % store
                    tmppos(:,:,p) = C;
                    
                end
                
                % remove individual trajectories with many NaNs completely
                nnan = sum(isnan(squeeze(tmppos(:,1,:))),2);
                %tmppos(nnan>npnts_avg/8,:,:) = NaN;
                tmppos(nnan>0,:,:) = NaN;
                
                % calculate a new average
                newavgpos = shiftdim(nanmean(tmppos,1));
                dconv(c) = mean(euclid(avgpos,newavgpos));
                avgpos = newavgpos;
                
            end
            
            % the average trajectory can be regarded as a normal vector of
            % an orthogonal plane. Find this normal vector by calculating
            % direction of trajectory as gradient in xyz
            avggra = gradient(avgpos);
            
            % store avgpos and gra in temporary cell
            tmp_avgpos{i} = avgpos;
            tmp_avggra{i} = avggra;
            
            if flg_ellipsoid
                % calculate the 3D confidence ellipsoid of the associated points
                CIellipsoid{i} = km_mdCI(tmppos,'minimal');
            end
            
            if flg_project
                % find the projection of the associated points on that plane
                pntsproject = km_intersectcurveplane(permute(tmppos,[2 3 1]),avggra,avgpos,'project');
                % now calculate the confidence ellipses of these points
                CIproject{i} = km_mdCI(permute(pntsproject,[2 1 3]),'minimal','nonposdef');
            end
            
            if flg_intersect
                % find the points along the trajectories that intersect with
                % the plane (for multiple crossings only the closest are taken)
                pntsintersect = km_intersectcurveplane(permute(tmppos,[2 3 1]),avggra,avgpos,{'pchip','extrap'});
                % now calculate the confidence ellipses of these points
                CIintersect{i} = km_mdCI(permute(pntsintersect,[2 1 3]),'minimal','nonposdef');
            end
            
        end
        
        % store average trajectory
        avgpos = permute(cat(3,tmp_avgpos{:}),[1 3 2]);
        data.traject.(fname).(srname).pos = avgpos;
        avggra = permute(cat(3,tmp_avggra{:}),[1 3 2]);
        data.traject.(fname).(srname).gra = avggra;
        
        if flg_ellipsoid
            % collect results in concatenated matrices
            CIellipsoid = collect_CIfields(CIellipsoid,'xyz','volume');
            if ~isfalse(cfg.traject.calc_pnts)
                % align vectors and calculate points to plot ellipse
                CIellipsoid = process_CI(CIellipsoid,avgpos,avggra,pnts);
            end
            data.traject.(fname).(srname).CIellipsoid = CIellipsoid;
        end
        if flg_project
            % collect results in concatenated matrices
            CIproject = collect_CIfields(CIproject,'xyz');
            if ~isfalse(cfg.traject.calc_pnts)
                % align vectors and calculate points to plot ellipse
                CIproject = process_CI(CIproject,avgpos,avggra,pnts);
            end
            data.traject.(fname).(srname).CIproject = CIproject;
        end
        if flg_intersect
            % collect results in concatenated matrices
            CIintersect = collect_CIfields(CIintersect,'xyz');
            if ~isfalse(cfg.traject.calc_pnts)
                % align vectors and calculate points to plot ellipse
                CIintersect = process_CI(CIintersect,avgpos,avggra,pnts);
            end
            data.traject.(fname).(srname).CIintersect = CIintersect;
        end
        
    end
    
    % store levels in configuration structure
    data.traject.(fname).lvl = lvl;
    
end

% evaluate experimentally specific function
if isfield(cfg.traject,'expfun') && ~isfalse(cfg.traject.expfun)
    data.traject = feval(cfg.traject.expfun,data.traject);
end


% update dataset
if isempty(regexpi(cfg.dataset,'_ana$','once'))
    if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
    cfg.dataset = [cfg.dataset '_ANA'];
end

% update processing directory
cfg = km_setcfg(cfg,'dirana');
cfg.dir.proc = cfg.dir.ana;
%--------------------------------------------------------------------------


% another idea...
% convolve with 3D gaussian and use cutoff to find cloud
% > isosurface, isonormal, smooth3

% another idea...
% calculate envelope of confidence ellipsoids for each distance point
% > use patch


function dati = interp2xi(tim,dat,timi,method)

% check input
if nargin < 4
    method = 'pchip';
end

% check size of original data
sizdat = size(dat);
if sizdat(end) < 4
    error('number of samples is too small for sensible interpolation');
end

% check size of original time
siztim = size(tim);
if sum(siztim(1:end-1)>1)>1
    error('dimensions of tim not unambiguous');
end
tim = reshape(tim,[prod(siztim(1:end-1)) siztim(end)]);
if numel(tim) == length(tim), tim = tim(:)'; end
siztim = size(tim);
if siztim(1)>sizdat(1)
    error('dimensions of tim do not match dimensions of dat');
end

% check size of interpolated time
if length(timi) == 1
    ntim = timi;
else
    siztimi = size(timi);
    if sum(siztimi(1:end-1)>1)>1
        error('dimensions of timi not unambiguous');
    end
    timi = reshape(timi,[prod(siztimi(1:end-1)) siztimi(end)]);
    if numel(timi) == length(timi), timi = timi(:)'; end
    siztimi = size(timi);
    if siztim(1:end-1) ~= siztimi(1:end-1)
        error('dimensions of timi do not match dimensions of tim');
    end
    ntim = siztimi(end);
end

% loop over first dimension of tim
dati = nan([sizdat(1:end-1) ntim]);
for i = 1:siztim(1)
    % do not allow NaNs
    if any(isnan(dat(i,:))), continue; end
        
    % get original time array
    x = tim(i,:);
    
    % get new time array
    if length(timi) == 1
        xi = linspace(min(x),max(x),timi);
    else
        xi = timi(i,:);
    end
    
    % if loop over tim is needed
    if siztim(1)>1 || length(sizdat)>2
        % get original data matrix
        y = nan(sizdat(2:end));
        y(:) = dat(i,:);

        % interpolate new data matrix and place in new matrix
        yi = interp1(x,y',xi,method,'extrap')';
        dati(i,:) = yi(:);
    
    % if no loop is needed
    else
        dati = interp1(x,dat',xi,method,'extrap')';
    end
    
end


function CI = process_CI(CI,pos,gra,pnts)
% align all eigenvariates iteratively to the middle point
[CI.vec, CI.val] = km_traject_vectors_align(CI.vec,CI.val,gra);
% calculate points by convolving a circle with the vector and values
CI.pnts = km_traject_pnts_conv(CI.vec,CI.val,pos,CI.k,pnts);
% align points
CI.pnts = km_traject_pnts_align(CI.pnts,pos);


function Sout = collect_CIfields(S,axfld,flg)
% restructure and collect fields from S
% and store in concatenated matrices

% check input
if nargin < 2
    axfld = 'xyz';
end
if nargin < 3
    flg = 'area';
end

nlvl = length(S);
nrpt = length(S{1});
k = nan(nlvl,nrpt);
vec = nan([size(S{1}(1).(axfld).eigvec) nlvl nrpt]);
val = nan([size(S{1}(1).(axfld).eigval) nlvl nrpt]);
for i = 1:nlvl
    for j = 1:nrpt
        k(i,j) = S{i}(j).(axfld).k;
        vec(:,:,i,j) = S{i}(j).(axfld).eigvec;
        val(:,:,i,j) = S{i}(j).(axfld).eigval;
    end
end

% store in output structure
Sout        = [];
Sout.k      = k;
Sout.vec	= vec;
Sout.val	= val;
if strcmpi(flg,'volume')
    [Sout.area Sout.volume] = km_traject_areavolume(val,k,flg);
else
    Sout.area = km_traject_areavolume(val,k,flg);
end
