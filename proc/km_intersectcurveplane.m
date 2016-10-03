function pnts = km_intersectcurveplane(curve,plane,origin,method,allownan)
% pnts = KM_INTERSECTCURVEPLANE(curve,plane);
%
% CURVE described by [x y z] along the first dimension, points along the
% second dimension, and different curves along the third dimension (when
% applicable)
%
% PLANE described by a normal vector in [x y z]. Multiple planes are listed
% on the second.
%
% INTERPARGIN are additional input arguments for interp1
%
% PNTS are the estimated intersection points of the given curve(s) with the
% specified plane. The first dimenion of pnts lists the number of planes,
% the second dimension lists the individual curves (points per plane).

% debug flg that plots the result (very slow)
flg_debug = false;

sizc = size(curve);
sizc = [sizc ones(1,3-length(sizc))];
sizp = size(plane);
if nargin < 3
    origin = zeros(sizp);
end
sizo = size(origin);
if nargin < 4
    if sizc(2) == 1
        method = 'project';
    else
        method = {'pchip','extrap'};
    end
end
if nargin < 5
    allownan = true;
end

if ~iscell(method) && ~strcmpi(method,'project')
    method = {method};
end
if ~all([sizc(1) sizp(1) sizo(1)]==3)
    error('KM:IntersectCurvePlane:inputsize','All inputs must be specified in [x y z] on the first dimension.');
end
if sizp(2) ~= sizo(2)
    error('KM:IntersectCurvePlane:inputsize','The input arguments ''plane'' and ''orgin'' must have the same number of elements on the second dimension.');
end
if isequal(method,'project') && sizc(2) ~= sizp(2)
    error('KM:IntersectCurvePlane:inputsize','If you want to project a point on the curve on the plane, the input arguments ''curve'' and ''plane'' must have the same number of elements on the second dimension.');
end

% loop over multiple planes/origins if requested
if sizp(2) > 1
    pnts = nan([sizc(1) sizc(3) sizp(2)]);
    for rpt = 1:sizp(2)
        if isequal(method,'project')
            pnts(:,:,rpt) = km_intersectcurveplane(curve(:,rpt,:),plane(:,rpt),origin(:,rpt),method,allownan);
        else
            pnts(:,:,rpt) = km_intersectcurveplane(curve,plane(:,rpt),origin(:,rpt),method,allownan);
        end
    end
    return
end

% initialize pnts to allow early return (with NaNs)
pnts = nan(sizc([1 3]));

% scale the normvector of the plane to 1
plane = plane/norm(plane);

% plot for debug
if flg_debug
    hl = nan(2);
    figure(1); hold off;
    plot3fix(curve,'linestyle','--'); hold on;
    xyzline = [origin origin + plane/10];
    hl(1,1) = line(xyzline(1,:),xyzline(2,:),xyzline(3,:),'linewidth',4,'color','k');
    hl(1,2) = plot3fix(xyzline(:,1),'linewidth',6,'color','k','marker','o');
    axis equal;
end

% translate
tcurve = bsxfun(@minus,curve,origin);

% construct an orthogonal rotation matrix that rotates the normal vector to
% the z-axis
i = plane(1);
j = plane(2);
k = plane(3);
rotz = [i j k]'; rotz = rotz/norm(rotz);
rotx = [-k 0 i]'; rotx = rotx/norm(rotx);
roty = [i*j -(i^2+k^2) j*k]'; roty = roty/norm(roty);
rot = [rotx roty rotz];

% rotate
for i = 1:sizc(3);
    tcurve(:,:,i) = rot\tcurve(:,:,i);
end

% plot for debug
if flg_debug
    tplane = rot\plane;
    figure(2); hold off;
    plot3fix(tcurve,'linestyle','--'); hold on;
    xyzline = [zeros(3,1) tplane/10];
    hl(2,1) = line(xyzline(1,:),xyzline(2,:),xyzline(3,:),'linewidth',4,'color','k');
    hl(2,2) = plot3fix(xyzline(:,1),'linewidth',6,'color','k','marker','o');
    axis equal;
end

% continue according to chosen method
if isequal(method,'project')
    % project on plane by removing z-value
    pnts = reshape(tcurve,sizc([1 3]));
    pnts(3,:) = 0;
    
else
    % find plane cross-point and interpolate
    
    % remove point repetitions along the z-axis
    dz = diff(shiftdim(tcurve(3,:,:),1));
    for i = 1:sizc(3);
        idx = [dz(:,i); 1]~=0;
        tmp = nan(size(tcurve,1),size(tcurve,2));
        tmp(:,1:sum(idx)) = tcurve(:,idx,i);
        tcurve(:,:,i) = tmp;
    end
    
    % identify point of plane crossing closest to the origin
    idxcross = findidxcross(tcurve,allownan);
    if all(isnan(idxcross(:))), return; end
    
    % limit window to n samples in each direction
    n = 6;
    idx = max(1,min(idxcross(1,:))-n):min(sizc(2),max(idxcross(2,:))+n);
    tcurve = tcurve(:,idx,:);
    
    % find the plane crossing again
    idxcross = findidxcross(tcurve,allownan);
    if all(isnan(idxcross(:))), return; end
    
    % check for direction through the z-plane
    for i = 1:sizc(3)
        if any(isnan(idxcross(:,i))), continue; end
        z = shiftdim(tcurve(3,:,i),1);
        if z(idxcross(2,i)) < z(idxcross(1,i))
            warning('KM:IntersectCurvePlane','Not all individual trajectories are pointing in the same direction as the average.')
        end
    end
    
    % use only a section with monotoneously increasing z-values
    for i = 1:sizc(3)
        if any(isnan(idxcross(:,i)))
            tcurve(:,:,i) = NaN;
        else
            z = shiftdim(tcurve(3,:,i),1);
            dz = [sign(diff(z)); 0];
            % try to find a positive section that includes the crossing
            idx_good = logic2idx(dz>0);
            idx_good(:,2) = idx_good(:,2)+1;
            if ~isempty(idx_good)
                idx_good = idx_good(idx_good(:,1)<=idxcross(1,i) & idx_good(:,2)>=idxcross(2,i),:);
            end
            if isempty(idx_good)
                % try to find a negative section that includes the crossing
                idx_good = logic2idx(dz<0);
                idx_good(:,2) = idx_good(:,2)+1;
                if ~isempty(idx_good)
                    idx_good = idx_good(idx_good(:,1)<=idxcross(1,i) & idx_good(:,2)>=idxcross(2,i),:);
                end
            end
            if isempty(idx_good)
                % okay, it is a point
                idx_good = idxcross(:,i);
            end
            % replace identified points by NaNs
            idx_bad = setdiff(1:size(tcurve,2),idx_good(1):idx_good(2));
            tcurve(:,idx_bad,i) = NaN;
        end
    end
    
    % loop over repetitions to interpolate
    z = shiftdim(tcurve(3,:,:),1);
    pnts = zeros(sizc([1 3]));
    for i = 1:sizc(3)
        % take only non-nan values along
        idx = ~isnan(z(:,i));
        
        % return nan if no point survived
        if ~any(idx)
            pnts(1:2,i) = nan(2,1);
            continue;
        % project instead of interpolate if only a single point survived
        elseif sum(idx) == 1
            pnts(1:2,i) = tcurve(1:2,idx,i);
            continue;
        end
        
        % check if data should be extrapolated
        if sign(z(idxcross(1,i),i)) == sign(z(idxcross(2,i),i))
            % if the local gradient is not similar to the average gradient,
            % don't use interpolation, but use projection.
            gra = nangradient(tcurve(:,:,i));
            [~,j] = min(abs(z(:,i)));
            theta = vectangle(gra(:,j),[0 0 1]');
            if theta < cos(pi/6)
                z(j,i) = 0;
            end
            
            % linear interpolate curve to zero and place in new matrix
            if length(method) > 1
                tmpargin = method(2:end);
            else
                tmpargin = {};
            end
            pnts(1:2,i) = interp1(z(idx,i),tcurve(1:2,idx,i)',0,'linear',tmpargin{:})';
            
        else
            % interpolate curve to zero using chosen method and place in new matrix
            pnts(1:2,i) = interp1(z(idx,i),tcurve(1:2,idx,i)',0,method{:})';
        end
    end
    
end

% if any axis is nan, make all nan
pnts(:,any(isnan(pnts),1)) = NaN;

% plot for debug
if flg_debug
    
    figure(3); hold off;
    plot3fix(tcurve,'linestyle','--'); hold on;
    xyzline = [zeros(3,1) tplane/10];
    line(xyzline(1,:),xyzline(2,:),xyzline(3,:),'linewidth',4,'color','k');
    plot3fix(xyzline(:,1),'linewidth',6,'color','k','marker','o');
    plot3fix(pnts,'linestyle','none','linewidth',2,'marker','o','markersize',6);
    axis equal;    
    
    figure(2); hold on; delete(hl(2,:));
    line(xyzline(1,:),xyzline(2,:),xyzline(3,:),'linewidth',4,'color','r');
    plot3fix(xyzline(:,1),'linewidth',6,'color','r','marker','o');
    plot3fix(pnts,'linestyle','none','linewidth',2,'marker','o','markersize',6);
    axis equal;
    
end

% rotate back
for i = 1:sizc(3);
    tcurve(:,:,i) = rot*tcurve(:,:,i);
end
pnts = rot*pnts;

% translate back
pnts = bsxfun(@plus,pnts,origin);

% plot for debug
if flg_debug
    
    figure(1); hold on; delete(hl(1,:));
    xyzline = [origin origin + (rot*tplane)/10];
    line(xyzline(1,:),xyzline(2,:),xyzline(3,:),'linewidth',4,'color','r');
    plot3fix(xyzline(:,1),'linewidth',6,'color','r','marker','o');
    plot3fix(pnts,'linestyle','none','linewidth',2,'marker','o','markersize',6);
    axis equal;
    
end


function [idxcross,typecross] = findidxcross(dat,allownan)

x = shiftdim(dat(1,:,:),1);
y = shiftdim(dat(2,:,:),1);
z = shiftdim(dat(3,:,:),1);

idxcross = nan(2,size(z,2));
typecross = true(1,size(z,2));
for i = 1:size(z,2)
    
    % detect all crossings
    zcross = diff(sign(z(:,i)));
    
    % convert full crossing and exact zero values (that result in
    % consequetive crossings) to indeces
    idx = [logic2idx(zcross==2)' logic2idx(zcross==-2)' logic2idx(zcross==1)' logic2idx(zcross==-1)'];
    
    % update offset, and sort according to onset
    idx(2,:) = idx(2,:)+1;
    [~,k] = sort(idx(1,:));
    idx = idx(:,k);
    idx(idx>size(z,1)) = NaN;
    
    % check gradient
    dz = z(idx(2,:),i)' - z(idx(1,:),i)';
    idx = idx(:,dz>0);
        
    % select the best of all detected cross-points (or find an alternative)
    if size(idx,2) == 1
        idxcross(:,i) = idx;
    elseif size(idx,2) > 1 %if isodd(length(idx))
        % linear interpolate all cross points to z = 0;
        xy = nan(2,size(idx,2));
        for j = 1:size(idx,2)
            xy(:,j) = interp1(z(idx(:,j),i),[x(idx(:,j),i) y(idx(:,j),i)],0);
        end
        % find closest point to origin
        d = euclid(xy,zeros(2,1));
        [~,j] = min(d);
        idxcross(:,i) = idx(:,j);
    end
    
    % no true cross-point could be unambiguously detected
    if any(isnan(idxcross(:,i)))
        % if there are no cross-points, or nothing else works, take the
        % point closest to the origin
        idx = findidxcross_closest(x(:,i),y(:,i),z(:,i));
        
        % and try to see if there are neigboroughing points with the right
        % gradient
        if size(z,1) == 1 || isnan(idx)
            idxcross(:,i) = [idx; idx];
        elseif idx == 1 || isnan(z(idx-1,i))
            idxcross(1,i) = idx;
            if z(idx+1,i) > z(idx,i)
                idxcross(2,i) = idx+1;
            else
                idxcross(2,i) = idx;
            end
        elseif idx == size(z,1) || isnan(z(idx+1,i))
            idxcross(2,i) = idx;
            if z(idx,i) > z(idx-1,i)
                idxcross(1,i) = idx-1;
            else
                idxcross(1,i) = idx;
            end
        elseif z(idx-1,i) < z(idx,i)
            idxcross(:,i) = [idx-1; idx];
        elseif z(idx+1,i) > z(idx,i)
            idxcross(:,i) = [idx; idx+1];
        else
            idxcross(:,i) = [idx; idx];
        end
        
        % log type of crossing as fake
        typecross(i) = false;
    end
    
    if ~allownan && isnan(idxcross(i))
        error('KM:IntersectCurvePlane:NoCrossing','No proper crossing or close crossing could be identified');
    end    
end

% this won't work if all idxcross are NaNs
if all(isnan(idxcross))
	error('KM:IntersectCurvePlane:NoCrossing','No proper crossing or close crossing could be identified');
end


function idx = findidxcross_closest(x,y,z)
d = euclid([x y z]',zeros(3,1));
[tst,idx] = min(d);
if isnan(tst), idx = NaN; end
if sum(d==tst) > 1
    idx = find(d==tst);
    dz = nangradient(z);
    tidx = [idx(dz(idx)>0) idx];
    idx = tidx(1);        
end
    