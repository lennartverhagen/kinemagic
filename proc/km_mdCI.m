function CI = km_mdCI(pos,varargin)
%--------------------------------------------------------------------------
% KM_MDCI takes multi-dimensional (2D/3D) position data (pos) and
% calculates the 2D 95% confidence ellipses of all orientations (e.g. XY,
% YZ, ZX) and the 3D 95% confidence ellipsoid for the 3D data. It
% calculates the eigen vectors and eigen values of the data and convolves
% these with the chi-squared probabilty density function. It reports
% on all these CI ellipses and ellipsoids in a structure called CI.
%
% FORMAT:   CI      = km_mdCI(pos,conf);
% INPUT:    pos     - input data with the different points on the first
%                     (vertical) dimension of the matrix:
%                     [ X1 Y1 Z1
%                       X2 Y2 Z2
%                       X3 Y3 Z3
%                       X4 Y4 Z4 ];
%           conf    - the desired confidence interval (1-alpha),
%                     default: 0.95
%           k       - scaling factor k, or scaling setting 'equal' or
%                     'unequal' (default) for all subsets of the dimensions
%           kdist   - distribution to use to determine the k: 'chisq'
%                     (default) or 'f'
%
% OUTPUT:   CI  - structure containing the
%           conf    - confidence interval
%           kdist   - distribution used to determine k
%           k       - the k scaling factor (nb: for all dimension levels,
%                     not for all dimensions individually)
%           mu      - the mean of the data
%           struct  - structure containing the
%               1D  (x, y, z)
%               2D  (xy, yz, zx)
%               3D  (xyz)
%                   - eigen vectors and values
%                   - k scaling factor
%                   - CI interval projected on original axes (x, y, z)
%                   - CI interval projected on principal axes (a, b, c)
%                       1D: length of the interval (a)
%                       2D: area of the interval ellipse (pi*ab)
%                       3D: volume of the interval ellipsoid ((3/4)*pi*abc)
%                   - logical indexing of data points falling outside the
%                     confidence interval
%                   - origin of the CI interval
%                   - points to plot
%                       1D: 2-by-3 points
%                       2D: 100-by-3 points
%                       3D: 20-by-20-by-3 points
%
% Chi-square distribution scaling taken from the function ERROR_ELLIPSE by
% AJ Johnson (2004).
% http://www.mathworkds.com/matlabcentral/fileexchange/4705
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010-2011, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2011-09-01
%--------------------------------------------------------------------------


%% Check input
%----------------------------------------
narginchk(1,6);
% set defaults
conf    = 0.95;
k       = 'equal';
kdist	= 'chisq';
output	= 'all';
warn_nonposdef = 'posdef';
n2Dpnts = 'n2Dpnts100';
n3Dpnts = 'n3Dpnts20';
% search for conf
idx = cellfun(@(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1,varargin);
if any(idx)
  conf = varargin{idx}; varargin = varargin(~idx);
end
% search for k as a numeric value
idx = cellfun(@(x) isnumeric(x) && any(x > 1),varargin);
if any(idx)
  k = varargin{idx}; varargin = varargin(~idx);
  conf = [];
  kdist = 'n/a';
else
  % search for k as 'equal' or 'unequal'
  idx = cellfun(@(x) ischar(x) && ~isempty(regexp(x,'^(un)?equal$','once')),varargin);
  if any(idx)
    k = varargin{idx}; varargin = varargin(~idx);
  end
  % search for kdist
  idx = cellfun(@(x) ischar(x) && ~isempty(regexp(x,'^chisq$|^f$','once')),varargin);
  if any(idx)
    kdist = varargin{idx}; varargin = varargin(~idx);
  end
  % search for output
  idx = cellfun(@(x) ischar(x) && ~isempty(regexp(x,'^all$|^medium$|^minimal$','once')),varargin);
  if any(idx)
    output = varargin{idx}; varargin = varargin(~idx);
  end
  % search for warn_nonposdef
  idx = cellfun(@(x) ischar(x) && ~isempty(regexp(x,'^posdef$|^nonposdef$','once')),varargin);
  if any(idx)
    warn_nonposdef = varargin{idx}; varargin = varargin(~idx);
  end
  % search for n2Dpnts
  idx = cellfun(@(x) ischar(x) && ~isempty(regexp(x,'^n2Dpnts','once')),varargin);
  if any(idx)
    n2Dpnts = varargin{idx}; varargin = varargin(~idx);
  end
  % search for n3Dpnts
  idx = cellfun(@(x) ischar(x) && ~isempty(regexp(x,'^n3Dpnts','once')),varargin);
  if any(idx)
    n3Dpnts = varargin{idx}; varargin = varargin(~idx);
  end
end
% warn if not all input arguments were used
if ~isempty(varargin)
  A = cellcont2str(varargin);
  B = repmat({' '},size(A)); B{end} = [];
  str = [A; B]; str = [str{:}];
  warning('KM:km_mdCI:NotAllInputUsed','Not all input arguments were used: %s',str);
end

% number of axes of position data
dims = size(pos);
naxes = size(pos,2);
if naxes > 3
  error('KM:km_mdCI:TooManyDimensions','Your data has more than 3 dimensions, although this should in principle not be a problem for the core mathematics, this function was designed to go up to only 3 dimensions, and won''t work properly with more.');
end

% allow repetitive calls to this function over the third dimension
if length(dims)>2
  for i = 1:prod(dims(3:end))
    CI(i) = km_mdCI(pos(:,:,i),conf,k,kdist,output,warn_nonposdef,n2Dpnts,n3Dpnts);
  end
  CI = reshape(CI,[dims(3:end) 1]);
  return;
end

% remove nans
idx_good = ~any(isnan(pos),2);
pos = pos(idx_good,:);

% the input should have at least one more element than the number of axes
if size(pos,1) < naxes+1
  % force the function to return all critical fields with NaNs;
  pos(:) = nan;
end

% number of points on the 2D ellipse and 3D ellipsoid
if ischar(n2Dpnts)
  n2Dpnts = str2double(regexp(n2Dpnts,'\d*$','once','match'));
end
if isempty(n3Dpnts) || ~isnumeric(n2Dpnts)
  n2Dpnts = 100;
end
if ischar(n3Dpnts)
  n3Dpnts = str2double(regexp(n3Dpnts,'\d*$','once','match'));
end
if isempty(n3Dpnts) || ~isnumeric(n3Dpnts)
  n3Dpnts = 20;
end


%% Covariance matrix
%----------------------------------------
% There are multiple ways to estimate the covariance matrix of a set of
% multivariate data points. Matlab uses a standard algortithm, but this
% approach is very unstable when using too few datapoints. I feel that this
% is the most vulnarable part of this code. It could be improved by
% implementing a penalized maximum likelyhood estimation, a Bayesian
% approach or including a shrinkage factor in the estimation.
% sigma_xyz = std(pos,1);
% rho = corrcoef(pos);
% sigma = [ sigma_x^2 rho*sigma_x*sigma_y rho*sigma_x*sigma_z
%           rho*sigma_x*sigma_y sigma_y^2 rho*sigma_y*sigma_z
%           rho*sigma_x*sigma_z rho*sigma_y*sigma_z sigma_z^2];
C = cov(pos);

% Make sure the matrix has positive eigenvalues,
% otherwise it's not a valid covariance matrix!
if strcmpi(warn_nonposdef,'posdef') && (any(~isfinite(C(:))) || any(eig(C)<=0))
  warning('KM:km_mdCI:NonPositiveDefinite','Probably no data points were entered, or they were all NaN. Some eigenvalues of the covariance matrix are non-positive, indicating that the covariance matrix is not positive-definite (which is required). Maybe this is due to the matrix being close to singlar which causes the eig algorithm to produce negative or zero values. Please check your data, but for now, as a dirty hack, only the absolute values are used.');
  %error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
end

%% Obtain/Compute quantile (k) for the desired percentile (conf)
%----------------------------------------
% To determine the confidance interval (for a specific alpha) the
% eigenvalues of the covariance matrix need to be scaled on the basis of an
% expected/appropriate statistical distribution (scaling factor k). One
% could use a chi-squared distribution, or an F-distribution (e.g.
% McIntyre, JNS, 1998), or as in the corrcoef function of Matlab use an
% asymptotic normal distribution: z = 0.5 * log((1+rv)./(1-rv)) with an
% approximate variance equal to 1/(n-3). Again, as with the estimation of
% the covariance matrix, also the scaling of the eigenvalues is only valid
% when using a large number of datapoints. Here I have chosen to apply a
% chi-squared distribution by default. You can also choose the
% F-distribution if you have a statistics toolbox licence, but in
% simulations I have found the chi-squared to be more accurate and more
% robust compared to the F-distribution.
% Here are some examples of the chi-squared:
% conf = [0.394  0.865  0.950  0.989];
% k    = [1.000 2.000  2.447  3.000];
%
% n is the number of samples (denominator degrees of freedom)
% r is the number of dimensions (numerator degrees of freedom)
[n r] = size(pos);
if isnumeric(k)
  if length(k) == 1
    k = ones(1,r) * k;
  elseif length(k) ~= r;
    error('KM:km_mdCI:kVector','if input k is specified as a vector, its length must match the number of dimensions');
  end
elseif ischar(k)
  switch k
    case 'equal'
      if strcmpi(kdist,'chisq')
        k = ones(1,r) * sqrt(qchisq(conf,r));
      elseif strcmpi(kdist,'f')
        k = ones(1,r) * finv(conf,r,n-r+1);
      end
    case 'unequal'
      k = nan(1,r);
      for i = 1:r
        if strcmpi(kdist,'chisq')
          k(i) = sqrt(qchisq(conf,i));
        elseif strcmpi(kdist,'f')
          k(i) = finv(conf,i,n-i+1);
        end
      end
  end
end

% store confidence interval, distribution used and k scaling factor
if strcmp(output,'all')
  CI.conf	= conf;
  CI.kdist = kdist;
  CI.k = k;
  CI.mu = mean(pos,1);
end

% Compute confidence interval, ellipse and ellipsoid
%----------------------------------------
% store the x, y and z 95%-CI eigen-values
if strcmp(output,'all')
  CI.x = store_eigs(pos,C,k(1),1,'x',idx_good,[],output);
  if naxes > 1
    CI.y    = store_eigs(pos,C,k(1),2,'y',idx_good,[],output);
    CI.xy   = store_eigs(pos,C,k(2),[1 2],'xy',idx_good,n2Dpnts,output);
  end
  if naxes > 2
    CI.z    = store_eigs(pos,C,k(1),3,'z',idx_good,[],output);
    CI.yz   = store_eigs(pos,C,k(2),[2 3],'yz',idx_good,n2Dpnts,output);
    CI.zx   = store_eigs(pos,C,k(2),[3 1],'zx',idx_good,n2Dpnts,output);
    CI.xyz	= store_eigs(pos,C,k(3),1:3,'xyz',idx_good,n3Dpnts,output);
  end
else
  if naxes == 1
    CI.x = store_eigs(pos,C,k(1),1,'x',[],idx_good,output);
  elseif naxes == 2
    CI.xy   = store_eigs(pos,C,k(2),[1 2],'xy',n2Dpnts,idx_good,output);
  else
    CI.xyz	= store_eigs(pos,C,k(3),1:3,'xyz',n3Dpnts,idx_good,output);
  end
end
%--------------------------------------------------------------------------


%% function store_eigs
%----------------------------------------
function CI = store_eigs(pos,C,k,idx_ax,axisname,idx_rpt,n,output)
% sort input
if nargin < 8,      output = 'all';                 end
if nargin < 7,      n = [];                         end
if nargin < 6,      idx_rpt = ones(size(pos,1),1);	end
if nargin < 5,      axisname = 'xyz';               end
if nargin < 4,      idx_ax = 1:3;                   end
if islogical(idx_ax),  idx_ax = find(idx_ax);       end
if length(axisname) > length(idx_ax)
  axisname = axisname(idx_ax);
elseif length(idx_ax) > length(axisname)
  error('the axis indeces [%s] do not match the axis names [%s]',num2str(idx_ax),axisname);
end

% select (a subset of) the data and calculate the means
pos = pos(:,idx_ax);
mu = mean(pos,1);

% select (a subset of) the covariance matrix
C = C(idx_ax,idx_ax);
naxes = length(idx_ax);

% calculate the eigen vectors and values
[eigvec,eigval] = calc_eigs(C);

% store eigen vectors and values
CI.eigvec	= eigvec;
CI.eigval	= eigval;
CI.k        = k;

% return if minimal output is requested
if strcmp(output,'minimal')
  return;
end

% get the projections of the eigenvalues on the x,y and z axis
% this part of the code should be checked and rechecked
xyzeigs = k * diag(eigvec * sqrt(eigval) * eigvec');
%xyzeigs = k * diag(sqrt(eigval) * eigvec');

% convolve eigen values with the chi-square scaling factor
abceigs = k * sqrt(diag(eigval));

% store the x, y and z 95%-CI eigen-values
CI.x = 0; CI.y = 0; CI.z = 0;
for i = 1:naxes
  CI.(axisname(i)) = xyzeigs(i);
end

% store multi-dimensional 95%-CI eigen-values from the largest to the smallest
abceigs = sort(abceigs,'descend');
CI.a = abceigs(1);
if naxes == 2
  CI.b = abceigs(2);
  CI.area = pi * prod(abceigs);
elseif naxes == 3
  CI.b = abceigs(2);
  CI.c = abceigs(3);
  CI.volume = (4/3) * pi * prod(abceigs);
end


%% decompose ellipsoid to test whether points are located inside or outside
% the code below works for 3 dimensions (ellipsoid), but is backwards
% compatible for ellipses (2D) and lines (1D)
%----------------------------------------
% subtract origin from data
d = bsxfun(@minus,pos,mu);
% devide by scaling factor k
d = d ./ k;
% find scaling and rotation matrix: product of sqrt(eigval) and eigvec'
X = (sqrt(eigval)*eigvec');
% if X-matrix is singular (eigvalue of zero), ignore spurious factor
X = X(~all(X==0,2),:);
% devide by transformation matrix (only if no NaNs present in X)
if ~any(isnan(X(:)))
  d = d/X;
end
% now we have rotated and scaled the data back to fit in a known ellipsoid
% with the radii a=b=c=1 (a sphere) which follows x.^2 + y.^2 + z.^2 = 1
d = sum(d.^2,2);
CI.d = nan(size(idx_rpt));
CI.idx_out = nan(size(idx_rpt));
CI.d(idx_rpt) = d;
CI.idx_out(idx_rpt) = d > 1;

% return if medium output is requested
if strcmp(output,'medium')
  return;
end

%% store points for interval, ellipse and ellipsoid plotting
%----------------------------------------
% store the x, y and z origins
CI.x0 = 0; CI.y0 = 0; CI.z0 = 0;
for i = 1:naxes
  CI.([axisname(i) '0']) = mu(i);
end

if naxes == 1
  % create the two extreme points and convolve with the eigen vectors and
  % values
  p = [-1;1] * sqrt(eigval) * eigvec';
  pnts.x = zeros(2,1); pnts.y = zeros(2,1); pnts.z = zeros(2,1);
  pnts.(axisname) = p;
  
  % convolve with the chi-square-distribution
  pnts.x = CI.x0 + k*pnts.x;
  pnts.y = CI.y0 + k*pnts.y;
  pnts.z = CI.z0 + k*pnts.z;
  
  % points to plot using plot/plot3
  CI.xyz = [pnts.x pnts.y pnts.z];
  
elseif naxes == 2
  % Points around a circle
  if isempty(n), n = 100; end
  p = linspace(0,2*pi,n)';
  
  % create a standard ellipse (circle) and convolve with the eigen
  % vectors and values
  p = [cos(p) sin(p)] * sqrt(eigval) * eigvec';
  pnts.x = zeros(n,1); pnts.y = zeros(n,1); pnts.z = zeros(n,1);
  for i = 1:naxes
    pnts.(axisname(i)) = p(:,i);
  end
  
  % convolve with the chi-square-distribution scaling
  pnts.x = CI.x0 + k*pnts.x;
  pnts.y = CI.y0 + k*pnts.y;
  pnts.z = CI.z0 + k*pnts.z;
  
  % points to plot using plot/plot3
  CI.xyz = [pnts.x pnts.y pnts.z];
  
elseif naxes == 3
  % create a standard ellipsoid and convolve with the eigen vectors
  if isempty(n), n = 20; end
  [x,y,z] = ellipsoid(0,0,0,1,1,1,n-1);
  xyz = [x(:) y(:) z(:)] * sqrt(eigval) * eigvec';
  
  % convolve with the chi-square-distribution scaling
  x(:) = CI.x0 + k*xyz(:,1);
  y(:) = CI.y0 + k*xyz(:,2);
  z(:) = CI.z0 + k*xyz(:,3);
  
  % points to plot using mesh/surf
  CI.xyz = nan([size(x) 3]);
  CI.xyz(:,:,1) = x;
  CI.xyz(:,:,2) = y;
  CI.xyz(:,:,3) = z;
end
%--------------------------------------------------------------------------


%% function calc_eigs
%----------------------------------------
function [eigvec,eigval] = calc_eigs(C)
% compute eigen-vectors and -values
if ~any(isnan(C(:)));
  [eigvec,eigval] = eig(C);
else
  eigvec = nan(size(C));
  eigval = nan(size(C));
end

% the eig algorithm can return rounding errors or even negative eigen
% values when the matrix is (close to) singular. The matrix will for
% example be singular when all points specified in 3D space lie on the same
% plane. Because this can be intentional, here is a dirty hack to ignore
% likely rounding errors and make all values positive.
eigval(abs(eigval)<=eps) = 0;
eigval = abs(eigval);
%--------------------------------------------------------------------------