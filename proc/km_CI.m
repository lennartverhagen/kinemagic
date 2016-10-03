function CI = km_CI(par,conf)
%--------------------------------------------------------------------------
% KM_CI calculates the 95% confidence interval of 1 dimensional data. It
% calculates the Eigen value of the data using VAR and convolves
% these with the chi-squared probabilty density function.
%
% FORMAT:   CI	= km_CI(param,conf);
% INPUT:    par - input data
%
% OUTPUT:   CI  - structure containing the
%                 some stuff
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% inspired by the function ERROR_ELLIPSE by AJ Johnson (2004).
% http:\\www.mathworkds.com\matlabcentral\fileexchange\4705-errorellipse
% pace@stny.rr.com
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% desired confidence percentile
if nargin < 2
    conf = 0.95;
end

% remove nans
par = par(~isnan(par));
par = par(:);

% the input should have at least one more elements than the degrees of
% freedom
if length(par) < 2
    % force the function to return all cirtical fields with NaNs;
    par(:) = nan;
end

% Covariance matrix
%----------------------------------------
% C = std(par)^2;
% C = cov(par);
mu = nanmean(par);
C = nanvar(par);

% Compute quantile (k) for the desired percentile (conf)
%----------------------------------------
% k is a scaling factor of the variance based on the chi-squared
% distribution. Below you will find some example:
% conf = [0.394  0.865  0.950  0.989];
% k    = [1.000  2.000  2.447  3.000];
k = sqrt(qchisq(conf,1));

CI.mean = mu;
CI.var  = C;
CI.CI   = k*C;
CI.ub   = mu + k*C;
CI.lb   = mu - k*C;

%--------------------------------------------------------------------------


%% function qchisq
%----------------------------------------
function x = qchisq(P,n)
% QCHISQ(P,N) - quantile of the chi-square distribution.
if nargin<2
  n=1;
end

s0 = P==0;
s1 = P==1;
s = P>0 & P<1;
x = 0.5*ones(size(P));
x(s0) = -inf;
x(s1) = inf;
x(~(s0|s1|s))=nan;

for ii=1:14
  dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
  x(s) = x(s)+dx;
  if all(abs(dx) < 1e-6)
    break;
  end
end
%--------------------------------------------------------------------------

%% function pchisq
%----------------------------------------
function F = pchisq(x,n)
% PCHISQ(X,N) - Probability function of the chi-square distribution.
if nargin<2
  n=1;
end
F=zeros(size(x));

if rem(n,2) == 0
  s = x>0;
  k = 0;
  for jj = 0:n/2-1;
    k = k + (x(s)/2).^jj/factorial(jj);
  end
  F(s) = 1-exp(-x(s)/2).*k;
else
  for ii=1:numel(x)
    if x(ii) > 0
      F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
    else
      F(ii) = 0;
    end
  end
end
%--------------------------------------------------------------------------

%% function dchisq
%----------------------------------------
function f = dchisq(x,n)
% DCHISQ(X,N) - Density function of the chi-square distribution.
if nargin<2
  n=1;
end
f=zeros(size(x));
s = x>=0;
f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));
%--------------------------------------------------------------------------
