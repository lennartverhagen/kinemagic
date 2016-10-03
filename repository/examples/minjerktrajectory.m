%% Minimum-jerk-trajectory
% Flash T, Hogan N (1985) The coordination of arm movements: an
% experimentally confirmed mathematical model. J neurosci 5:1688-1703
%
% See also http://www.shadmehrlab.org/book/minimum_jerk/minimumjerk.html
%
% Finding a trajectory for a fixed distance and time can be
% operationalized as finding the trajectory with the minimal summed squared
% jerk. This can be quantified as:
%
% x(t) = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5;
%
% with derivatives
%
% v(t) = c1 + 2*c2*t + 3*c3*t^2 + 4*c4*t^3 + 5*c5*t^4;
% a(t) = 2*c2 + 6*c3*t + 12*c4*t^2 + 20*c5*t^3;
%
% Smeets and Brenner (1999) proposed to regard grasping movements as the
% minimum-jerk-trajectory of two digits moving independently with the
% following contraints:
%
% x(0) = 0;     x(mt) = d;
% v(0) = 0;     v(mt) = 0;
% a(0) = 0;     a(mt) = ap/mt^2
%
% legend:
% x = position; v = velocity; a = acceleration
% d = distance; mt = movement time; ap = "approach parameter"
%
% this results in the following parameter values
%
% c0 = 0;
% c1 = 0;
% c2 = 0;
% c3 = (10*d+ap) / (2*mt^3);
% c4 = (-15*d-ap) / (mt^4);
% c5 = (12*d+ap) / (2*mt^5);
%
% for ease of analysis, a normalized time ratio is introduced
% tr = t/mt;
%
% the new constrained and substituted equation now yields
%
% x(tr) = (0.5*ap*(tr-1)^2 + d*(6*tr^2 - 15*tr + 10)) * tr^3;


%% Example
% A disc at distance d with radius r is grasped at two opposing points,
% located at an angle phi from the x-axis

% settings
d = 0.4;
r = 0.03;
phi = -30*(pi/180);
ap = 1;
tr = 0:0.01:1;

% goal
xg = r.*sin(0:0.01:2*pi);
yg = r.*cos(0:0.01:2*pi) + d;

% thumb
xt = -cos(phi) .* (0.5.*ap.*(tr-1).^2 + r.*(6.*tr.^2 - 15.*tr + 10)) .* tr.^3;
yt = (-0.5.*ap.*sin(phi).*(tr-1).^2 + (d-r.*sin(phi)) .* (6.*tr.^2 - 15.*tr + 10)) .* tr.^3;

% index finger
xi = cos(phi) .* (0.5.*ap.*(tr-1).^2 + r.*(6.*tr.^2 - 15.*tr + 10)) .* tr.^3;
yi = (0.5.*ap.*sin(phi).*(tr-1).^2 + (d+r.*sin(phi)) .* (6.*tr.^2 - 15.*tr + 10)) .* tr.^3;

% calculate grip aperture in Euclidian distance (using Pythagoras)
apt = [xi-xt; yi-yt];
apt = squeeze(sqrt(sum(apt.^2,1)));

% plot
figure;
subplot(2,1,1); hold off;
plot(xg,yg,'k'); hold on;
plot(xt,yt,'r');
plot(xi,yi,'g'); hold off;
axis equal;

% plot grip aperture
subplot(2,1,2); hold off;
plot(tr,apt,'b');
[mga,maxt] = max(apt);
text(tr(maxt),mga+0.02,sprintf('%0.3f',mga));