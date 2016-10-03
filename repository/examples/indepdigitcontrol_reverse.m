function indepdigitcontrol_reverse
%% Smeets-Brenner model of independent digit control in grasping
%--------------------------------------------------------------------------

%% set target and goal
%----------------------------------------
nf = 40;
n = 40;
% target front side
tf = [ -0.01*ones(nf,1) linspace(-0.03,0.03,nf)' ];
btf = pi*ones(nf,1);
% target top side
tt = [ linspace(-0.01,0.01,n)' 0.03*ones(n,1) ];
btt = pi/2*ones(n,1);
% target back side
tbb = [ 0.01*ones(nf,1) -linspace(-0.03,0.03,nf)' ];
btbb = zeros(nf,1);
% target bottom side
tb = [ -linspace(-0.01,0.01,n)' -0.03*ones(n,1) ];
btb = 3/2*pi*ones(n,1);
% construct target
itarget = [tf; tt; tbb; tb];
ibt = [btf; btt; btbb; btb];
ttarget = [tf; tt; tbb; tb];
tbt = [btf; btt; btbb; btb];
% goal positions
igoal = [0 0.03];
tgoal = [0 -0.03];

% remove the goal positions from the target otherwise the repulsive force
% from the goal will overcome the attraction at very small distances
rg = 0.002; % radius of digit approaching goal
gtdiff = bsxfun(@minus,igoal,itarget);
dg = euclid(gtdiff);
idx = dg<rg;
itarget(idx,:) = NaN;
ibt(idx,:) = NaN;
gtdiff = bsxfun(@minus,tgoal,ttarget);
dg = euclid(gtdiff);
idx = dg<rg;
ttarget(idx,:) = NaN;
tbt(idx,:) = NaN;

% rotate target and goal around reference
rot = -60/180*pi;
d = euclid(itarget);
beta = radangle(itarget);
itarget = xycomp(d,beta+rot);
ibt = ibt + rot;
d = euclid(ttarget);
beta = radangle(ttarget);
ttarget = xycomp(d,beta+rot);
tbt = tbt + rot;
d = euclid(igoal);
beta = radangle(igoal);
igoal = xycomp(d,beta+rot);
d = euclid(tgoal);
beta = radangle(tgoal);
tgoal = xycomp(d,beta+rot);

% target reference point in meters [x y]
ref = [0.25 0.30];
ntp = size(itarget,1);
itarget = itarget + repmat(ref,ntp,1);
ntp = size(ttarget,1);
ttarget = ttarget + repmat(ref,ntp,1);
igoal = igoal + ref;
tgoal = tgoal + ref;

% starting position
startpos = [0 0];

%% set obstable
%----------------------------------------
n = 200;
%obstacle = [ linspace(0.15,goal_i(1)+0.15,n)' linspace(0,goal_i(2),n)' ];
obstacle = [ linspace(0.20,ref(1)+0.10,n)' linspace(0.05,ref(2)-0.05,n)' ];
bo = radangle(obstacle(end,:) - obstacle(1,:)) + pi/2;

% starting position and velocity
ip = igoal + xycomp(rg,pi/2+rot);
iv = [0 0];
tp = tgoal + xycomp(rg,3/2*pi+rot);
tv = [0 0];

% set time steps
dt = 0.001;  % in 10 ms steps
t = 0:dt:1.4; % from 0 to 1 second 

% set scaling factors
As = 0.7;
Ag = 0.01;
Rt = 0.1;
Ro = 0.001;
K  = 0.05;
Ds = 0.3;
Dg = 0.1;

% initialize position, velocity and forces in memory
ipos    = nan(length(t),2); tpos    = nan(length(t),2);
ivel    = nan(length(t),2);	tvel    = nan(length(t),2);
iF      = nan(length(t),2); tF      = nan(length(t),2);
iFss    = nan(length(t),2); tFss    = nan(length(t),2);
iFsg    = nan(length(t),2); tFsg    = nan(length(t),2);
iFrt    = nan(length(t),2); tFrt    = nan(length(t),2);
iFro    = nan(length(t),2); tFro    = nan(length(t),2);
iFb     = nan(length(t),2); tFb     = nan(length(t),2);
iFds    = nan(length(t),2); tFds    = nan(length(t),2);
iFdg    = nan(length(t),2); tFdg    = nan(length(t),2);
ipos(1,:) = ip;              tpos(1,:) = tp;
ivel(1,:) = iv;              tvel(1,:) = tv;

%% loop over time steps 
%----------------------------------------
j = 1;
while j < length(t)
        
    % position and velocity
    ip = ipos(j,:);
    iv = ivel(j,:);
    tp = tpos(j,:);
    tv = tvel(j,:);
    
    %% index finger
    %----------------------------------------
    % get forces
    iFss(j,:)	= As .* Fss_fun(ip,startpos,rg*3);       % startpos attroactor
    iFsg(j,:)	= Ag .* Fsg_fun(ip,igoal,rg*3);          % goal repulsor
    iFrt(j,:)   = Rt .* Frt_fun(ip,iv,itarget,ibt);      % target repulsor
    iFro(j,:)   = Ro .* Fro_fun(ip,iv,obstacle,bo);      % obstacle repulsor
    iFb(j,:)	= K  .* Fb_fun(ip,tp);                   % digits belong to one hand
    iFds(j,:)	= Ds .* Fd_fun(ip,iv,startpos,tp,startpos,rg*3);	% damping and simultaneous arrival
    iFdg(j,:)	= Dg .* Fdg_fun(ip,iv,igoal,tp,tgoal,rg*3);       	% damping and simultaneous arrival
    
    % combine all forces acting on the digit
    iF(j,:) = iFss(j,:) + iFsg(j,:) + iFrt(j,:) + iFro(j,:) + iFb(j,:) + iFds(j,:) + iFdg(j,:);
    
    % calculate new position
    % the digits have an idealized mass, making the acceleration equal to
    % the force
    ivel(j+1,:) = iv + iF(j,:).*dt;
    ipos(j+1,:) = ip + ivel(j+1,:).*dt;
    
    %% thumb
    %----------------------------------------
    % get forces
    tFss(j,:)	= As .* Fss_fun(tp,startpos,rg*3);
    tFsg(j,:)	= Ag .* Fsg_fun(tp,tgoal,rg*3);
    tFrt(j,:)   = Rt .* Frt_fun(tp,tv,ttarget,tbt);
    tFro(j,:)   = Ro .* Fro_fun(tp,tv,obstacle,bo);
    tFb(j,:)	= K  .* Fb_fun(tp,ip);
    tFds(j,:)	= Ds .* Fd_fun(tp,tv,startpos,ip,startpos,rg*3);
    tFdg(j,:)	= Dg .* Fdg_fun(tp,tv,tgoal,ip,igoal,rg*3);
    
    % combine all forces acting on the digit
    tF(j,:) = tFss(j,:) + tFsg(j,:) + tFrt(j,:) + tFro(j,:) + tFb(j,:) + tFds(j,:) + tFdg(j,:);
    
    % calculate new position
    % the digits have an idealized mass, making the acceleration equal to
    % the force
    tvel(j+1,:) = tv + tF(j,:).*dt;
    tpos(j+1,:) = tp + tvel(j+1,:).*dt;
    
%     % make object solid
%     if pos(j+1,1) > min(target(:,1)) && pos(j+1,1) < max(target(:,1)) && ...
%        pos(j+1,2) > min(target(:,2)) && pos(j+1,2) < max(target(:,2))
%         if pos(j,1) <= pos(j+1,1)
%             pos(j+1,2) = min(target(:,2));
%         else
%             pos(j+1,1) = max(target(:,1));
%         end
%         vel(j+1,1) = 0;
%         vel(j+1,2) = 0;
%         break
%     end
    
    % stop when the startpos is almost reached to avoid x/0
    if euclid(startpos-ipos(j+1,:)) < rg || ...
       euclid(startpos-tpos(j+1,:)) < rg
        ivel(j+1,:) = [NaN NaN];
        tvel(j+1,:) = [NaN NaN];
        break
    end
    
%     % stop if the velocity is unstable
%     if j > 2 && ...
%        ((sign(ivel(j,1)) ~= sign(ivel(j-1,1)) && sign(ivel(j,1)) ~= sign(ivel(j+1,1))) || ...
%        (sign(ivel(j,2)) ~= sign(ivel(j-1,2)) && sign(ivel(j,2)) ~= sign(ivel(j+1,2))) || ...
%        (sign(tvel(j,1)) ~= sign(tvel(j-1,1)) && sign(tvel(j,1)) ~= sign(tvel(j+1,1))) || ...
%        (sign(tvel(j,2)) ~= sign(tvel(j-1,2)) && sign(tvel(j,2)) ~= sign(tvel(j+1,2))))
%         ivel(j-1,:) = [NaN NaN]; ivel(j,:) = [NaN NaN]; ivel(j+1,:) = [NaN NaN];
%         ipos(j-1,:) = [NaN NaN]; ipos(j,:) = [NaN NaN]; ipos(j+1,:) = [NaN NaN];
%         tvel(j-1,:) = [NaN NaN]; tvel(j,:) = [NaN NaN]; tvel(j+1,:) = [NaN NaN];
%         tpos(j-1,:) = [NaN NaN]; tpos(j,:) = [NaN NaN]; tpos(j+1,:) = [NaN NaN];
%         break
%     end
    
%     % stop if the velocity is extreme
%     if any(ivel(j+1,:)>1.5) || ...
%        any(tvel(j+1,:)>1.5)
%         ivel(j+1,:) = [NaN NaN];
%         tvel(j+1,:) = [NaN NaN];
%         break
%     end

    % loop increment
    j = j + 1;
    
end

% plot target, goal and obstacle
figure;
subplot(2,2,1); hold on;
plot(itarget(:,1),itarget(:,2),'k','linewidth',2);
plot(ttarget(:,1),ttarget(:,2),'k','linewidth',2);
plot(igoal(1),igoal(2),'ko','linewidth',2);
plot(tgoal(1),tgoal(2),'ko','linewidth',2);
plot(obstacle(:,1),obstacle(:,2),'r--','linewidth',2);

% plot position
subplot(2,2,1);
plot(ipos(:,1),ipos(:,2),'b');
plot(tpos(:,1),tpos(:,2),'g');
axis equal; hold off;
subplot(4,2,5); hold on;
plot(t,ipos(:,1),'b');
plot(t,tpos(:,1),'g');
subplot(4,2,7); hold on;
plot(t,ipos(:,2),'b');
plot(t,tpos(:,2),'g');

% plot grip aperture
subplot(4,2,2); hold on;
gripapt = euclid(ipos-tpos);
plot(t,gripapt,'k');
[mga,maxt] = max(gripapt);
text(t(maxt),mga+0.02,sprintf('%0.3f',mga));

% plot velocity
subplot(4,2,4); hold on;
plot(t,euclid(ivel),'b');
plot(t,euclid(tvel),'g');
subplot(4,2,6); hold on;
plot(t,ivel(:,1),'b');
plot(t,tvel(:,1),'g');
subplot(4,2,8); hold on;
plot(t,ivel(:,2),'b');
plot(t,tvel(:,2),'g');


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

function Fss = Fss_fun(p,startpos,rg)
%% move digits back to start positions
%----------------------------------------
% Fss = A./dg
%
% Fss: attracting force
% A: free scaling parameter
% dg: distance of digit to start
A = 1;
gsdiff = startpos-p;
ds = euclid(gsdiff);

% the attracting force should not go to infinity with close approaches
if ds < rg
    ds = rg;
end

Fss = A./ds;
% get x and y components of force
Fss = xycomp(Fss,radangle(gsdiff));

function Fsg = Fsg_fun(p,goal,rg)
%% move digits back from goal positions
%----------------------------------------
% Fsg = A./dg
%
% Fsg: repulsing force
% A: free scaling parameter
% dg: distance of digit to goal
A = 1;
gpdiff = goal-p;
dg = euclid(gpdiff);

% the force works opposite the velocity
beta = radangle(gpdiff)+pi;
beta(isnan(beta)|isinf(beta)) = 0;
% describe beta between 0 and 2pi
beta = beta + 2*(beta<0)*pi;

% the attracting force should not go to infinity with close approaches
if dg < rg
    dg = rg;
end

Fsg = A./dg;
% get x and y components of force
Fsg = xycomp(Fsg,beta);

function Frt = Frt_fun(p,v,target,targetbeta)
%% avoid hitting target
%----------------------------------------
% Frt = Rt * v * Z[(cos(bt-bv) + 1).^2)/dt.^2]ds
%
% Frt: repulsive force from target surface
% Rt: free scaling parameter
% v: velocity
% Z[...]ds: intregral over target surface points, implemented as the mean
% bt: angle from digit to target surface point
% bv: angle of velocity
% dt: distance to target surface point
Rt = 1;
tpdiff = target-repmat(p,size(target,1),1);
dt = euclid(tpdiff);
bt = radangle(tpdiff);
bv = radangle(v);
if isnan(bv)||isinf(bv), bv = 0;   end
ds = 1/length(dt);  % here the target surface is normalized to 1. In effect, a large target is as repellent as a small one.

% determine field of view based on the assumption that the target is
% defined clockwise. ignore target surface that is out of field of view and
% far away
%idx = gradient(bt)<=0 & dt>0.03;
%targetbeta(idx) = NaN;

flg = 'targetbeta';
switch flg
    case 'mean(bt)'
        Z = nansum(((cos(bt-bv)+1).^2)./dt.^2).*ds;
        Frt = Rt * euclid(v) * Z;
        % the force works away from the mean angle between the digit and surface
        beta = mean(bt)+pi;
        beta(isnan(beta)|isinf(beta)) = 0;
        % describe beta between 0 and 2pi
        beta = beta + 2*(beta<0)*pi;
        % get x and y components of force
        Frt = xycomp(Frt,beta);
    case 'targetbeta'
        % the direction of force can be different for every surface point
        Z = ((cos(bt-bv)+1).^2)./dt.^2.*ds;
        Frt = Rt * euclid(v) * Z;
        % get x and y components of force
        Frt = xycomp(Frt,targetbeta);
        Frt = nansum(Frt,1);
end

function Fro = Fro_fun(p,v,obstacle,obsbeta)
%% avoid hitting obstacle
%----------------------------------------
% Fro = Ro * v * Z[(cos(bo) + 1).^2)/do.^4]ds
%
% Fro: repulsive force from object surface
% Ro: free scaling parameter
% v: velocity
% Z[...]ds: intregral over target surface points
% bo: angle between direction of movement and object surface point
% do: distance to object surface point
Ro = 1;
opdiff = obstacle-repmat(p,size(obstacle,1),1);
do = euclid(opdiff);
bo = radangle(opdiff);
bv = radangle(v);
if isnan(bv)||isinf(bv), bv = 0;   end
ds = 1/length(do);  % here the target surface is normalized to 1. In effect, a large target is as repellent as a small one.

flg = 'obsbeta';
switch flg
    case 'mean(bo)'
        Z = sum(((cos(bo-bv)+1).^2)./do.^4).*ds;
        Fro = Ro * euclid(v) * Z;
        % the force works away from the mean angle between the digit and surface
        beta = mean(bo)+pi;
        beta(isnan(beta)|isinf(beta)) = 0;
        % describe beta between 0 and 2pi
        beta = beta + 2*(beta<0)*pi;
        % get x and y components of force
        Fro = xycomp(Fro,beta);
    case 'obsbeta'
        % the direction of force can be different for every surface point
        Z = ((cos(bo-bv)+1).^2)./do.^4.*ds;
        Fro = Ro * euclid(v) * Z;
        % get x and y components of force
        Fro = xycomp(Fro,obsbeta);
        Fro = sum(Fro,1);
end

function Fb = Fb_fun(p,po)
%% limit distance between digits
%----------------------------------------
% Fb = K * (E-db)
%
% Fb: grip aperture limiting force
% K: free scaling parameter
% E: free scaling parameter
% db: distance between digits
K = 1;
E = 0.02;
db = euclid(po-p);
Fb = K .* (E-db);
% the force works in the direction of the other digit
beta = radangle(po-p);
beta(isnan(beta)|isinf(beta)) = 0;
% describe beta between 0 and 2pi
beta = beta + 2*(beta<0)*pi;
% get x and y components of force
Fb = xycomp(Fb,beta);

function Fd = Fd_fun(p,v,goal,po,goalo,rg)
%% add damping and support simultaneous arrival
%----------------------------------------
% Fd = D * (dgo./dg).^2 * (v./dg).^2
%
% Fd: damping force
% D: free scaling parameter
% dgo: distance from other digit to goal
% dg: distance from digit to goal
% v: velocity
D = 1;
gpdiff = goal-p;
dg = euclid(gpdiff);

% the dampling force should not go to infinity with close approaches
if dg < rg
    dg = rg;
end

gopodiff = goalo-po;
dgo = euclid(gopodiff);
Fd = D .* (dgo./dg).^2 .* (euclid(v)./dg).^2;
% % the force works away from the goal
% beta = radangle(gpdiff)+pi;
% the force works opposite the velocity
beta = radangle(v)+pi;
beta(isnan(beta)|isinf(beta)) = 0;
% describe beta between 0 and 2pi
beta = beta + 2*(beta<0)*pi;
% get x and y components of force
Fd = xycomp(Fd,beta);

function Fdg = Fdg_fun(p,v,goal,po,goalo,rg)
%% add damping and support simultaneous return from goal
%----------------------------------------
% Fd = D * (dgo./dg).^2 * (v./dg).^2
%
% Fd: damping force
% D: free scaling parameter
% dgo: distance from other digit to goal
% dg: distance from digit to goal
% v: velocity
D = 1;
gpdiff = goal-p;
dg = euclid(gpdiff);

% the dampling force should not go to infinity with close approaches
if dg < rg
    dg = rg;
end

gopodiff = goalo-po;
dgo = euclid(gopodiff);
Fdg = D .* (dgo./dg).^2; % .* (euclid(v)./dg).^2;
% % the force works away from the goal
% beta = radangle(gpdiff)+pi;
% the force works in the direction of the velocity
beta = radangle(v);
beta(isnan(beta)|isinf(beta)) = 0;
% describe beta between 0 and 2pi
beta = beta + 2*(beta<0)*pi;
% get x and y components of force
Fdg = xycomp(Fdg,beta);

function d = euclid(d)
%% Euclidian distance
%----------------------------------------
d = squeeze(sqrt(sum(d.^2,2)));

function beta = radangle(d)
%% Angle in radians
%----------------------------------------
a = d(:,1);
b = d(:,2);
% calculate angle of force
beta = atan(b./a);
beta(isnan(beta)|isinf(beta)) = 0;
% flip angle if negative distance in x direction
beta = beta - (a<0)*pi;
% describe beta between 0 and 2pi
beta = beta + 2*(beta<0)*pi;

function d = xycomp(d,beta)
%% describe a vector in x and y components
%----------------------------------------
dx = d .* cos(beta);
dy = d .* sin(beta);
d = [dx dy];
