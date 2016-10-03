% startup('km');

% defaults
scrsz = get(0,'ScreenSize');
dir_root = fullfile(filesep,'home','action','lenver','Data','TMSEEG');
dir_report = fullfile('Group','kinematics','report');

% select subject, session and trial
dataset = 'MASKpchip_FILTlp15_MOV_ANA';
subj = 'subj27';
sess = 'vertex';
trlsel  = 28;
dir_proc = fullfile(dir_root,subj,'kinematics','ana');

% load configuration and parameters
dcfg = load(fullfile(dir_proc,sprintf('%s_%s_%s_CFG',subj,sess,dataset)));
dcfg = dcfg.cfg;
% load data
tdata = load(fullfile(dir_proc,sprintf('%s_%s_%s_DATA',subj,sess,dataset)));
dat = struct('time',{[]},'dimord',{[]},'marker',{[]},'axis',{[]});

% select movement
mov = dcfg.movement;
% make the movemements a little bit longer...
mov = cellfun(@(x) x + [-0.5 1],mov,'UniformOutput',false);
% make sure that all trials have the same number of movements
mov = km_setequalnmov(mov);

% select time-series
tcfg                = [];
tcfg.pos.marker    	= {'MCP','t','i'};
tcfg.pos.axis    	= {'x','y','z'};
tcfg.vel.marker   	= {'MCP','t','i'};
tcfg.vel.axis     	= {'yz'};
tcfg.gripapt.marker	= 'gripti';
tcfg.gripapt.axis	= 'yz';
tcfg.gripaptvel.marker	= 'gripti';
tcfg.gripaptvel.axis 	= 'yz';
tcfg.gripori.marker	= 'gripti';
tcfg.gripori.axis	= 'yz';

% get time
movnr = km_getmovidx(tcfg,'vel');
tim = cellfun(@(x) x(movnr,:),mov,'UniformOutput',false);
dat.time = cell(1,length(tim));
for t = 1:length(tim)
    it = km_time2logic(tim{t},tdata.data.time{t});
    dat.time{t} = tdata.data.time{t}(it);
end
    
% loop over time series
tseries = {'pos','vel','gripapt','gripaptvel','gripori'};
for ts = 1:length(tseries)
    tsname = tseries{ts};
    % get time series
    dat.(tsname) = km_getmovdat(tcfg,tdata.data,mov,tsname);
    % store dimension order
    nd = ndims(dat.(tsname){1});
    if nd == 2
        dat.dimord.(tsname) = 'axis_time';
    elseif nd == 3
        dat.dimord.(tsname) = 'sensor_axis_time';
    else
        dat.dimord.(tsname) = 'unknown';
    end
    % store axis names
    dat.axis.(tsname) = tcfg.(tsname).axis;
    dat.marker.(tsname) = tcfg.(tsname).marker;
end

% get time
t = dat.time{trlsel};
t = round(t*1e9)/1e9;
nt = length(t);
mi = t >= 0.668 & t <= 1.584;
ti = t >= 0.668 & t <= 1.336;
ai = t >= 1.336 & t <= 1.584;
si = t >= 1.584 & t <= 1.584+0.6;

% get time-points
%mip = km_nearest(t,max(t(ai))-(max(t(ai))-min(t(ai)))*(4:-1:0));
mip = km_nearest(t,linspace(min(t(ti)),max(t(ai)),5));
tip = km_nearest(t,linspace(min(t(ti)),max(t(ti)),5));
aip = [find(ai,1,'first') find(ai,1,'last')];
sip = [find(si,1,'first') find(si,1,'last')];

% get parameters
x = reshape(dat.pos{trlsel}(:,1,:),[3 nt]);
y = reshape(dat.pos{trlsel}(:,2,:),[3 nt]);
z = reshape(dat.pos{trlsel}(:,3,:),[3 nt]);
v = reshape(dat.vel{trlsel},[3 nt]);
ga = dat.gripapt{trlsel};
go = dat.gripori{trlsel};
gv = dat.gripaptvel{trlsel};

% set time relative to movement onset
t = t - t(find(ti,1,'first'));

% reset position origin
x = x-min(x(:));
y = y-min(y(:));
z = z-min(z(:));

% plot time-points in y-z (flipped on y)
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
plot(-y(1:3,mip(1:3))',z(1:3,mip(1:3))','linestyle','none','marker','o','markersize',6,'markerfacecolor','b','markeredgecolor','none');
hold on;
plot(-y(1:3,mip(4:5))',z(1:3,mip(4:5))','linestyle','none','marker','o','markersize',6,'markerfacecolor','r','markeredgecolor','none');

% plot whole time-series in y-z (flipped on y)
plot(-y(1,ti)',z(1,ti)','b-.');
plot(-y(2,ti)',z(2,ti)','b:');
plot(-y(3,ti)',z(3,ti)','b--');
plot(-y(1:3,ai)',z(1:3,ai)','r');

% set axis
axis equal;
ax = axis;
axis([-0.3 0 0 0.45]);

% wrap-up and save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
figname = fullfile(dir_root,dir_report,'AIPpaper_pos_whole');
export_fig(figname,'-pdf');

% restrict data to zoom windows
xtmp = x; ytmp = y; ztmp = z;
xlim = [0  0.05];
ylim = [0.15 0.27];
zlim = [0.26 0.42];
x(x<min(xlim)) = NaN; x(x>max(xlim)) = NaN;
y(y<min(ylim)) = NaN; y(y>max(ylim)) = NaN;
z(z<min(zlim)) = NaN; z(z>max(zlim)) = NaN;

% zoom in on approach phase
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
plot3(xtmp(1:3,aip(2))',ytmp(1:3,aip(2))',ztmp(1:3,aip(2))','linestyle','none','marker','o','markersize',6,'markerfacecolor','r','markeredgecolor','none');
hold on; axis equal;
%plot3(x(1,ti)',y(1,ti)',z(1,ti)','b-.');
%plot3(x(2,ti)',y(2,ti)',z(2,ti)','b:');
%plot3(x(3,ti)',y(3,ti)',z(3,ti)','b--');
plot3(x(1:3,ti)',y(1:3,ti)',z(1:3,ti)','b-');
plot3(x(1:3,ai)',y(1:3,ai)',z(1:3,ai)','r');
plot3(x(1:3,si)',y(1:3,si)',z(1:3,si)','color',[102 102 102]/255);

% set axis
axis equal; view(-60,20);
% ax = axis;
axis([xlim ylim zlim]);

% save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
figname = fullfile(dir_root,dir_report,'AIPpaper_pos_zoom');
export_fig(figname,'-pdf');


% draw the object

% use patches to make cube
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
xp = [0 1 1 0 0 1 1 0];
yp = [0 0 1 1 0 0 1 1];
zp = [0 0 0 0 1 1 1 1];
vmat = [xp'*6 yp'*2 zp'*6];
fmat = [1 2 6 5
        2 3 7 6
        3 4 8 7
        4 1 5 8
        1 2 3 4
        5 6 7 8];
fvcmat = [0.2 0.2 0.2
          0.2 0.2 0.2
          0.2 0.2 0.5
          0.2 0.2 0.2
          0.2 0.2 0.2
          0.2 0.2 0.2];
%hp = patch('vertices',vmat,'faces',fmat,...
%           'FaceVertexCdata',fvcmat,'FaceColor','flat');
hp = patch('vertices',vmat,'faces',fmat);

% rotate 45 deg
xdir = [1 0 0];
center = [2 1 3];
rotate(hp,xdir,-45,center);
axis equal; view(-60,20);

% set axis
axis([0 6 -2 4 0 6]);

% save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
figname = fullfile(dir_root,dir_report,'AIPpaper_object');
export_fig(figname,'-pdf');


% plot velocity
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
plot(t(ti)',v(1,ti)','b-.');
hold on;
plot(t(ti)',v(2,ti)','b:');
plot(t(ti)',v(3,ti)','b--');
plot(t(ai)',v(1,ai)','r-.');
plot(t(ai)',v(2,ai)','r:');
plot(t(ai)',v(3,ai)','r--');
plot(t(si)',v(1,si)','color',[102 102 102]/255,'linestyle','-.');
plot(t(si)',v(2,si)','color',[102 102 102]/255,'linestyle',':');
plot(t(si)',v(3,si)','color',[102 102 102]/255,'linestyle','--');

% set axis
axis([0 1 0 1.5]);

% save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
figname = fullfile(dir_root,dir_report,'AIPpaper_vel');
export_fig(figname,'-pdf');


% plot grip aperture
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
plot(t(ti)',ga(ti)','b','linewidth',1);
hold on;
plot(t(ai)',ga(ai)','color',[255 153 51]/255);
plot(t(si)',ga(si)','color',[102 102 102]/255);

% set axis
axis([0 1 0.02 0.14]);

% save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
figname = fullfile(dir_root,dir_report,'AIPpaper_gripapt');
export_fig(figname,'-pdf');


% plot grip aperture orientation
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
plot(t(ti)',go(ti)','b');
hold on;
plot(t(ai)',go(ai)','color',[255 153 51]/255);
plot(t(si)',go(si)','color',[102 102 102]/255);

% set axis
axis([0 1 0 180]);

% save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
figname = fullfile(dir_root,dir_report,'AIPpaper_gripaptori');
export_fig(figname,'-pdf');


% plot grip aperture velocity
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
line([0 1],[0 0],'color',[102 102 102]/255,'linestyle','--');
hold on;
plot(t(ti)',gv(ti)','b','linewidth',1);
plot(t(ai)',gv(ai)','color',[255 153 51]/255);
plot(t(si)',gv(si)','color',[102 102 102]/255);

% set axis
axis([0 1 -0.8 0.4]);

% save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
figname = fullfile(dir_root,dir_report,'AIPpaper_gripaptvel');
export_fig(figname,'-pdf');