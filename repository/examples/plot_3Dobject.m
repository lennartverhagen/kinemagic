% draw a prism slanted from 0 to 90 degrees in 15 degree steps
rot = [0:15:90];

% use patches to make cube
figure;
%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
xp = [0 1 1 0 0 1 1 0];
yp = [0 0 1 1 0 0 1 1];
zp = [0 0 0 0 1 1 1 1];
sc = [6 2 6];
dy = repmat(max(sc)*2,size(yp));
fmat = [1 2 6 5
        2 3 7 6
        3 4 8 7
        4 1 5 8
        1 2 3 4
        5 6 7 8];
fvcmat = [0.8 0.8 0.8
          0.8 0.8 0.8
          0.8 0.8 0.8
          0.8 0.8 0.8
          0.8 0.8 0.8
          0.8 0.8 0.8];

% plot and rotate prism
xdir = [1 0 0];
center = [3 1 3];
for r = 1:length(rot)
    yshift = dy*(r-1);
    vmat = [xp'*sc(1) yshift'+yp'*sc(2) zp'*sc(3)];
    hp = patch('vertices',vmat,'faces',fmat);
    %hp = patch('vertices',vmat,'faces',fmat,...
    %           'FaceVertexCdata',fvcmat,'FaceColor','flat');
    hold on;
    rotate(hp,xdir,-rot(r),center+[0 yshift(1) 0]);
end
axis equal; view(60,20);

% set axis
axis([0 6 -2 (yshift(1)+dy(1)/2) -2 8]);

% save the figure
set(gcf,'Color','w'); set(gca,'Box','off');
set(get(gca,'YLabel'),'FontSize',14);
set(get(gca,'XLabel'),'FontSize',14);
set(gca,'FontSize',12);
dir_root = fullfile(filesep,'home','action','lenver','Data','TMSEEG');
dir_report = fullfile('Group','kinematics','report');
figname = fullfile(dir_root,dir_report,'AIPpaper_object');
export_fig(figname,'-pdf');
