figure;
r = 1;
dat = data.pos{r};
hold off;
m = 1;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'r');
hold on;
m = 2;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'g');
m = 3;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'b');
m = 4;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'k');
axis equal

figure;
m = 1;
plot(squeeze(dat(m,1,:)),'r')

figure;
r = 1;
dat = data.pos{r};
hold off;
m = 2;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'r');
hold on;
m = 3;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'g');
m = 4;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'k');
m = 5;
plot3(squeeze(dat(m,1,:)),squeeze(dat(m,2,:)),squeeze(dat(m,3,:)),'b');
axis equal
