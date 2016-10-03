% simulate data
x = 0:0.01:25;
%y = abs(sin(x));
y = sin(x);

% keep consequitive movements separated
conseqflg = 'both';   % settings: {'no','onset','offset','both'}

% choose absolute data or not...
absflg = false;
if absflg, y = abs(y); end

plot(x,y);
clear data;
data(1,1,:) = y;

% simulate initial movement definitions
sel = (x>7 & x<8.5 & abs(y)>0.9) | (x>10.5 & x<12 & abs(y)>0.9) | (x>16.5 & x<18 & abs(y)>0.9);
mov = x(logic2idx(sel));

% plot initial selection
v = axis;
hold on;
for i = 1:size(mov,1)
    fill(mov(i,[1 2 2 1]),v([3 3 4 4]),'g','LineStyle','none','FaceAlpha',0.3);
end
hold off;
clear select;
select(1,1,:) = sel;

% and slide...
% examples: crit = {-Inf -0.6 0 0.6};
crit = -Inf;
select = km_slidesel(data,select,crit,conseqflg);

% plot new selection
mov = x(logic2idx(squeeze(select)));
hold on;
for i = 1:size(mov,1)
    fill(mov(i,[1 2 2 1]),v([3 3 4 4]),'r','LineStyle','none','FaceAlpha',0.2);
end
hold off;


% test what happens if selection is only top point of wave...

%%%% SECOND DATASET

% time points
x = -10:0.1:10;

% cosine movement
y = -cos(x)';

% normal distribution
ydist = pdf('norm',x,0,2)';

% create data by combining the normal distribution with the sine wave
dat = y .* ydist;
absdat = abs(dat);

% select movements
crit = 0.1;
sel = absdat > crit;

% test
slidecrit = 0.001;
newsel = km_slidesel(dat,sel,slidecrit,'merge');

% plot data
figure; hold on;
plot(x,absdat,'r');
plot(x,dat,'g');
v = axis;

% plot selection windows
mov = x(logic2idx(squeeze(sel)));
hold on;
for i = 1:size(mov,1)
    fill(mov(i,[1 2 2 1]),v([3 3 4 4]),'r','LineStyle','none','FaceAlpha',0.2);
end
hold off;

% plot selection windows
mov = x(logic2idx(squeeze(newsel)));
hold on;
for i = 1:size(mov,1)
    fill(mov(i,[1 2 2 1]),v([3 3 4 4]),'g','LineStyle','none','FaceAlpha',0.2);
end
hold off;