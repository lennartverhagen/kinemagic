function CP_plot_trajectory(cfg)
%--------------------------------------------------------------------------
% This function plots trajectory dynamics (see km_plot files)
% and performs cluster statistics (to do so
% stop the file before plotting).
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------

% set configuration
cfg = km_setcfg(cfg,'plot');


%% Get data for plotting
%----------------------------------------
% Load raw dataset(s)
[cfg,Sparam,SCI,Sepe,Stseries] = km_plot_getdata(cfg);

% Construct contrasts
%[cfg,Stseries] = km_plot_getcontr(cfg,Stseries);

% Get descriptives
[cfg,param,tseries] = km_plot_descriptives(cfg,Sparam,SCI,Sepe,Stseries);
clear Sparam Stseries;

% which factor?
%fact = 'actionxsel_P1';
%fact = 'actionxsel_P1';
fact = 'actionxaddresseexsel_P1';
% which sensor?
sns_all = {'i'};
nsns = length(sns_all);
% which type?
fldname = 'CIproject';
% how many points on the circle?
npnt = 51;
q = linspace(0,2*pi,npnt)';
pnts = zeros(length(q),3);
pnts(:,[2 3]) = [cos(q) sin(q)];
    
%initialize
all_area_all = cell(1,nsns);
all_area_all_mean = cell(1,nsns);
all_area_all_sem = cell(1,nsns);
all_pnts_all = cell(1,nsns);
all_pnts_all_mean = cell(1,nsns);
all_pos_all_mean = cell(1,nsns);
all_area_comminstr_mean = cell(1,nsns);
all_area_comminstr_sem = cell(1,nsns);
for isns = 1:nsns
sns = sns_all{isns};

% collect data
dims = size(cfg.traject{1}{1}.(fact).(sns).pos);
nax  = dims(1);
nlvl = dims(2);
nrpt = dims(3);
nsubj = length(cfg.traject);
area_all = nan(nlvl,nrpt,nsubj);
pnts_all = nan(npnt,nax,nlvl,nrpt,nsubj);
pos_all = nan(nax,nlvl,nrpt,nsubj);
gra_all = nan(nax,nlvl,nrpt,nsubj);
vec_all = nan(nax,nax,nlvl,nrpt,nsubj);
val_all = nan(nax,nax,nlvl,nrpt,nsubj);
k_all = nan(nlvl,nrpt,nsubj);
%subj_bad = [1:7 9:10];
subj_bad = [];

for s = 1:nsubj
    if ismember(s,subj_bad); continue; end
    area_all(:,:,s) = cfg.traject{s}{1}.(fact).(sns).(fldname).area;
    vec_all(:,:,:,:,s) = cfg.traject{s}{1}.(fact).(sns).(fldname).vec;
    val_all(:,:,:,:,s) = cfg.traject{s}{1}.(fact).(sns).(fldname).val;
    pnts_all(:,:,:,:,s) = cfg.traject{s}{1}.(fact).(sns).(fldname).pnts;
    k_all(:,:,s) = cfg.traject{s}{1}.(fact).(sns).(fldname).k;
    pos_all(:,:,:,s) =  cfg.traject{s}{1}.(fact).(sns).pos;
    gra_all(:,:,:,s) =  cfg.traject{s}{1}.(fact).(sns).gra;
end

% % restructure eigvec and eigval of each ellipsoid to match its neighbour,
% % starting with the middle point
% [vec_all, val_all] = km_traject_vectors_align(vec_all,val_all,gra_all);
% 
% % convolve a circle (pnts) with the principle axis of the ellipse
% pnts_all = km_traject_pnts_conv(vec_all,val_all,pos_all,k_all,pnts);
% 
% % align points by setting the origin of each ellipsoide to the point with highest x value
% pnts_all = km_traject_pnts_align(pnts_all,pos_all);
% 
% % calculate new area
% area_all = km_traject_areavolume(val_all,k_all);
% 
% % The average gradient can be considered a normal vector describing a plane
% % that intersects with the ellipsoid. This intersection is an ellipse
% % orthogonal to the direction of the movement.
% for s = 1:nsubj
%     if ismember(s,subj_bad); continue; end
%     for slnt = 1:7
%         % calculate direction of trajectory as gradient in xyz
%         gra = gradient(squeeze(avgpos_all(slnt,:,:,s))');
%         for i = 1:nrpt
%             % retrieve current eigenvalues
%             eigval_curr = diag(squeeze(eigval_all(slnt,i,:,:,s)));
%             % retrieve current eigenvectors
%             eigvec_curr = squeeze(eigvec_all(slnt,i,:,:,s));
%             % retrieve normal vector (gradient of average trajectory)
%             normvect_curr = gra(:,i);
%             % find the intersecting ellipse
%             [eigval_new, eigvec_new] = planeEllipsoidIntersection(eigval_curr,normvect_curr,eigvec_curr);
%             % match structure of eigenvalues (ascending order)
%             [eigval_new,ni] = sort(eigval_new);
%             eigvec_new = eigvec_new(:,ni);
%             % place back
%             eigval_all(slnt,i,:,:,s) = diag(eigval_new);
%             eigvec_all(slnt,i,:,:,s) = eigvec_new;
%         end
%     end
% end

% for s = 3
%     for slnt = [1 6]
%         i = 90:nrpt;
%         figure(10*s+slnt); hold on;
%         xyz = squeeze(pnts_all(slnt,i,:,:,s));
%         p1 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3)); %view(3);
%         p2 = plot3fix(squeeze(avgpos_all(slnt,i,:,s))','k','linewidth',2);
%         colormap gray; axis tight; axis equal; view(-70,10);
%     end
% end

% sort area over factors
switch fact
    case {'actionxsel_P1','actionxgoal_P1'}
        area_comm = mean(area_all([1 3 5],:,:),1);
        area_instr = mean(area_all([2 4 6],:,:),1);
        area_comminstr = cat(1,area_comm,area_instr);
    case {'actionxaddresseexsel_P1','actionxaddresseexgoal_P1','actionxtime_halfxsel_P1'}
        area_commleft = mean(area_all([1:4:12],:,:),1);
        area_commright = mean(area_all([3:4:12],:,:),1);
        area_instrleft = mean(area_all([2:4:12],:,:),1);
        area_instrright = mean(area_all([4:4:12],:,:),1);
        area_comminstr = cat(1,area_commleft,area_commright,area_instrleft,area_instrright);
end

%store
all_area_all{isns} = area_all;
all_pnts_all{isns} = pnts_all;
if nsubj == 1
    all_area_all_mean{isns} = area_all;
    all_pnts_all_mean{isns} = pnts_all;
    all_pos_all_mean{isns} = pos_all;
    all_area_comminstr_mean{isns} = area_comminstr;
else
    all_area_all_mean{isns} = nanmean(area_all,3);
    all_area_all_sem{isns} = nansem(area_all,0,3);
    all_pnts_all_mean{isns} = nanmean(pnts_all,5);
    all_pos_all_mean{isns} = nanmean(pos_all,4);
    all_area_comminstr_mean{isns} = nanmean(area_comminstr,3);
    all_area_comminstr_sem{isns} = nansem(area_comminstr,0,3);
end

end

switch fact
    case {'actionxsel_P1','actionxgoal_P1'}
        % plot average trajectory
        h_f1 = figure; hold on;
        plot3fix(permute(all_pos_all_mean{1}(:,[1 3 5],:),[1 3 2]),'r');
        plot3fix(permute(all_pos_all_mean{1}(:,[2 4 6],:),[1 3 2]),'b');
        axis([-0.05 0.3 -0.7 -0.4 0.1 0.5]);
        daspect([1,1,1]); view(-70,10);
        
        % plot trajectory variability
        h_f1 = figure; hold on;
        xyz = permute(squeeze(all_pnts_all_mean{1}(:,:,1,:)),[1 3 2]);
        p11 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','red','EdgeColor','none');
        xyz = permute(squeeze(all_pnts_all_mean{1}(:,:,3,:)),[1 3 2]);
        p12 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','red','EdgeColor','none');
        xyz = permute(squeeze(all_pnts_all_mean{1}(:,:,5,:)),[1 3 2]);
        p13 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','red','EdgeColor','none');
        xyz = permute(squeeze(all_pnts_all_mean{1}(:,:,2,:)),[1 3 2]);
        p21 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','blue','EdgeColor','none');
        xyz = permute(squeeze(all_pnts_all_mean{1}(:,:,4,:)),[1 3 2]);
        p22 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','blue','EdgeColor','none');
        xyz = permute(squeeze(all_pnts_all_mean{1}(:,:,6,:)),[1 3 2]);
        p23 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','blue','EdgeColor','none');
        axis([-0.05 0.3 -0.7 -0.4 0.1 0.5]);
        daspect([1,1,1]); view(-70,10);
        h_cl1 = camlight(20,40); lighting phong;
        
        
        % plot area
        figure; hold on;
        x = linspace(0,100,nrpt);
        y = all_area_comminstr_mean{1}(1,:);
        plot(x,y,'linewidth',2,'color',str2rgb('red'));
        if nsubj > 1
            yerr = all_area_comminstr_sem{1}(1,:);
            ploterrorcloud(x,y,yerr,str2rgb('red'));
        end
        y = all_area_comminstr_mean{1}(2,:);
        plot(x,y,'linewidth',2,'color',str2rgb('blue'));
        if nsubj > 1
            yerr = all_area_comminstr_sem{1}(2,:);
            ploterrorcloud(x,y,yerr,str2rgb('blue'));
        end
        axis([0 100 0 100]);
        %figname = fullfile(dir_report,'area_MPJ');
        %export_fig(figname,'-pdf');
        
        
        
        
    case {'actionxaddresseexsel_P1','actionxaddresseexgoal_P1','actionxtime_halfxsel_P1'}
        % color setttings
        cl = str2rgb('robk');
        
        % plot average trajectory
        h_f1 = figure; hold on;
        h(1,:) = plot3fix(permute(all_pos_all_mean{1}(:,[1:4:12],:),[1 3 2]),'color',cl(1,:),'linewidth',2);
        h(2,:) = plot3fix(permute(all_pos_all_mean{1}(:,[3:4:12],:),[1 3 2]),'color',cl(2,:),'linewidth',2);
        h(3,:) = plot3fix(permute(all_pos_all_mean{1}(:,[2:4:12],:),[1 3 2]),'color',cl(3,:),'linewidth',2);
        h(4,:) = plot3fix(permute(all_pos_all_mean{1}(:,[4:4:12],:),[1 3 2]),'color',cl(4,:),'linewidth',2);
        axis([-0.05 0.3 -0.7 -0.4 0.1 0.5]);
        legend(h(:,1),'comm-left','comm-right','instr-left','instr-right','location','NorthOutside');
        xlabel('x'); ylabel('y'); zlabel('z');
        daspect([1,1,1]); view(-70,10);
        
        % plot trajectory variability
        h_f1 = figure; hold on;
        for lvl = 1:4
            for target = lvl:4:12
                xyz = permute(squeeze(all_pnts_all_mean{1}(:,:,target,:)),[1 3 2]);
                h(lvl) = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor',cl(lvl,:),'EdgeColor','none');
            end
        end
        axis([-0.05 0.3 -0.7 -0.4 0.1 0.5]);
        legend(h(:),'comm-left','comm-right','instr-left','instr-right','location','NorthOutside');
        daspect([1,1,1]); view(-70,10);
        h_cl1 = camlight(20,40); lighting phong;
        
        % plot area
        figure; hold on;
        x = linspace(0,100,nrpt);
        for lvl = 1:4
            y = all_area_comminstr_mean{1}(lvl,:);
            h(lvl) = plot(x,y,'linewidth',2,'color',cl(lvl,:));
            if nsubj > 1
                yerr = all_area_comminstr_sem{1}(lvl,:);
                ploterrorcloud(x,y,yerr,cl(lvl,:));
            end
        end
        axis([0 100 0 100]);
        legend(h(:),'comm-left','comm-right','instr-left','instr-right','location','NorthOutside');
        %figname = fullfile(dir_report,'area_MPJ');
        %export_fig(figname,'-pdf');
        
end

return

% save
dir_report = fullfile(filesep,'home','action','ankmur','CommPointing','tempplots_kin');
if ~exist(dir_report,'dir'), mkdir(dir_report); end
figname = fullfile(dir_report,'EPE');
export_fig(figname,'-pdf');
%export_fig(figname,'-png','-r600');


% hack plot
for lvl = 1:7
xyz = squeeze(all_pnts_all{3}(lvl,:,1:2:25,:,:));
for s = 1:23
    figure(s);
    dat = squeeze(xyz(1:10:end,:,:,s));
    surf(dat(:,:,1),dat(:,:,2),dat(:,:,3));
    axis equal;
    axis([0.1 0.3 -0.5 -0.1 0.05 0.55]);
    view(-70,10); %colormap gray;
end
end


% stats
p = 1:nrpt;
for i = 1:nrpt
    %a = sqrt(squeeze(area_verthori(1,i,:)));
    %b = sqrt(squeeze(area_verthori(2,i,:)));
    a = sqrt(squeeze(area_comm(1,i,:)));
    b = sqrt(squeeze(area_instr(1,i,:)));
    %a = squeeze(pos_all(1,10,i,:)); % (x,leftcommleft,allpoints,subj)
    %b = squeeze(pos_all(1,12,i,:)); % (x,leftcommright,allpoints,subj)
    [~,p(i)] = ttest(a,b);
end

%tmpHACK effect size calc
%Collapse data for action variability effect --> (whole cluster)
%cluster)
i = [129:183];
a = [squeeze(mean(area_comm(1,i,:),2)) squeeze(mean(area_instr(1,i,:),2))];
% save
header = 'variability_effectsize.txt';
save(header,'a','-ascii');
%Collapse data for action variability effect --> 153/300 (smallest p first
%cluster)
i = [201:300];
a = [squeeze(mean(area_comm(1,i,:),2)) squeeze(mean(area_instr(1,i,:),2))];
% save
header = 'variability2_effectsize.txt';
save(header,'a','-ascii');
%END tmpHACK

% cluster statistics
subjsel = setdiff(1:nsubj,subj_bad);
rptsel = 1:nrpt;
nrptsel = length(rptsel);
nsubjsel = length(subjsel);
dataA.dimord = 'subj_chan_time';
dataA.label = {'verthori'};
dataA.time = 1:nrptsel;
dataB = dataA;
dataA.avg = sqrt(reshape(squeeze(area_verthori(1,rptsel,subjsel))',[nsubjsel 1 nrptsel]));
dataB.avg = sqrt(reshape(squeeze(area_verthori(2,rptsel,subjsel))',[nsubjsel 1 nrptsel]));
% --> choose 2 conditions to contrast
% Left target (1 = commleft, 2 = instrleft, 3 = commright, 4 = instrright)
% Middle target (5 = commleft, 6 = instrleft, 7 = commright, 8 = instrright)
% Right target (9 = commleft, 10 = instrleft, 11 = commright, 12 = instrright)
dataA.avg = sqrt(reshape(squeeze(pos_all(1,1,rptsel,subjsel))',[nsubjsel 1 nrptsel]));
dataB.avg = sqrt(reshape(squeeze(pos_all(1,2,rptsel,subjsel))',[nsubjsel 1 nrptsel]));
%dataA.avg = reshape(squeeze(pos_all(1,10,rptsel,subjsel))',[nsubjsel 1 nrptsel]);
%dataB.avg = reshape(squeeze(pos_all(1,12,rptsel,subjsel))',[nsubjsel 1 nrptsel]);
dataA.avg = squeeze(pos_all(1,1,rptsel,subjsel))';
dataB.avg = squeeze(pos_all(1,2,rptsel,subjsel))';

% group design
ftcfg = [];
ftcfg.neighbours(1).label = '1';
ftcfg.neighbours(1).neighblabel = {};
ftcfg.neighbours(2).label = '2';
ftcfg.neighbours(2).neighblabel = {};
ftcfg.neighbours(3).label = '3';
ftcfg.neighbours(3).neighblabel = {};
ftcfg.method = 'montecarlo';
ftcfg.numrandomization = 10000;
ftcfg.correctm = 'cluster';
ftcfg.clusterstatistic = 'maxsum';
ftcfg.clusterthreshold = 'parametric';
ftcfg.clusteralpha = 0.05;
ftcfg.alpha = 0.05;
ftcfg.parameter = 'avg';
ftcfg.statistic	= 'depsamplesT';
ftcfg.design(1,1:2*nsubjsel) = [ones(1,nsubjsel) 2*ones(1,nsubjsel)];
ftcfg.design(2,1:2*nsubjsel) = [1:nsubjsel 1:nsubjsel];
ftcfg.ivar = 1;
ftcfg.uvar = 2;

% set path
addpath /home/common/matlab/fieldtrip
ft_defaults

% plotting
stat = ft_timelockstatistics(ftcfg,dataA.avg,dataB.avg);
%stat = ft_freqstatistics(ftcfg,dataA,dataB);
figure; 
plot(1:nrpt,p); 
hold on; 
plot(rptsel,stat.mask,'k');
plot(rptsel, stat.posclusterslabelmat,'r');
ylim([-0.1,1.1]);

% % first test for size
% ftcfg.clusterstatistic = 'maxsize';
% stat = ft_timelockstatistics(ftcfg,dataA,dataB);
% probsize = stat.prob;
% % then test for peak
% ftcfg.clusterstatistic = 'max';
% stat = ft_timelockstatistics(ftcfg,dataA,dataB);
% probpeak = stat.prob;
% % combine
% prob = min(probsize,probpeak)*2;
% prob(prob>1) = 1;
% 
% % plotting
% figure; plot(1:nrpt,area_verthori_p); hold on; plot(rptsel,prob,'b');
% 
% % compare that to the default
% ftcfg.clusterstatistic = 'maxsum';
% stat = ft_timelockstatistics(ftcfg,dataA,dataB);
% hold on; plot(rptsel,stat.prob,'k');


% calculate new points for verthori from eig variates
q = linspace(0,2*pi,length(q))';
rpt = zeros(length(q),3);
rpt(:,[2 3]) = [cos(q) sin(q)];
new_pnts_verthori_mean = nan(2,nrpt,length(q),3);
for i = 1:2
    for j = 1:nrpt
        tmp = rpt*sqrt(squeeze(eigval_verthori_mean(i,j,:,:)))*squeeze(eigvec_verthori_mean(i,j,:,:))';
        %new_pnts_verthori_mean(i,j,:,:) = bsxfun(@plus,[0 j/100000 j/100000],tmp);
        tmp_offset = squeeze(avgpos_verthori_mean(i,j,:));
        new_pnts_verthori_mean(i,j,:,:) = bsxfun(@plus,tmp_offset',tmp);
    end
end

% average points over subjects


% plot surface of all subjects
n = 1:nrpt;
for s = 21
    if ismember(s,subj_bad); continue; end
    for lvl = 1:7;
        figure(10*s+lvl); hold on;
        xyz = squeeze(pnts_all(lvl,n,:,:,s));
        p1 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3));
        axis equal; axis tight; colormap gray; view(-70,20); %camlight;
    end
end


figure; hold on;
for i = 1:9:nrpt
    plot3fix(pnts_verthori_mean(1,i,:,:),'linewidth',2,'color','r');
end
for i = 1:9:nrpt
    plot3fix(pnts_verthori_mean(2,i,:,:),'linewidth',2,'color','g');
end
legend('vert','hori');

%

figure; hold on;
xyz = squeeze(pnts_verthori_mean(1,1:2:96,:,:));
%p = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3));
p1 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','red','EdgeColor','none');
%xyz = squeeze(pnts_verthori_mean(2,1:2:96,:,:));
%p2 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','green','EdgeColor','none');
h_cl = camlight(40,10); lighting phong; daspect([1,1,1]); view(3);
%colormap gray; camlight; axis tight; axis equal; %alpha(0.3);
h_cl = camlight(h_cl,40,10);



figure; hold on;
xyz = squeeze(new_pnts_verthori_mean(1,1:2:96,:,:));
%p = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3));
p1 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','red','EdgeColor','none');
%xyz = squeeze(pnts_verthori_mean(2,1:2:96,:,:));
%p2 = surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'FaceColor','green','EdgeColor','none');
h_cl = camlight(40,10); lighting phong; daspect([1,1,1]); view(3);




% p2 = patch(isosurface(xyz,.5),...
%     'FaceColor','red','EdgeColor','none');
% isonormals(data,p2)
% view(3); daspect([1 1 1]); axis tight
% camlight;  camlight(-80,-10); lighting phong; 


function [tmp,ni] = matchvectors(tmp,tmp0,gra,nvect)
if nargin < 3
    nvect = 3;
end

% calculate angle between reference and current eigenvectors
theta = [vectangle(tmp0(:,1),tmp)
         vectangle(tmp0(:,2),tmp)
         vectangle(tmp0(:,3),tmp)];
         
if nvect == 3
    % find the smallest angle (largest cosine)
    T = abs(theta);
    ni = nan(1,3);
    while any(T(:))
        % find matches for all three vectors
        [m,n] = max(T,[],2);
        % appoint closest match first
        [~,nm] = max(m);
        ni(nm) = n(nm);
        % then continue to search
        T(nm,:) = 0; T(:,n(nm)) = 0;
    end
    % check for overlap
    if length(unique(ni)) < 3 || any(isnan(ni));
        error('this does not work')
    end
    new_sign = sign([theta(1,ni(1)) theta(2,ni(2)) theta(3,ni(3))]);
elseif nvect == 2
    % calculate angle between reference and current eigenvectors
    theta2 = [0 theta(2,2:3)];
    theta3 = [0 theta(2,2:3)];
    % find the smallest angle (largest cosine)
    [m2,n2] = max(abs(theta2));
    [m3,n3] = max(abs(theta3));
    % do not permit different vectors to be mapped to the same
    if m2 > m3
        n3 = setdiff(2:3,n2);
    else
        n2 = setdiff(2:3,n3);
    end
    ni = [1 n2 n3];
    new_sign = sign([1 theta2(n2) theta3(n3)]);
end

% map referenced vectors
tmp = bsxfun(@times,new_sign,tmp(:,ni));

% check if a vector matches the direction of movement very closely
gra_theta = vectangle(gra,tmp);
[m,n] = max(abs(gra_theta));
% if this is not the first vector, place it up front
if n ~= 1 && ( m>0.8 || all(m > abs(theta(:,n))) )
    tmpi = 1:3; tmpi(1) = n; tmpi(n) = 1;
    ni = ni(tmpi);
    tmp = tmp(:,tmpi);
    % check if the matched components match in direction and flip sign if
    % necessary
    tmp = bsxfun(@times,sign(vectangle(tmp0,tmp)),tmp);
end