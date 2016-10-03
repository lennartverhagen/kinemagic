%% Basics
%----------------------------------------
dir_root	= fullfile(filesep,'home','action','lenver','Data','TMSEEG');
subjects	= {      '',      '',      '',      '',      '',...
                     '','subj07','subj08','subj09','subj10',...
               'subj11','subj12','subj13','subj14','subj15',...
               'subj16','subj17','subj18','subj19','subj20',...
               'subj21','subj22','subj23','subj24','subj25',...
               'subj26','subj27','subj28','subj29','subj30',...
               'subj31','subj32','subj33','subj34','subj35'};
sessions    = {'vertex','V6A','AIP'};

allsubj     = 8:35;
%badcapfit   = [8 10];
%badstereo   = [6 15 16];
%badvergence = [23 27];
%nogoodERD   = [10 13 22 23];
%badEEGart    = [4 9 18 31 32];
badnodata   = [8 35];
bad2fast    = [9 13 16];
badsubj     = unique([badnodata bad2fast]);
goodsubj    = setdiff(allsubj,badsubj);

subjsel	= goodsubj;
sesssel	= 1:3;

cfg = [];
cfg.dataset = '_MASKpchip_FILTlp15_ANA';
param = 'mga';
vision = [1 2];
slant = 0:15:90;
fname = 'slant';
marker = 'mean_thumbindex';

epe = nan(length(subjsel),length(sesssel),length(slant));
for s = 1:length(subjsel)
    cfg.subj = subjects{subjsel(s)};
    for ss = 1:length(sesssel)
        cfg.sess = sessions{sesssel(ss)};
        cfg.dir.proc = fullfile(dir_root,cfg.subj,'kinematics','ana');
        
        load(fullfile(cfg.dir.proc,sprintf('%s_%s%s_CFG',cfg.subj,cfg.sess,cfg.dataset)),'cfg');
        
        trl = cfg.trl;
        vars = cfg.vars;
        CI = cfg.eposerr.(fname).(marker);
        
        % end-point-error
        epe(s,ss,:) = CI.xyz.x * 1e2;
        %epe(s,ss,:) = CI.yz.area * 1e4;
        %epe(s,ss,:) = CI.xyz.volume * 1e6;
        
    end
end

%figure(s+1);
figure;
hold off
avg = squeeze(mean(epe));
err = squeeze(sem(epe));
errorbar(-1:15:89,avg(1,:),err(1,:),'ro--','LineWidth',1); hold on
errorbar(0:15:90,avg(2,:),err(2,:),'go--','LineWidth',1);
errorbar(1:15:91,avg(3,:),err(3,:),'bo--','LineWidth',1);
avg = squeeze(mean(mean(epe,2)));
err = squeeze(sem(mean(epe,2)));
errorbar(2:15:92,avg,err,'ko--','LineWidth',2);





