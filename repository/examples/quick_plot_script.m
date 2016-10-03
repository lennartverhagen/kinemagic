%param = {'rt','mt','td','mv','pv','tpv','rtpv','mga','tmga','rtmga','ega','ego'};

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

subjsel	= 34;
sesssel	= 1;

cfg = [];
cfg.dataset = '_MASKpchip_FILTlp15_ANA';
param = 'mt';
vision = [1 2];
slant = 0:15:90;

plot_effect = {'vision x slant'};
%plot_effect = {'slant'};

dat = nan(length(subjsel),length(sesssel),length(vision),length(slant));
for s = 1:length(subjsel)
    cfg.subj = subjects{subjsel(s)};
    for ss = 1:length(sesssel)
        cfg.sess = sessions{sesssel(ss)};
        cfg.dir.proc = fullfile(dir_root,cfg.subj,'kinematics','ana');
        
        load(fullfile(cfg.dir.proc,sprintf('%s_%s%s_CFG',cfg.subj,cfg.sess,cfg.dataset)),'cfg');
        
        trl = cfg.trl;
        vars = cfg.vars;
        
        % indeces
        for i = 1:length(vision)
            iv{i} = cfg.trl(:,8)==vision(i);
        end
        for i = 1:length(slant)
            io{i} = cfg.trl(:,9)==slant(i);
        end
        ip = ismember(vars,param);
        
        for i = 1:length(iv)
            for j = 1:length(io)
                if strcmpi(param,'ego')
                    dat(s,ss,i,j) = nanmean(cfg.trl(iv{i}&io{j},ip)) - slant(j);
                    %dat(s,ss,i,j) = nanmean(cfg.trl(iv{i}&io{j},ip));% - slant(j);
                else
                    dat(s,ss,i,j) = nanmean(cfg.trl(iv{i}&io{j},ip));
                end
            end
        end
        
    end
end


% m_bino = nan(size(bino,1),7);
% m_mono = nan(size(mono,1),7);
% for s = 1:size(bino,1)
%     m_bino(s,:) = nanmean([bino{s,1};bino{s,2};bino{s,3}]);
%     m_mono(s,:) = nanmean([mono{s,1};mono{s,2};mono{s,3}]);
%     if strcmpi(param,'ego')
%         m_bino(s,:) = m_bino(s,:) - (0:15:90);
%         m_mono(s,:) = m_mono(s,:) - (0:15:90);
%     end
%     %figure(s);
%     %hold off
%     %plot(0:15:90,m_bino(s,:),'ro--')
%     %hold on;
%     %plot(0:15:90, m_mono(s,:),'go--')
% end

nsubj = length(subjsel);
nsess = length(sesssel);
nvision = 2;
nslant = 7;

bino = reshape(dat(:,:,1,:),[nsubj nsess nslant]);
mono = reshape(dat(:,:,2,:),[nsubj nsess nslant]);
b_avg = reshape(mean(bino,1),[nsess nslant]);
m_avg = reshape(mean(mono,1),[nsess nslant]);
b_err = reshape(sem(bino-mono,0,1),[nsess nslant]);    % reshape(sem(bino),[nsubj nsess nslant])
m_err = reshape(sem(bino-mono,0,1),[nsess nslant]);    % reshape(sem(mono),[nsubj nsess nslant])
if nsess > 1
    % session x vision x slant
    figure(s+1);
    hold off
    errorbar(-1:15:89,b_avg(1,:),b_err(1,:),'ko-','LineWidth',1);	hold on
    errorbar(-1:15:89,b_avg(2,:),b_err(2,:),'bo-','LineWidth',1);
    errorbar(-1:15:89,b_avg(3,:),b_err(3,:),'mo-','LineWidth',1);	hold on
    errorbar(1:15:91,m_avg(1,:),m_err(1,:),'k*--','LineWidth',1);
    errorbar(1:15:91,m_avg(2,:),m_err(2,:),'b*--','LineWidth',1);	hold on
    errorbar(1:15:91,m_avg(3,:),m_err(3,:),'m*--','LineWidth',1);
    
    
    % session x slant
    con = reshape(mean(dat(:,1,:,:),3),[nsubj nslant]);
    V6A = reshape(mean(dat(:,2,:,:),3),[nsubj nslant]);
    AIP = reshape(mean(dat(:,3,:,:),3),[nsubj nslant]);
    c_avg = squeeze(mean(con,1));
    v_avg = squeeze(mean(V6A,1));
    a_avg = squeeze(mean(AIP,1));
    c_err = squeeze(sem(con,0,1));
    v_err = squeeze(sem(V6A,0,1));
    a_err = squeeze(sem(AIP,0,1));
    figure(s+2);
    hold off
    errorbar(-1:15:89,c_avg,c_err,'ko--','LineWidth',1);	hold on
    errorbar(0:15:90,v_avg,v_err,'bo--','LineWidth',1);
    errorbar(1:15:91,a_avg,a_err,'mo--','LineWidth',1);	hold on
end

if any(strcmpi(plot_effect,'vision x slant'))
    bino = reshape(mean(dat(:,:,1,:),2),[nsubj nslant]);
    mono = reshape(mean(dat(:,:,2,:),2),[nsubj nslant]);
    b_avg = squeeze(mean(bino,1));
    m_avg = squeeze(mean(mono,1));
    b_err = squeeze(sem(bino-mono,0,1));    % squeeze(sem(bino))
    m_err = squeeze(sem(bino-mono,0,1));    % squeeze(sem(mono))
    figure(s+3);
    hold off
    errorbar(-1:15:89,b_avg,b_err,'ro--','LineWidth',1);	hold on
    errorbar(1:15:91,m_avg,m_err,'go--','LineWidth',1);
end




% bm_avg = mean(m_bino-m_mono);
% bm_err = sem(m_bino-m_mono);
% figure(s+2);
% hold off
% errorbar(0:15:90,bm_avg,bm_err,'bo--','LineWidth',2)
% hold on
% line([-20 100],[0 0],'LineWidth',2,'Color','k','LineStyle','--');



