%--------------------------------------------------------------------------
% This script gives an overview of behavioral outcome measures.
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------

% what subject to analyse
subjsel    	= [1:15];
nsubj       = length(subjsel);

% settings
dataset = 'FILTlp15_MOV_ANA';
dir_proc = fullfile(filesep,'home','action','ankmur','CommPointing','Data-analysis','KIN','ana');
subjects        = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15'};

% code to take the right reactionvaliable
code = [11 12 21 22];
meanRT_P2 = nan(nsubj,length(code));
meanMT_P2 = meanRT_P2 ; meanHT_P2 = meanRT_P2; meanBT_P2  = meanRT_P2 ;

% loop over subjects
for s = 1:nsubj
    
    % select the correct file
    subj = subjects{subjsel(s)};
    fname_cfg = fullfile(dir_proc,sprintf('%s_%s_CFG',subj,dataset));
    fname_idxr = fullfile(dir_proc,sprintf('%s_%s_IDXR',subj,dataset));
    
    % load the cfg files with the trl matrix
    cfg = load(fname_cfg);
    cfg = cfg.cfg;
    
    % specify the conditions
    vars = cfg.vars;
    [~,idxcode] = ismember({'action','addressee'},vars);
    [~,trlvar] = ismember({'ht_P1','rt_P2','mt_P2','ht_P2','bt_P2'},vars);
    
    % code for sorting the conditions action and addressee
    cond_code = (cfg.trl(:,idxcode));
    cond_code = 10*cond_code(:,1) + cond_code(:,2);
    
    % select movement variables
    trl = (cfg.trl(:,trlvar));
    
    % add extra columns in trl and vars
    trl = [trl nan(size(trl,1),1)];
    
    % shift median split coding one down
    trl(:,6) = [NaN; trl(1:end-1,6)];
    
    % load the trial rejection indices
    idxr = load(fname_idxr);
    idxr = idxr.cfg;
    
    % specify which trials to reject
    idxr = idxr.corr_P1 | idxr.corr_P2 | idxr.toolate | idxr.first4trials | idxr.rt_1 | idxr.nomov | idxr.outlier | idxr.blockstart;
    idxr = any(idxr,2);
    
    % apply additional criteria to reject
    idxr = any(idxr | trl(:,1)<=0 | trl(:,2)<50 | trl(:,3)<50,2);
    
    % temporarilly replace reject trials with NaN
    tmptrl = trl;
    tmptrl(idxr,:) = NaN;
    
    for c = 1:4
        % find median split
        tmp = nanmedian(tmptrl(cond_code==code(c),3));
        % code trials lower and higher than median split
        trl(cond_code==code(c) & trl(:,3)<tmp,6) = 1;
        trl(cond_code==code(c) & trl(:,3)>=tmp,6) = 2;
    end
    
    % code for sorting the conditions action and median split
    [~,idxcode] = ismember({'action'},vars);
    cond_code = [cfg.trl(:,idxcode) trl(:,6)];
    cond_code = 10*cond_code(:,1) + cond_code(:,2);
    
    % replace reject trials with NaN
    trl(idxr,:) = NaN;
    
    for c = 1:4
        med_HT_P1(s,c) = nanmedian(trl(cond_code==code(c),1));
    end
  
end