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
subjsel    	= [1:12 14];
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
    
    % load the trial rejection indices
    idxr = load(fname_idxr);
    idxr = idxr.cfg;
    
    % specify which trials to reject
    idxr = idxr.corr_P1 | idxr.corr_P2 | idxr.toolate | idxr.first4trials | idxr.rt_1 | idxr.nomov | idxr.outlier;
    idxr = any(idxr,2);
    
    % apply additional criteria to reject
    idxr = any(idxr | trl(:,1)<=0 | trl(:,2)<50 | trl(:,3)<50,2);
    
    % replace reject trials with NaN
    trl(idxr,:) = NaN;
    
    % initiate matrixes
    medHT_P1 = nan(nsubj,length(code));
    medRT_P2 = nan(nsubj,length(code));
    medMT_P2 = nan(nsubj,length(code));
    
    %% relative holding time P1(left-right) to relative RT P2(left-right)
    % separate according to different conditions
    for c = 1:length(code)
        % median holdingtime communicator
        medHT_P1(s,c) = nanmedian(trl(cond_code==code(c),1));
        % median reactiontime addressees
        medRT_P2(s,c) = nanmedian(trl(cond_code==code(c),2));
        medMT_P2(s,c) = nanmedian(trl(cond_code==code(c),3));
        medHT_P2(s,c) = nanmedian(trl(cond_code==code(c),4));
        medBT_P2(s,c) = nanmedian(trl(cond_code==code(c),5));
    end
    
    % separate according to different conditions
    for c = 1:length(code)
        % median holdingtime communicator
        medHT_P1(s,c) = nanmedian(trl(cond_code==code(c),1));
        % median reactiontime addressees
        meanRT_P2(s,c) = nanmean(trl(cond_code==code(c),2));
        meanMT_P2(s,c) = nanmean(trl(cond_code==code(c),3));
        meanHT_P2(s,c) = nanmean(trl(cond_code==code(c),4));
        meanBT_P2(s,c) = nanmean(trl(cond_code==code(c),5));
    end
    
    % calc relative median HT_P1 Left-Right
    rel_HT_Comm_P1(s) = medHT_P1(s,1) - medHT_P1(s,2);
    rel_HT_Instr_P1(s) = medHT_P1(s,3) - medHT_P1(s,4);
    % calc relative median RT_P2 Left-Right
    rel_RT_Comm_P2(s) = medRT_P2(s,1) - medRT_P2(s,2);
    rel_RT_Instr_P2(s) = medRT_P2(s,3) - medRT_P2(s,4);
    % calc relative median MT_P2 Left-Right
    rel_MT_Comm_P2(s) = medMT_P2(s,1) - medMT_P2(s,2);
    rel_MT_Instr_P2(s) = medMT_P2(s,3) - medMT_P2(s,4);
    % calc mean Reaction time & Movement time
    meanRTMT_Comm_P2(s) = (rel_RT_Comm_P2(s) + rel_MT_Comm_P2(s))/2;
    meanRTMT_Instr_P2(s) = (rel_RT_Instr_P2(s) + rel_MT_Instr_P2(s))/2;
    
    
    %% make two groups for HT_P1 and RT_P2 of all trials
    
     % initiate matrixes
    trl_HT_P1 = nan(nsubj,length(code),60);
    trl_RT_P2 = nan(nsubj,length(code),60);
    lowRT_P2 = nan(nsubj,length(code),60);
    highRT_P2 = nan(nsubj,length(code),60);
   
    % reaction times all trials rt_1 according to condition
    tmp = trl(cond_code==code(c),1);
    trl_HT_P1(s,c,1:length(tmp)) = tmp;
    tmp = trl(cond_code==code(c),2);
    trl_RT_P2(s,c,1:length(tmp)) = tmp;
    
    % sort absolute values
    [sortedValues,idx] = sort(trl_RT_P2(s,c,:));
    
    % initiate matrixes for low and high values
    highValues = sortedValues;
    lowValues = sortedValues;
    
    % take out high values for LowRT
    mu = highValues > medRT_P2(s,c);
    ind_low = idx(highValues > medRT_P2(s,c));
    lowValues(mu) = [];
    lowRT_P2(s,c,1:length(lowValues)) = lowValues;
    
    % take out low values for HighRT
    mu = lowValues <= medRT_P2(s);
    highValues(mu) = [];
    HighRT_P2(s,c,1:length(highValues)) = highValues;
    
    % clear tmp variables
    %clearvars n sortedValues;
    %clearvars n mu;
    
    
    
end

figure; hold off;
[R,P] = plot_corr(rel_RT_Comm_P2,rel_HT_Comm_P1,'r',false)
hold on
[R,P] = plot_corr(rel_RT_Instr_P2,rel_HT_Instr_P1,'b',true)
