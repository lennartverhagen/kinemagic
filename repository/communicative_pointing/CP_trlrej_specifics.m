%--------------------------------------------------------------------------
% This script gives the conditions in which 
% trials were rejected.
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

trl_rej_type = 1;   % show errors per target:sel_P1
%trl_rej_type = 2;    % show errors per condition:action,addresse

% when type = 1 then select one of the following
%trl_rej_spec = 2;   % After first trialrejection, before epe
trl_rej_spec = 1;  % After epe analysis, outliers rejected

% settings
dataset = 'FILTlp15_MOV_ANA';
dir_proc = fullfile(filesep,'home','action','ankmur','CommPointing','Data-analysis','KIN','ana');
%dir_proc = fullfile(filesep,'CommPointing','Data-analysis','KIN','ana');
subjects        = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15'};

commleft = 0;
commright = 0;
instrleft = 0;
instrright = 0;

if trl_rej_type == 1

    % loop over subjects
    for s = 1:nsubj
        
        code = [114 115 116 124 125 126 214 215 216 224 225 226]';
        mu = nan([length(code),nsubj]);
        ni = mu; nr = mu; mu2 = mu; ni2 = ni; nr2 = nr;
        su(s) = 0;  st(s) = 0; su2(s) = 0;  st2(s) = 0;
    
         % select the correct file
        subj = subjects{s};
        fname_cfg = fullfile(dir_proc,sprintf('%s_%s_CFG',subj,dataset));
        fname_idxr = fullfile(dir_proc,sprintf('%s_%s_IDXR',subj,dataset));
        
        % load the correct files
        cfg = load(fname_cfg);
        cfg = cfg.cfg;
        idxr = load(fname_idxr);
        idxr = idxr.cfg;
        
        %% HACK - gather specifics
        % save data
        fprintf('\n\nsubject: %s\n',subj);
        tmp = idxr.outlier | idxr.apriori | idxr.corr_P1 | idxr.corr_P2 | idxr.toolate | idxr.first2trials | idxr.rt_1 | idxr.nomov;
        fprintf('first2trials:%3.0f\t\n',sum(idxr.first2trials & tmp));
        tmp(idxr.first2trials) = 0;
        fprintf('errall:%3.0f\t\n',sum((idxr.corr_P2 | idxr.corr_P1) & tmp));
        fprintf('err_P2:%3.0f\t\n',sum(idxr.corr_P2& tmp));
        fprintf('err_P1:%3.0f\t\n',sum(idxr.corr_P1& tmp));
        fprintf('overlap:%3.0f\t\n',sum((idxr.corr_P2 | idxr.corr_P1) & (idxr.corr_P2 & idxr.corr_P1)));
        tmp(idxr.corr_P2 | idxr.corr_P1) = 0;
        fprintf('nomov:%3.0f\t\n',sum(idxr.nomov & tmp));
        tmp(idxr.nomov) = 0;
        fprintf('toolate:%3.0f\t\n',sum(idxr.toolate & tmp));
        tmp(idxr.toolate) = 0;
        fprintf('rt_1:%3.0f\t\n',sum(idxr.rt_1 & tmp));
        tmp(idxr.rt_1) = 0;
        fprintf('outlier:%3.0f\t\n',sum(idxr.outlier & tmp));
        tmp(idxr.outlier) = 0;
        fprintf('totrej:%3.0f\t\n',sum(idxr.apriori));
        fprintf('ntot:%3.0f\t\n',length(idxr.apriori));
        %END
       
        %matrix = [idxr.outlier,idxr.corr_P1,idxr.corr_P2,idxr.toolate,idxr.corr_P1,idxr.first2trials,idxr.rt_1,idxr.nomov];
        %[I,J,V] = find(matrix)
        % specify what rejection reason you want to evaluate
        idxr = idxr.outlier | idxr.apriori | idxr.corr_P1 | idxr.corr_P2 | idxr.toolate | idxr.first2trials | idxr.rt_1 | idxr.nomov;
        %idxr = idxr.nomov;
        idxr = any(idxr,2);
        
        % specify in which conditions the rejections are shown -->kan uit loop?
        vars = cfg.vars;
        [~,idxvar] = ismember({'action','addressee','sel_P1'},vars);
        %[~,idxvar] = ismember({'action','time_half','sel_P1'},vars);
        
        % first rejection, before epe
        if trl_rej_spec == 1
            
            % show conditionnumbers for all trials
            trl = (cfg.trl(:,idxvar));
            % merge to form a code for all trails
            trl = 100*trl(:,1) + 10*trl(:,2) + trl(:,3);
            
            % calc survivors and rejections
            for c = 1:length(code)
                % absolute value of trials that survived
                ni(c,s) = sum(~idxr(trl==code(c)));
                % absolute value of trials that were rejected
                nr(c,s) = sum(idxr(trl==code(c)));
                % percentage of trials that were rejected per code
                mu(c,s) = mean(idxr(trl==code(c)))*100;
                % calc sum of rejections
                su(s) = sum(su(s)) + nr(c,s);
                % which was the lowest trlnr in one condition
                
                % calc mean of survived trials per condition
                st(s) = sum(st(s)) + ni(c,s);
            end
            
            % calc mean survived trials per subject per condition
            st(s) = st(s)/length(code);
            sc(s) = min(ni(:,s));
            
            % print this in command window
            %fname = fullfile(filesep,'CommPointing','Data-analysis','samplesize','rejectioncrits.txt');
            %fid = fopen(fname,'w');
            %dlmwrite(fname,meanRTall,'-append','delimiter','\t');

            % print this in command window
            fprintf('\n\nsubject: %s\n',subj);
            fprintf('meansurvivors:%3.0f\t\n',st(s));
            fprintf('sumrejections:%3.0f\t\n',su(s));
            fprintf('min1condition:%3.0f\t\n\n',sc(s));
            fprintf('%3.0f\t%3.0f\t%3.0f\t%3.0f\n',[code ni(:,s) nr(:,s) mu(:,s)]');
            fprintf('%3.0f\t%3.0f\t%3.0f\t%3.0f\n',[code ni(:,s) nr(:,s) mu(:,s)]');
            
            
            
        elseif trl_rej_spec == 2
            
            % show how many trials were rejected in which condition
            trl = (cfg.trl(:,idxvar));
            trl = 100*trl(:,1) + 10*trl(:,2) + trl(:,3);
            for c = 1:length(code)
                % absolute value of trials that survived
                ni2(c,s) = sum(~idxr(trl==code(c)));
                % absolute value of trials that were rejected
                nr2(c,s) = sum(idxr(trl==code(c)));
                % percentage of trials that were rejected per code
                mu2(c,s) = mean(idxr(trl==code(c)))*100;
                % calc sum of rejections
                su2(s) = sum(su2(s)) + nr2(c,s);
                % TO DO: calc mean of survived trials per condition
                st2(s) = sum(st2(s)) + ni2(c,s);
            end
            
            % calc mean survived trials per subject per condition
            st2(s) = st2(s)/length(code);
            
            % print this in command window
            fprintf('\n\nsubject: %s\n',subj);
            %fprintf('\n\nmeansurvivors:%3.0f\t\n',st2(s));
            fprintf('\n\nsumrejections:%3.0f\t\n',su2(s));
            %fprintf('%3.0f\t%3.0f\t%3.0f\t%3.0f\n',[code ni2(:,s) nr2(:,s) mu2(:,s)]');
            
            commleft = commleft + nr2(1,s);
            commright = commleft + nr2(2,s);
            instrleft = commleft + nr2(3,s);
            instrright = commleft + nr2(4,s);
            
            disp(commleft);
            
        end;
        % store and retrieve variables based on trl_rej_spec
        if trl_rej_spec == 1
            save('temp.mat','mu','ni','nr');
        else
            load('temp.mat'); %delete('temp.mat');
            % trials rejected in end point error CI:
            % code-survivors1-survivors2-difference-difference in percentages
            for s = 1:nsubj
                tmp = [ni(:,s) - ni2(:,s), mu(:,s) - mu2(:,s)];
                subj = subjects{subjsel(s)};
                fprintf('\n\nsubject: %s\n',subj);
                fprintf('%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\n',[code ni(:,s) ni2(:,s) tmp ]');
            end
        end
    end
end

%% calculate errors in each condition: action and addressee
if trl_rej_type == 2

   mu = nan([length(code),nsubj]);
    ni = mu; nr = mu; mu2 = mu; ni2 = ni; nr2 = nr;
    su(s) = 0;  st(s) = 0; su2(s) = 0;  st2(s) = 0;
    
    % loop over subjects
    for s = 1:nsubj
        
    code = [11 12 21 22]';
   
    
         % select the correct file
        subj = subjects{subjsel(s)};
        fname_cfg = fullfile(dir_proc,sprintf('%s_%s_CFG',subj,dataset));
        fname_idxr = fullfile(dir_proc,sprintf('%s_%s_IDXR',subj,dataset));
        
        % load the correct files
        cfg = load(fname_cfg);
        cfg = cfg.cfg;
        idxr = load(fname_idxr);
        idxr = idxr.cfg;
        
        % specify what rejection reason you want to evaluate
        %idxr = idxr.outlier | idxr.toolate | idxr.corr_P1 | idxr.corr_P2 | idxr.first2trials |  idxr.nomov;
        idxr = idxr.corr_P1;
        idxr = any(idxr,2);
        
        % specify in which conditions the rejections are shown -->kan uit loop?
        vars = cfg.vars;
        [~,idxvar] = ismember({'action','addressee'},vars);
        
        % first rejection, before epe
        if trl_rej_spec == 1
            
            % show conditionnumbers for all trials
            trl = (cfg.trl(:,idxvar));
            % merge to form a code for all trails
            trl = 10*trl(:,1) + 1*trl(:,2);
            
            % calc survivors and rejections
            for c = 1:length(code)
                % absolute value of trials that survived
                ni(c,s) = sum(~idxr(trl==code(c)));
                % absolute value of trials that were rejected
                nr(c,s) = sum(idxr(trl==code(c)));
                % percentage of trials that were rejected per code
                mu(c,s) = mean(idxr(trl==code(c)))*100;
                % calc sum of rejections
                su(s) = sum(su(s)) + nr(c,s);
                % TO DO: calc mean of survived trials per condition
                st(s) = sum(st(s)) + ni(c,s);
            end
            
            % calc mean survived trials per subject per condition
            st(s) = st(s)/length(code);
            
            % print this in command window
            fprintf('\n\nsubject: %s\n',subj);
            fprintf('\n\nmeansurvivors:%3.0f\t\n',st(s));
            fprintf('\n\nsumrejections:%3.0f\t\n',su(s));
            fprintf('%3.0f\t%3.0f\t%3.0f\t%3.0f\n',[code ni(:,s) nr(:,s) mu(:,s)]');
            
            
            commleft = commleft + nr(1,s);
            commright = commright + nr(2,s);
            instrleft = instrleft + nr(3,s);
            instrright = instrright + nr(4,s);
            
            disp(nr);
            
            %save('temp.mat','mu','ni','nr');
            %mean_all(s,1) = st(s);
            %mean_err(s,1) = su(s);
        end
    end
end

return


nrpt = 1;
p = 1:nrpt;
for i = 1:nrpt
    %a = sqrt(squeeze(area_verthori(1,i,:)));
    %b = sqrt(squeeze(area_verthori(2,i,:)));
    %a = sqrt(squeeze(area_all(1,i,:)));
    %b = sqrt(squeeze(area_all(2,i,:)));
    a = nr(1,:)' + nr(3,:)'; %
    b = nr(2,:)' + nr(4,:)'; % 
    [H,p(i),CI,stats] = ttest(a,b);
end


totaltrials = 240*13;

CL = (commleft/totaltrials)*100;
CR = (commright/totaltrials)*100;
IL = (instrleft/totaltrials)*100;
IR = (instrright/totaltrials)*100;
totalC = CL + CR;
totalI = IL + IR;
totalL = CL + IL;
totalR = CR + IR;

% mean survivors
mean_surv = mean(st)* 4;
perc_surv = (mean_surv/240)*100;

% mean rejections
mean_rej = mean(su);
perc_rej = (mean_rej/240)*100;

% P1 total errors
sum_sur = [57 54 54 56 48 55 52 52 55 56 55 55 48];
sum_sur = sum_sur * 4;
mean_sur = mean(sum_sur);
perc_sur = (mean_sur/240)*100;


% calc mean errors

% P1 total errors
sum_err_P1 = [2 5 3 9 21 5 23 12 9 5 12 9 39];
mean_err_P1 = mean(sum_err_P1);
mean_err_P1 = mean_err_P1;
perc_err_P1 = (mean_err_P1/240)*100;

% Too late
perc_err_late = (21/240)*100;

% P2 total error
sum_err_P2 = [10 15 6 13 15 10 13 10 10 9 13 16 42];
mean_err_P2 = mean(sum_err_P2);
perc_err_P2 = (mean_err_P2/240)*100;

% outliers end-points
sum_outlier = [29 35 35 19 55 36 36 45 27 23 46 30 28];
mean_out = mean(sum_outlier);
corr_trial = 240 - (0.108*240);
perc_out = (mean_out/corr_trial)*100;


