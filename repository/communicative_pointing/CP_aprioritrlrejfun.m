function idx = CP_aprioritrlrejfun(cfg,~,idx)
%--------------------------------------------------------------------------
% This function selects trials to be rejected based on
% behavioural parameters
%
% See also KM_TRLREJ
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------

% load data
dir_proc = fullfile(filesep,'CommPointing','Data-analysis','KIN','ana');
fname_data = fullfile(cfg.dir.proc,sprintf('%s%s%s_DATA',cfg.subj,cfg.sess,cfg.dataset));
load(fname_data);

% get trldata and variable names
vars = cfg.vars;
trl = cfg.trl;

% reject trials a priori
ntrl = size(trl,1);
idx.apriori = false(ntrl,1);

% exclude first trials of all blocks
action = trl(:,strcmp(vars,'action'));
addressee = trl(:,strcmp(vars,'addressee'));
comb = action*10 + addressee;
blocks = [logic2idx(comb==11); logic2idx(comb==12); logic2idx(comb==21); logic2idx(comb==22)];
idx.blockstart = false(ntrl,1);
idx.blockstart(sort(blocks(:,1))) = true;

% exclude first 2 trials
if isfield(idx,'first4trials'), idx = rmfield(idx,'first4trials'); end
trlnr = trl(:,strcmp(vars,'trl_nr'));
idx.first2trials = trlnr<3;

% exlude blocks
%idx.first4blocks = ;

% exclude incorrect trials
corr_P1 = trl(:,strcmp(vars,'corr_P1'));
corr_P2 = trl(:,strcmp(vars,'corr_P2'));
idx.corr_P1 = corr_P1<1 | corr_P1>1;
idx.corr_P2 = corr_P2<1 | corr_P2>1;

% exclude too late trials
rt_P1 = trl(:,strcmp(vars,'rt_P1'));
rt_P2 = trl(:,strcmp(vars,'rt_P2'));
idx.toolate = rt_P1<0 | rt_P2<0;

% exclude trials without (correctly identified) movements
idx.nomov = cellfun(@(x) any(isnan(x(:))),data.movement)';

% exclude trials lying out of the confidence ellipse
if isfield(data,'eposerr')
%     fname = 'sel_P1';
%     sname = 'i';
%     idx_out_zx = data.eposerr.(fname).(sname).zx.idx_out;
%     idx_out_xy = data.eposerr.(fname).(sname).xy.idx_out;
%     goal = trl(:,strcmp(vars,fname));
%     idx.outlier = false(size(goal));
%     idx.outlier(goal==4) = (idx_out_zx{1}~=0 | idx_out_xy{1}~=0);
%     idx.outlier(goal==5) = (idx_out_zx{2}~=0 | idx_out_xy{2}~=0);
%     idx.outlier(goal==6) = (idx_out_zx{3}~=0 | idx_out_xy{3}~=0);
    fname = 'actionxsel_P1';
    sname = 'i';
    idx_out_zx = data.eposerr.(fname).(sname).zx.idx_out;
    idx_out_xy = data.eposerr.(fname).(sname).xy.idx_out;
    goal = trl(:,strcmp(vars,'sel_P1'));
    action = trl(:,strcmp(vars,'action'));
    idx.outlier = false(size(goal));
    idx.outlier(goal==4 & action==1) = (idx_out_zx{1}~=0 | idx_out_xy{1}~=0);
    idx.outlier(goal==4 & action==2) = (idx_out_zx{2}~=0 | idx_out_xy{2}~=0);
    idx.outlier(goal==5 & action==1) = (idx_out_zx{3}~=0 | idx_out_xy{3}~=0);
    idx.outlier(goal==5 & action==2) = (idx_out_zx{4}~=0 | idx_out_xy{4}~=0);
    idx.outlier(goal==6 & action==1) = (idx_out_zx{5}~=0 | idx_out_xy{5}~=0);
    idx.outlier(goal==6 & action==2) = (idx_out_zx{6}~=0 | idx_out_xy{6}~=0);
else
    idx.outlier = idx.zeros;
end

% combine exclusion criteria
idx.apriori(idx.corr_P1 | idx.corr_P2 | idx.first2trials | idx.toolate | idx.nomov | idx.outlier) = true;
%idx.apriori(idx.apriori | idx.toolate) = true;