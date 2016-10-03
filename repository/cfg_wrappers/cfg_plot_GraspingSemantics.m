function cfg_plot_GraspingSemantics
%--------------------------------------------------------------------------
% CFG_PLOT_GRASPINGSEMANTICS is a wrapper function for plotting the data
% of the Grasping Semantics experiment
%
% Evaluate the code below to add the correct paths for this experiment:
% startup('Kinemagic');
% cdata('tmseeg');
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------


%% Basics
%----------------------------------------
subjects	= {'subj01','subj02','subj03','subj04','subj05',...
               'subj06','subj07','subj08','subj09','subj10',...
               'subj11','subj12','subj13','subj14','subj15',...
               'subj16','subj17','subj18','subj19','subj20',...
               'pilot99'};

subjsel = 1:15;

%% Set task configuration
%----------------------------------------
% basic configuration
cfg             = [];
cfg.save        = 'no';
cfg.load        = 'no';
cfg.exp         = 'GraspingSemantics';
cfg.dataset     = 'MASKpchip_FILTlp15_MOV_ANA';
cfg.dir.root	= fullfile(filesep,'home','action','lenver','Data','GraspingSemantics');
cfg.dir.ana     = 'ana';
cfg.subj        = subjects(subjsel);

% parameter and time-series selection
%----------------------------------------

%cfg.plot.param = {'rt','mt','td','mv','pv','rtpv','pa','rtpa','pd','rtpd','mga','rtmga','ega','dega','epe'};
%cfg.plot.param = {'rt','mt','td','mv','pv','rtpv','mga','rtmga','ega','dega'};
cfg.plot.param = {'rt','mt'};
%cfg.plot.param = {'pa','rtmga'};
%cfg.plot.param  = {'epe'};
cfg.plot.tseries = 'no';

% factors/condition selection
%----------------------------------------
% VISION
%cfg.plot.fact   = 'vision';
% SIZE
%cfg.plot.fact   = 'size';
% COLOR
%cfg.plot.fact   = 'color';
% WORD
%cfg.plot.fact   = 'word';
% VISION X SIZE
%cfg.plot.fact   = 'vision x size';
% SIZE X WORD
%cfg.plot.fact   = 'size x word';
%cfg.plot.fact   = 'word x size';
% VISION X WORD
%cfg.plot.fact   = 'vision x word';
% COLOR X WORD
%cfg.plot.fact   = 'color x word';
%cfg.plot.fact   = 'word x color';
% VISION X SIZE X WORD
%cfg.plot.fact   = 'vision x word x size';
% SIZE X WORD X RUN
%cfg.plot.fact   = 'word x size x run';
% COLOR X WORD X RUN
%cfg.plot.fact   = 'word x color x run';
% CONGRUENCY
cfg.plot.fact   = 'size x colorcongr';

% factors/condition definition
%----------------------------------------
cfg.plot.lvl.run        = [1:4];        % run number
cfg.plot.lvl.vision     = [1 2];        % 1:bino, 2:right, 3:left
cfg.plot.lvl.size       = [1 2];        % 1:small, 2:large, 3:middle
cfg.plot.lvl.color      = [1 2];        % 1:yellow, 2:blue, 3:white
cfg.plot.lvl.word       = [1 2];        % 1:klein, 2:groot, 3:geel, 4:blauw
cfg.plot.lvl.sizecongr  = [1 2];        % 1:congruent, 2:incongruent
cfg.plot.lvl.colorcongr = [1 2];        % 1:congruent, 2:incongruent

cfg.plot.lvltick.run	= {'1','2','3','4'};
cfg.plot.lvltick.vision	= {'bino','mono','left'};
cfg.plot.lvltick.size   = {'small','large','middle'};
cfg.plot.lvltick.color  = {'yellow','blue','white'};
cfg.plot.lvltick.word   = {'KLEIN','GROOT','GEEL','BLAUW'};
cfg.plot.lvltick.sizecongr  = {'congr','incongr'};
cfg.plot.lvltick.colorcongr = {'congr','incongr'};

% factors/condition selection
%----------------------------------------
cfg.plot.factexcl       = {'size','vision','word'};
cfg.plot.lvlexcl.run    = [1 2];
cfg.plot.lvlexcl.vision = 3;
cfg.plot.lvlexcl.size   = 3;
cfg.plot.lvlexcl.color  = 3;
cfg.plot.lvlexcl.word   = [1 2];

% end-point-error settings
%----------------------------------------
cfg.plot.epe.marker	= {'thumb','index','MPJ'};
cfg.plot.epe.axis  	= 'xyz';
cfg.plot.epe.param 	= 'volume';

% confidence interval settings
%----------------------------------------
cfg.plot.CI         = 'yes';

% trial rejection
%----------------------------------------
cfg.plot.trlrej.dataset  = 'MASKpchip_FILTlp15_MOV_ANA_IDXR';
%cfg.plot.trlrej.mode    = 'apriori';
%cfg.plot.reptcontrib	= 10;   % main effects
%cfg.plot.reptcontrib	= 6;    % interaction
cfg.plot.reptcontrib	= 2;    % test
cfg.plot.subjcontrib	= 1;

% statistics
%----------------------------------------
%cfg.paramreport.write   = {'groupuni'};
%cfg.paramreport.vars    = {'all','-trlbeg','-trlend','-trloff','-RT','-TMT','-early','-late','-still','-pause','-t_*'};
%cfg.paramreport.reject  = {'reject'};
cfg.paramreport         = 'no';
cfg.stats               = 'no';

% selecting plot types
%----------------------------------------
if ~isempty(regexp(cfg.plot.fact,'slant','once'))
    cfg.paramplot.type      = 'errorbar';
else
    cfg.paramplot.type      = 'bar';
end
% descriptives: mean, median, logmean, meanlog, logmedian, medianlog
cfg.paramplot.descrip.group	= 'mean';
cfg.paramplot.descrip.subj	= 'mean';
%cfg.paramplot.descrip.subj	= 'median';
%cfg.paramplot.descrip.subj	= 'meanlog';
%cfg.paramplot.descrip.subj	= 'logmean';
%cfg.paramplot.descrip.subj	= 'logmedian';
cfg.paramplot.descrip.err	= 'sem';

% plotting configuration
%----------------------------------------
% Colors
Cblack      = [ 0   0   0 ];
Cgrey       = [102 102 102];
Cred        = [255  0   0 ];    % bino
Cgreen      = [ 0  255  0 ];    % mono
Cblue       = [ 0   0  255];
Corange     = [255 153  0 ];
Cmagenta	= [204  0  204];
Ccyan       = [ 0  204 204];
%cfg.plot.clspec     = [Cred; Cgreen]/255;          % vision
cfg.plot.clspec     = [Cred; Cgreen; Corange; Cblue]/255;  % word
%cfg.plot.clspec     = [Cred; Cgreen; Cred; Cgreen]/255;    % vision x size
%cfg.plot.clspec     = [Cmagenta; Ccyan; Cmagenta; Ccyan; Cmagenta; Ccyan; Cmagenta; Ccyan]/255;     % size x word
%cfg.plot.clspec     = [Corange; Cblue; Corange; Cblue; Corange; Cblue; Corange; Cblue]/255;         % color x word
%cfg.plot.lnstyle	= {'-','- ','-'};


% execute plotting wrapper
%----------------------------------------
km_plot(cfg);
