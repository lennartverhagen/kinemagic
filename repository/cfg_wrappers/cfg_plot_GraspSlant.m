function cfg_plot_GraspSlant
%--------------------------------------------------------------------------
% CFG_PLOT_GRASPSLANT is a wrapper function for kinematic data plotting for
% the GraspSlant experiment (TMS over vertex).
%
% Evaluate the code below to add the correct paths for this experiment:
% startup('Kinemagic');
% cdata('slant');
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------


%% Basics
%----------------------------------------
subjects	= {      '',      '',      '',      '',      '',...
                     '','subj07','subj08','subj09','subj10',...
               'subj11','subj12','subj13','subj14','subj15',...
               'subj16','subj17','subj18','subj19','subj20',...
               'subj21','subj22','subj23','subj24','subj25',...
               'subj26','subj27','subj28','subj29','subj30',...
               'subj31','subj32','subj33','subj34','subj35'};
sessions    = {'vertex','V6A','AIP'};

allsubj     = 8:35;
badstereo   = [6 15 16];
badvergence = [23 27];
bad2fast    = [9 17];
bad2slow    = [18 22];
bad2var     = [11 19 22];
badsubj     = unique([bad2slow bad2var]);
goodsubj    = setdiff(allsubj,badsubj);

subjsel	= goodsubj;
%subjsel	= allsubj;
sesssel	= 1;


%% Set task configuration
%----------------------------------------
% basic configuration
cfg             = [];
cfg.save        = 'no';
cfg.load        = 'no';
cfg.exp         = 'GraspSlant';
cfg.dataset     = 'MASKpchip_FILTlp15_MOV_ANA';
if ispc
    cfg.dir.root	= fullfile('M:','Data','GraspSlant');
else
    cfg.dir.root	= fullfile(filesep,'home','action','lenver','Data','GraspSlant');
end
cfg.dir.ana     = 'ana';
cfg.dir.report  = 'report';
cfg.subj        = subjects(subjsel);
cfg.sess        = sessions(sesssel);

% parameter and time-series selection
%----------------------------------------
% plots for slant paper
cfg.plot.param  = {'rt','mt','tl','rtpv','rtpga','mv','pv','pga'};
%cfg.plot.tseries = {'pos','vel','gripapt','gripori'};
% other
%cfg.plot.param  = {};
cfg.plot.tseries = 'no';

% RUN
%cfg.plot.fact   = 'run';
% VISION
%cfg.plot.fact   = 'vision';
% SLANT
%cfg.plot.fact   = 'slant';
% VISION x SLANT
cfg.plot.fact   = 'vision x slant';
% CONTRAST and main effect subtraction
%cfg.plot.contr  = {};
%cfg.plot.factrem= {};

% factors/condition selection
%----------------------------------------
cfg.plot.lvl.run        = [1 2];    % run number
cfg.plot.lvl.vision     = [1 2];    % 1:bino, 2:mono
cfg.plot.lvl.slant      = 0:15:90;  % clockwise deviatation from the vertical plane
cfg.plot.lvl.verthori	= [1 2];    % 1:vert, 2:hori
cfg.plot.lvl.congr      = [1 2];    % 1:congr, 2:incongr
cfg.plot.lvltick.run	= {'1','2'};
cfg.plot.lvltick.vision	= {'bino','mono','left'};
cfg.plot.lvltick.slant	= {'0','15','30','45','60','75','90'};
cfg.plot.lvltick.verthori	= {'vert','hori'};
cfg.plot.lvltick.congr	= {'congr','incongr'};

% factors/condition selection
%----------------------------------------
cfg.plot.factexcl       = {};
cfg.plot.lvlexcl.run    = 1;
cfg.plot.lvlexcl.vision = 3;
cfg.plot.lvlexcl.slant  = 3;

% end-point-error settings
%----------------------------------------
cfg.plot.epe.marker	= {'t','i','meanti','MCP'};
cfg.plot.epe.axis  	= 'yz';
cfg.plot.epe.param 	= 'area';
%cfg.plot.epe.axis  	= 'xyz';
%cfg.plot.epe.param 	= 'volume';

% confidence interval settings
%----------------------------------------
%cfg.plot.CI         = 'yes';
cfg.plot.CI.param    = {'rt','mt'};

% time series settings
%----------------------------------------
trlrejlst	= {'apriori','rt','mt','tl','mv','pv','pga','ga_E','dgo_E_Re'};
movname     = 'movement';
%cfg.tseries.pos.marker      = {'j','t','i','mti'};
%cfg.tseries.pos.axis        = {'x','y','z'};
cfg.tseries.pos.marker      = {'t','i'};
cfg.tseries.pos.axis        = {'y','z'};
cfg.tseries.pos.trlrej      = trlrejlst;
cfg.tseries.pos.movname     = movname;
%cfg.tseries.pos.basecorr.win        = @(x) [x(2)+0.1 x(2)+0.3];
%cfg.tseries.pos.basecorr.val        = 0.06;
%cfg.tseries.pos.basecorr.descrip    = 'nanmedian';
%cfg.tseries.pos.basecorr.corrflg    = 'sess';
cfg.tseries.vel.marker      = {'j','t','i','mti'};
cfg.tseries.vel.axis        = {'yz'};
cfg.tseries.vel.trlrej      = trlrejlst;
cfg.tseries.vel.movname     = movname;
cfg.tseries.gripapt.marker	= {'gti','gjt','gji'};
cfg.tseries.gripapt.axis	= 'xyz';
cfg.tseries.gripapt.trlrej	= trlrejlst;
cfg.tseries.gripapt.movname = movname;
% cfg.tseries.gripapt.basecorr.win        = @(x) [x(2)+0.1 x(2)+0.3];
% cfg.tseries.gripapt.basecorr.val        = 0.06;
% cfg.tseries.gripapt.basecorr.descrip    = 'nanmedian';
% cfg.tseries.gripapt.basecorr.corrflg    = 'sess';
cfg.tseries.gripaptvel.marker	= {'gti','gjt','gji'};
cfg.tseries.gripaptvel.axis     = 'yz';
cfg.tseries.gripaptvel.trlrej	= trlrejlst;
cfg.tseries.gripaptvel.movname	= movname;
cfg.tseries.gripori.marker	= {'gti','gjt','gji'};
cfg.tseries.gripori.axis	= 'yz';
cfg.tseries.gripori.trlrej	= trlrejlst;
cfg.tseries.gripori.movname = movname;

% process timeseries before plotting
%----------------------------------------
cfg.tseries.proc.task               = {'grip'};
cfg.tseries.proc.grip.marker        = {{'t','i'},{'mti','o'},{'j','t'},{'j','i'}};
cfg.tseries.proc.grip.label     	= {'gti','d2t','gjt','gji'};
cfg.tseries.proc.grip.apt.axis    	= {'xyz','yz'};
cfg.tseries.proc.grip.ori.axis    	= {'yz','z'};
cfg.tseries.proc.grip.ori.dir   	= 'cw';
cfg.tseries.proc.grip.ori.offset	= 90;

% trial rejection
%----------------------------------------
cfg.plot.trlrej.listwise = trlrejlst;
cfg.plot.trlrej.pairwise = 'no';
cfg.plot.trlrej.dataset  = 'MASKpchip_FILTlp15_MOV_ANA_IDXR';
cfg.plot.trlrej.mt     	 = {'mt','mt_T','mt_A'};
%cfg.plot.reptcontrib	= 10;   % main effects
%cfg.plot.reptcontrib	= 6;    % interaction
cfg.plot.reptcontrib	= 3;    % test
cfg.plot.subjcontrib	= 6;


% report and statistics
%----------------------------------------
cfg.paramreport.write   = {'groupuni'};
cfg.paramreport.vars    = {'all','-trlbeg','-trlend','-trloff','-RT','-MT','-GT','-HRT','-early','-late','-still','-pause','-t_*'};
cfg.paramreport.reject  = {'all'};
cfg.paramreport.exclude = 'cases';
%cfg.paramreport         = 'no';
cfg.stats               = 'no';


% selecting plot types
%----------------------------------------
if ~isfalse(cfg.plot.tseries),  cfg.tseriesplot = 'yes';    end
if ~isempty(regexp(cfg.plot.fact,'slant','once'))
    %cfg.paramplot.type      = 'errorbar';
    cfg.paramplot.type      = 'errorcloud';
else
    cfg.paramplot.type      = 'bar';
end
% descriptives: mean, median, logmean, meanlog, logmedian, medianlog
cfg.plot.descrip.group	= 'mean';
cfg.plot.descrip.subj	= 'mean';
%cfg.plot.descrip.subj	= 'var';
%cfg.plot.descrip.subj	= 'median';
%cfg.plot.descrip.subj	= 'meanlog';
%cfg.plot.descrip.subj	= 'logmean';
%cfg.plot.descrip.subj	= 'logmedian';
cfg.plot.descrip.err	= 'sem';
cfg.plot.descrip.contrerr= 'no';


% plotting configuration
%----------------------------------------
% Colors
cfg.plot.flcl	= str2rgb('rgh');
%cfg.plot.lncl	= jetpack(7);
%cfg.plot.lnstyle	= {'-','- ','-'};


% execute plotting wrapper
%----------------------------------------
clc;
km_plot(cfg);
