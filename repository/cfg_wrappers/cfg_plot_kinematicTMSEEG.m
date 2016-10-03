function cfg_plot_kinematicTMSEEG
%--------------------------------------------------------------------------
% CFG_PLOT_KINEMATICTMSEEG is a wrapper function for EEG-TMS kinematic data
% plotting for the TMS-EEG Grasping experiment (TMS over V6A, AIP and
% vertex).
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
subjects	= {      '',      '',      '',      '',      '',...
                     '','subj07','subj08','subj09','subj10',...
               'subj11','subj12','subj13','subj14','subj15',...
               'subj16','subj17','subj18','subj19','subj20',...
               'subj21','subj22','subj23','subj24','subj25',...
               'subj26','subj27','subj28','subj29','subj30',...
               'subj31','subj32','subj33','subj34','subj35'};
sessions    = {'vertex','V6A','AIP'};

allsubj     = 8:35;
badcapfit   = [8 10];
badstereo   = [6 15 16];
badvergence = [23 27];
%nogoodERD   = [10 13 22 23];
%badEEGart    = [4 9 18 31 32];
badmissrun  = [8 35];
badfewtrls  = 35;
bad2fast    = [9 17];
bad2slow    = [18 22];
badEPE      = [14 24];
badsubj     = unique([badcapfit badmissrun badfewtrls bad2slow]);
goodsubj    = setdiff(allsubj,badsubj);

subjsel	= goodsubj;
%subjsel	= allsubj;
sesssel	= [1 3 2];


%% Set task configuration
%----------------------------------------
% basic configuration
cfg             = [];
cfg.save        = 'no';
cfg.load        = 'no';
cfg.exp         = 'kinematicTMSEEG';
cfg.dataset     = 'MASKpchip_FILTlp15_MOV_ANA';
if ispc
    cfg.dir.root	= fullfile('M:','Data','TMSEEG');
else
    cfg.dir.root	= fullfile(filesep,'home','action','lenver','Data','TMSEEG');
end
cfg.subj        = subjects(subjsel);
cfg.sess        = sessions(sesssel);

% parameter and time-series selection
%----------------------------------------
%cfg.plot.param	 = {'rt','mt','tt','rtt','ttt','rttt','itt','ritt','itdiff','ritdiff','gt','rgt','td','tcd','gcd','mv','mtv','mgv','pv','tpv','rtpv','pa','tpa','rtpa','pd','tpd','rtpd','nvc','tnvc','gnvc','pgv','npgv','tgv','mggv','gcgd','mga','tmga','rtmga','tga','dtga','ega','dega','tgo','dtgo','adtgo','ego','dego','adego'};	% 'epos'
%cfg.plot.param  = {'epe'};
%cfg.plot.param  = {'CI'};
% VERTEX
%cfg.plot.param  = {'mt','td','mv','pv','rtpv','mga','rtmga','ega','ego'};   % effect of slant
%cfg.plot.param  = {'rt','mt','tt','gt','td','tcd','gcd','mv','mtv','mgv'};  % effect of vision
%cfg.plot.param  = {'td','mga','tga','tgv','dtgo'};    % effect of vision x slant
% V6A
%cfg.plot.param  = {'rtt','rgt'};    % effect of vision
%cfg.plot.param  = {'mv','mtv','pv','rtpv','mga','tga'};    % effect of vision
% AIP
%cfg.plot.param  = {'rt','tt','tcd'};    % main effect of stimulation and effect of vision
%cfg.plot.param  = {'rgt','gcd','mgv','mga','tga','tgv','dtgo'}; % size x vision x slant
% plots for AIP paper
%cfg.plot.param  = {'rt','mt','tt','rgt','td','tcd','gcd','mv','mtv','mgv','pv','rtpv','mga','rtmga','tga','dtgo','tgv'};
%cfg.plot.param  = {'rt','mt','tt','rgt','td_meanti','tcd_meanti','gcd_meanti','mv_meanti','mtv_meanti','mgv_meanti','pv_meanti','rtpv_meanti','mga','rtmga','tga','dtgo','tgv'};
%cfg.plot.param  = {'tga','gcgd','dtgo'};
% plots for slant paper
cfg.plot.param  = {'rt','mt','tt','gt','rtt','rgt','td','tcd','gcd','mv','mtv','mgv','pv','tpv','rtpv','mga','tmga','rtmga','tga','dtgo'};
% other
%cfg.plot.param  = {};
cfg.plot.tseries = 'no';
%cfg.plot.tseries = {'pos','vel','gripapt','gripaptvel','gripori'};
%cfg.plot.tseries = {'gripapt'};

% RUN
%cfg.plot.fact   = 'run';
% TMS TIME
%cfg.plot.fact   = 'TMSel';
% TMS SITE
%cfg.plot.fact   = 'site';
% VISION
%cfg.plot.fact   = 'vision';
% VISION X SITE
%cfg.plot.fact   = 'vision x site';
% SLANT
%cfg.plot.fact   = 'slant';
% SITE X SLANT
%cfg.plot.fact   = 'site x slant';
%cfg.plot.fact   = 'verthori x site';
cfg.plot.fact   = 'site x verthori';
% VISION x SLANT
%cfg.plot.fact   = 'slant x vision';
%cfg.plot.fact   = 'vision x slant';
%cfg.plot.fact   = 'vision x verthori';
% SITE X VISION x SLANT
%cfg.plot.fact   = 'site x vision x slant';
%cfg.plot.fact   = 'site x vision x verthori';
%cfg.plot.fact   = 'vision x verthori x site';
% SITE X VISION x SLANT x RUN
%cfg.plot.fact   = 'vision x verthori x site x run';
% SITE X CONGR
%cfg.plot.addvar.newvar = 'congr';
%cfg.plot.addvar.oldvar = {'vision','verthori'};
%cfg.plot.addvar.fun = '2-((tmp(:,1)==1 & tmp(:,2)==1) | (tmp(:,1)==2 & tmp(:,2)==2))';
%cfg.plot.fact   = 'congr x site';
% CONTRAST and main effect subtraction
cfg.plot.contr  = {};
%cfg.plot.contr  = {'site'};
%cfg.plot.contr  = {'vision','verthori'};
%cfg.plot.factrem= {};
%cfg.plot.factrem= {'verthori','site'};

% factors/condition selection
%----------------------------------------
cfg.plot.lvl.run        = [1 2];    % run number
cfg.plot.lvl.TMSel      = [1 2];    % 1:early, 2:late
cfg.plot.lvl.vision     = [1 2];    % 1:bino, 2:mono
cfg.plot.lvl.slant      = 0:15:90;  % clockwise deviatation from the vertical plane
cfg.plot.lvl.verthori	= [1 2];    % 1:vert, 2:hori
cfg.plot.lvl.congr      = [1 2];    % 1:congr, 2:incongr
cfg.plot.lvl.site       = [1 2 3];  % 1:vertex, 2:V6A, 3:AIP
cfg.plot.lvl.site       = cfg.plot.lvl.site(sesssel);
cfg.plot.lvltick.run	= {'1','2'};
cfg.plot.lvltick.TMSel	= {'early','late'};
cfg.plot.lvltick.site	= {'vertex','V6A','AIP'};
cfg.plot.lvltick.vision	= {'bino','mono','left'};
cfg.plot.lvltick.slant	= {'0','15','30','45','60','75','90'};
cfg.plot.lvltick.verthori	= {'vert','hori'};
cfg.plot.lvltick.congr	= {'congr','incongr'};

% factors/condition selection
%----------------------------------------
cfg.plot.factexcl       = {};
cfg.plot.lvlexcl.run    = 1;
cfg.plot.lvlexcl.TMSel  = 2;
cfg.plot.lvlexcl.vision = 3;
cfg.plot.lvlexcl.slant  = 3;
cfg.plot.lvlexcl.site   = 2;

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
cfg.plot.CI.param    = {'rt','mt','td','gcd','tgv','mga','tga','dtga','tgo','dtgo','ega','dega','ego','dego'};

% time series settings
%----------------------------------------
movname = 'movement';
cfg.tseries.pos.marker      = {'MCP','t','i'};
cfg.tseries.pos.axis        = {'x','y','z'};
%cfg.tseries.pos.marker      = {'meanti'};
%cfg.tseries.pos.axis        = {'y','z'};
cfg.tseries.pos.trlrej      = {'rt','mt','mga','rtmga','tga'};
cfg.tseries.pos.movname     = movname;
cfg.tseries.vel.marker      = {'MCP','t','i'};
cfg.tseries.vel.axis        = {'yz'};
cfg.tseries.vel.trlrej      = {'rt','mt','mga','rtmga','tga'};
cfg.tseries.vel.movname     = movname;
cfg.tseries.gripapt.marker	= {'gripti','gripMCPt','gripMCPi'};
cfg.tseries.gripapt.axis	= 'xyz';
cfg.tseries.gripapt.trlrej	= {'rt','mt','mga','rtmga','tga'};
cfg.tseries.gripapt.movname = movname;
%cfg.tseries.gripapt.basecorr.win    = @(x) [x(2)+0.1 x(2)+0.3];
%cfg.tseries.gripapt.basecorr.val    = 0.06;
%cfg.tseries.gripapt.basecorr.descrip	= 'nanmedian';
%cfg.tseries.gripapt.basecorr.corrflg	= 'sess';
cfg.tseries.gripaptvel.marker	= {'gripti','gripMCPt','gripMCPi'};
cfg.tseries.gripaptvel.axis     = 'yz';
cfg.tseries.gripaptvel.trlrej	= {'rt','mt','mga','rtmga','tga','tgv'};
cfg.tseries.gripaptvel.movname	= movname;
cfg.tseries.gripori.marker	= {'gripti','gripMCPt','gripMCPi'};
cfg.tseries.gripori.axis	= 'yz';
cfg.tseries.gripori.trlrej	= {'rt','mt','mga','rtmga','tga'};
cfg.tseries.gripori.movname = movname;

% process timeseries before plotting
%----------------------------------------
cfg.tseries.proc.task               = {'grip'};
cfg.tseries.proc.grip.marker        = {{'t','i'},{'meanti','obj'},{'MCP','t'},{'MCP','i'}};
cfg.tseries.proc.grip.label     	= {'gripti','d2t','gripMCPt','gripMCPi'};
cfg.tseries.proc.grip.apt.axis    	= {'xyz','yz'};
cfg.tseries.proc.grip.ori.axis    	= {'yz','z'};
cfg.tseries.proc.grip.ori.dir   	= 'cw';
cfg.tseries.proc.grip.ori.offset	= 90;

% trial rejection
%----------------------------------------
cfg.plot.trlrej.listwise = {'apriori','rt'};
cfg.plot.trlrej.pairwise = cfg.plot.param;
cfg.plot.trlrej.dataset  = 'MASKpchip_FILTlp15_MOV_ANA_IDXR';
cfg.plot.trlrej.tt     	 = {'tt','rtt'};
cfg.plot.trlrej.rtt      = {'tt','rtt'};
cfg.plot.trlrej.ttt      = {'ttt','rttt'};
cfg.plot.trlrej.rttt     = {'ttt','rttt'};
cfg.plot.trlrej.itt      = {'itt','ritt'};
cfg.plot.trlrej.ritt     = {'itt','ritt'};
cfg.plot.trlrej.itdiff   = {'itdiff','ritdiff'};
cfg.plot.trlrej.ritdiff  = {'itdiff','ritdiff'};
cfg.plot.trlrej.gt     	 = {'gt','rgt'};
cfg.plot.trlrej.rgt      = {'gt','rgt'};
cfg.plot.trlrej.pv     	 = {'pv','tpv','rtpv'};
cfg.plot.trlrej.tpv      = {'pv','tpv','rtpv'};
cfg.plot.trlrej.rtpv     = {'pv','tpv','rtpv'};
cfg.plot.trlrej.pa     	 = {'pa','tpa','rtpa'};
cfg.plot.trlrej.tpa      = {'pa','tpa','rtpa'};
cfg.plot.trlrej.rtpa     = {'pa','tpa','rtpa'};
cfg.plot.trlrej.pd     	 = {'pd','tpd','rtpd'};
cfg.plot.trlrej.tpd      = {'pd','tpd','rtpd'};
cfg.plot.trlrej.rtpd     = {'pd','tpd','rtpd'};
cfg.plot.trlrej.mga      = {'mga','tmga','rtmga'};
cfg.plot.trlrej.tmga     = {'mga','tmga','rtmga'};
cfg.plot.trlrej.rtmga    = {'mga','tmga','rtmga'};
cfg.plot.trlrej.tga      = {'tga','dtga'};
cfg.plot.trlrej.dtga     = {'tga','dtga'};
cfg.plot.trlrej.ega      = {'ega','dega'};
cfg.plot.trlrej.dega     = {'ega','dega'};
cfg.plot.trlrej.tgo      = {'tgo','dtgo','adtgo'};
cfg.plot.trlrej.dtgo     = {'tgo','dtgo','adtgo'};
cfg.plot.trlrej.adtgo    = {'tgo','dtgo','adtgo'};
cfg.plot.trlrej.ego      = {'ego','dego','adego'};
cfg.plot.trlrej.dego     = {'ego','dego','adego'};
cfg.plot.trlrej.adego    = {'ego','dego','adego'};
%cfg.plot.reptcontrib	= 10;   % main effects
%cfg.plot.reptcontrib	= 6;    % interaction
cfg.plot.reptcontrib	= 4;    % test
cfg.plot.subjcontrib	= 1;


% report and statistics
%----------------------------------------
%cfg.paramreport.write   = {'groupuni'};
%cfg.paramreport.vars    = {'all','-trlbeg','-trlend','-trloff','-RT','-MT','-GT','-HRT','-early','-late','-still','-pause','-t_*'};
%cfg.paramreport.reject  = {'all'};
%cfg.paramreport.exclude = 'cases';
cfg.paramreport         = 'no';
cfg.stats               = 'no';


% selecting plot types
%----------------------------------------
if ~isfalse(cfg.plot.tseries),  cfg.tseriesplot = 'yes';    end
if ~isempty(regexp(cfg.plot.fact,'slant','once'))
    cfg.paramplot.type      = 'errorbar';
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


% plotting configuration
%----------------------------------------
% Colors
%cfg.plot.flcl	= str2rgb('rgh');
%cfg.plot.flcl	= str2rgb('boh');
%cfg.plot.flcl	= str2rgb('hhhhhh');
cfg.plot.flcl	= str2rgb('hmchmc');
%cfg.plot.flcl	= str2rgb('bobobo');
%cfg.plot.flcl	= str2rgb('rgrgrg');
cfg.plot.lncl	= str2rgb('kkkkkk');
%cfg.plot.lncl	= str2rgb('bboo');
%cfg.plot.lncl	= str2rgb('rgbcmyo');
%cfg.plot.lnstyle	= {'-','- ','-'};


% execute plotting wrapper
%----------------------------------------
km_plot(cfg);
