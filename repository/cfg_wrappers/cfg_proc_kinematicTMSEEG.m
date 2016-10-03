function cfg_proc_kinematicTMSEEG
%--------------------------------------------------------------------------
% CFG_PROC_KINEMATICTMSEEG is a wrapper function for EEG-TMS kinematic data
% processing and analyses for the Grasping Bino-Mono EEG-TMS experiment.
%
% Evaluate the code below to add the correct paths for this experiment:
% startup('Kinemagic');
% cdata('tmseeg');
%
% credits: the code is original and my own, but the code organization is
% based on that of FieldTrip: many thanks to the FieldTrip developers!
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

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

% subj07 has no object marker
%subjsel = 8:14;
%subjsel = 15:21;
%subjsel = 22:28;
%subjsel = 29:35;
subjsel = 8:35;
sesssel	= 1:3;

settask = {'ana','trlrej'};
%settask = {'read','ana','trlrej','CI'};

%% processing
%----------------------------------------
% loop over task sets
if ~iscell(settask),    settask = {settask};    end
for st = 1:length(settask)
currsettask = settask{st};

% basic configuration
cfg             = [];
cfg.exp        	= 'kinematicTMSEEG';
cfg.save        = 'yes';
       
switch lower(currsettask)
    
    %% read
    case 'read'
        cfg.dataset	= '';
        cfg.task	= {'read','chunk','hemi','axes','calc','interp','readlog','alignlog','artfctdef','maskart','filter','grip','deriv'};
        
        % read in data
        cfg.read.type           = 'polhemus';  	% 'polhemus', 'minibird'
        cfg.read.expfun        	= cfg.exp;
        % cut data in chunks at breaks
        cfg.cutbreaks.maxbreak	= '3*samp';     % maximum of 3 missing samples
        cfg.cutbreaks.minchunk	= 4;            % minimum of 4 second chunks after cutting
        % hemifield crossing correction
        cfg.hemifieldcross      = 'loose';      % 'no','loose','strict'
        % axes transformation
        cfg.axes.order          = 'xyz';        % 'no','xyz','xzy', etc
        cfg.axes.mirror         = 'no';         % 'no','x','y','z';
        cfg.axes.rotate         = 'no';         % 'no' or rotation in degrees over [x y z] axes
        % calculate extra markers
        cfg.calcmarker.type     = 'mean';
        cfg.calcmarker.marker	= {'t','i'};
        cfg.calcmarker.name     = 'meanti';
        % interpolate
        cfg.interp.method       = 'pchip';
        cfg.interp.freq         = 'max';
        % read in the logfile
        cfg.readlog             = 'yes';
        % align the logfile with the data
        %cfg.alignlog            = 'no';
        cfg.alignlog.evtcode    = 1;
        cfg.alignlog.logvar     = 't_stim';        
        % artifact detection
        cfg.artfctdef.type          = {'log','tms'};
        cfg.artfctdef.log.art       = 't_TMS';
        cfg.artfctdef.log.prepad	= -0.005;
        cfg.artfctdef.log.postpad	= 0.020;
        % artifact masking
        cfg.maskart           	= 'pchip';
        % filter
        cfg.filter.mode         = 'lowpass';
        cfg.filter.type         = 'but';
        cfg.filter.freq         = 15;
        cfg.filter.ord          = 6;
        % grip
        cfg.grip.marker         = {{'t','i'},{'meanti','obj'}};
        cfg.grip.label          = {'gripti','d2t'};
        cfg.grip.apt.axis       = {'xyz','yz'};
        cfg.grip.ori.axis       = {'yz','z'};
        cfg.grip.ori.dir        = 'cw';
        cfg.grip.ori.offset     = 90;
        % derivatives: velocity and acceleration
        cfg.deriv.method        = 'gradient';
        cfg.deriv.axis         	= {'xyz','yz','x'};
        
    %% movdef
    case 'ana'
        cfg.dataset	= 'MASKpchip_FILTlp15';
        %cfg.task    = {'movdef','movcheck','movplot','trldef','cuttrials','movparam'};
        cfg.task    = {'movdef','movcheck','trldef','cuttrials','trlcheck','movparam'};
        
        % define movement on- and off-sets
        %cfg.movdef.param            = 'vel';
        %cfg.movdef.vel.axis         = 'yz';
        %cfg.movdef.vel.marker       = 'meanti';
        cfg.movdef.msi                  = 'yes';
        % movement onset
        cfg.movdef(1).twin              = 1;
        cfg.movdef(1).param             = {'pos','vel','vel','gripapt','gripvel','logfile'};
        cfg.movdef(1).pos.marker        = {'t','i'};
        cfg.movdef(1).pos.refmarker     = 'obj';
        cfg.movdef(1).pos.crit          = '(x>-0.1) & (x<0.1) & (y<-0.15) & (z<-0.2)';
        cfg.movdef(1).pos.fun           = 'mov';
        cfg.movdef(1).pos.twin          = [-1 0];
        cfg.movdef(1).pos.comb          = 'all';
        cfg.movdef(1).vel{1}.marker     = {'t','i'};
        cfg.movdef(1).vel{1}.axis       = 'yz';
        cfg.movdef(1).vel{1}.crit       = 0.02;
        cfg.movdef(1).vel{1}.slide      = 'yes';
        cfg.movdef(1).vel{1}.conseqflg  = 'merge';
        cfg.movdef(1).vel{1}.fun        = 'mov';
        cfg.movdef(1).vel{1}.comb       = 'any';
        cfg.movdef(1).vel{2}.marker     = {'t','i'};
        cfg.movdef(1).vel{2}.axis       = 'yz';
        cfg.movdef(1).vel{2}.crit       = -Inf;
        cfg.movdef(1).vel{2}.slide      = 'no';
        cfg.movdef(1).vel{2}.fun        = '(1 - (dat./maxdat)) .* mov';
        cfg.movdef(1).vel{2}.comb       = 'mean';
        cfg.movdef(1).gripapt.marker	= 'd2t';
        cfg.movdef(1).gripapt.axis      = 'yz';
        cfg.movdef(1).gripapt.crit      = 'dat>0.3';
        cfg.movdef(1).gripapt.fun       = '(dat./maxabsdat) .* mov';
        cfg.movdef(1).gripapt.comb      = 'mean';
        cfg.movdef(1).gripvel.marker    = 'gripti';
        cfg.movdef(1).gripvel.axis      = 'yz';
        cfg.movdef(1).gripvel.abs       = 'yes';
        cfg.movdef(1).gripvel.crit      = 0.01;
        cfg.movdef(1).gripvel.slide     = 'yes';
        cfg.movdef(1).gripvel.conseqflg = 'merge';
        cfg.movdef(1).gripvel.fun       = '(1 - (dat./(2*maxdat)))  .* mov';
        cfg.movdef(1).gripvel.comb      = 'mean';
        %cfg.movdef(1).gripacc.marker   = 'dist2target';
        %cfg.movdef(1).gripacc.axis     = 'yz';
        %cfg.movdef(1).gripacc.crit     = 'dat<=0'; % yes... zero or lower
        %cfg.movdef(1).gripacc.fun      = 'mov';
        cfg.movdef(1).logfile.var       = '^t_grasp';
        cfg.movdef(1).logfile.twin      = 4;
        cfg.movdef(1).logfile.conv      = 'betapdfmaxcentered';
        cfg.movdef(1).logfile.betapdf.a = 7;
        cfg.movdef(1).logfile.betapdf.b = 3;
        cfg.movdef(1).logfile.fun       = 'mov/maxmov';
        % movement offset
        cfg.movdef(2).twin              = 1;
        cfg.movdef(2).param             = {'pos','pos','vel','vel','vel','gripapt','gripapt','gripvel','gripvel','gripvel','probwin'};
        %cfg.movdef(2).param             = {'pos','pos','vel','vel','vel','gripapt','gripapt','gripvel','gripvel','gripvel'};
        cfg.movdef(2).pos{1}.marker     = {'t','i'};
        cfg.movdef(2).pos{1}.refmarker  = 'obj';
        cfg.movdef(2).pos{1}.crit       = '(x>-0.07) & (x<0.07) & (y>-0.07) & (y<0.07) & (z>-0.07) & (z<0.07)';
        cfg.movdef(2).pos{1}.twin       = [0 0.2];
        cfg.movdef(2).pos{1}.fun        = 'mov';
        cfg.movdef(2).pos{1}.comb       = 'all';
        cfg.movdef(2).pos{2}.marker     = {'obj'};
        cfg.movdef(2).pos{2}.axis       = 'x';
        cfg.movdef(2).pos{2}.crit       = -Inf;
        cfg.movdef(2).pos{2}.fun        = '1 - (dat-meddat)./0.02';
        cfg.movdef(2).pos{2}.comb       = 'mean';
        cfg.movdef(2).vel{1}.marker     = {'t','i'};
        cfg.movdef(2).vel{1}.axis       = 'yz';
        cfg.movdef(2).vel{1}.crit       = -Inf;
        cfg.movdef(2).vel{1}.slide      = 'no';
        cfg.movdef(2).vel{1}.fun        = '(1 - (dat./maxdat)) .* mov';
        cfg.movdef(2).vel{1}.comb       = 'mean';
        cfg.movdef(2).vel{2}.marker     = {'t','i','obj'};
        cfg.movdef(2).vel{2}.axis       = 'x';
        cfg.movdef(2).vel{2}.crit       = -Inf;
        cfg.movdef(2).vel{2}.slide      = 'no';
        cfg.movdef(2).vel{2}.fun        = '(1 - (dat./(2*maxdat))) .* mov';
        cfg.movdef(2).vel{2}.comb       = 'mean';
        cfg.movdef(2).vel{3}.marker     = {'obj'};
        cfg.movdef(2).vel{3}.axis       = 'x';
        cfg.movdef(2).vel{3}.crit       = 'dat<0.04';
        cfg.movdef(2).vel{3}.slide      = 'no';
        cfg.movdef(2).vel{3}.twin       = [-1 0];
        cfg.movdef(2).vel{3}.fun        = '(1 - (dat./maxdat)) .* mov';
        cfg.movdef(2).vel{3}.comb       = 'mean';
        cfg.movdef(2).gripapt{1}.marker	= 'd2t';
        cfg.movdef(2).gripapt{1}.axis   = 'yz';
        cfg.movdef(2).gripapt{1}.crit   = 'dat<0.04';
        cfg.movdef(2).gripapt{1}.fun    = '(1 - (dat./(maxdat/2))) .* mov';
        cfg.movdef(2).gripapt{1}.comb   = 'max';
        cfg.movdef(2).gripapt{2}.marker	= 'gripti';
        cfg.movdef(2).gripapt{2}.axis   = 'yz';
        cfg.movdef(2).gripapt{2}.crit   = 'dat<0.11 & dat>0.075';
        cfg.movdef(2).gripapt{2}.fun    = '(3 - (dat./(0.075/2))) .* mov';
        cfg.movdef(2).gripapt{2}.comb   = 'max';
        cfg.movdef(2).gripvel{1}.marker = 'gripti';
        cfg.movdef(2).gripvel{1}.axis   = 'yz';
        cfg.movdef(2).gripvel{1}.abs    = 'yes';
        cfg.movdef(2).gripvel{1}.crit   = 'dat<0.05';
        cfg.movdef(2).gripvel{1}.slide  = 'no';
        cfg.movdef(2).gripvel{1}.fun    = 'mov';
        cfg.movdef(2).gripvel{1}.twin   = [0 0.2];
        cfg.movdef(2).gripvel{1}.comb   = 'any';
        cfg.movdef(2).gripvel{2}.marker = 'gripti';
        cfg.movdef(2).gripvel{2}.axis   = 'yz';
        cfg.movdef(2).gripvel{2}.abs    = 'no';
        cfg.movdef(2).gripvel{2}.crit   = 'dat>-0.05 & dat<0';
        cfg.movdef(2).gripvel{2}.slide  = 'no';
        cfg.movdef(2).gripvel{2}.fun    = '(1 - (dat./mindat)) .* mov';
        cfg.movdef(2).gripvel{2}.comb   = 'max';
        cfg.movdef(2).gripvel{3}.marker = 'gripti';
        cfg.movdef(2).gripvel{3}.axis   = 'yz';
        cfg.movdef(2).gripvel{3}.abs    = 'yes';
        cfg.movdef(2).gripvel{3}.crit   = 0.005;
        cfg.movdef(2).gripvel{3}.slide  = 'yes';
        cfg.movdef(2).gripvel{3}.conseqflg = 'merge';
        cfg.movdef(2).gripvel{3}.fun    = 'mov';
        cfg.movdef(2).gripvel{3}.comb   = 'any';
        %cfg.movdef(2).gripacc.marker   = 'gripti';
        %cfg.movdef(2).gripacc.axis     = 'yz';
        %cfg.movdef(2).gripacc.crit     = 'dat>=0'; % yes... zero or higher
        %cfg.movdef(2).gripacc.fun      = 'mov';
        %cfg.movdef(2).gripacc.comb     = 'max';
        cfg.movdef(2).probwin.winidx    = [1 5];
        cfg.movdef(2).probwin.fun       = 'betapdf((x/2)+0.03,2,5)/10+0.75';
        %cfg.movdef(2).probwin.fun      = 'betapdf((x+1)/3,2,3)/1.5';
        %cfg.movdef(2).probwin.fun      = '1-0.2*x';
        % check movements and keep only that pass
        cfg.movcheck.param          = {'dur','peakvel','ega','pos','bpos'};
        cfg.movcheck.dur.on         = [0.4 3];
        cfg.movcheck.dur.off        = 0.05;
        cfg.movcheck.peakvel.axis	= 'yz';
        cfg.movcheck.peakvel.marker	= 'meanti';
        cfg.movcheck.peakvel.crit   = 0.5;
        cfg.movcheck.ega.axis       = 'yz';
        cfg.movcheck.ega.marker     = 'd2t';
        cfg.movcheck.ega.crit       = [-Inf 0.05];
        cfg.movcheck.pos.marker     = 'i';
        cfg.movcheck.pos.refmarker  = 'obj';
        cfg.movcheck.pos.crit       = '(z>-4*y) & (z<-0.5*y)';
        cfg.movcheck.pos.type       = 'excl';
        cfg.movcheck.bpos.marker    = 'meanti';
        cfg.movcheck.bpos.refmarker = 'obj';
        cfg.movcheck.bpos.crit      = '(y<-0.15) & (z<-0.2)';
        cfg.movcheck.bpos.type      = 'incl';
        % visually inspect the movements
        cfg.movplot.param           = {'vel','gripapt','gripaptvel','pos'};
        cfg.movplot.pos.marker      = 'meanti';
        cfg.movplot.pos.axis        = 'z';
        cfg.movplot.vel.marker      = 'meanti';
        cfg.movplot.vel.axis        = 'yz';
        cfg.movplot.gripapt.marker  = 'gripti';
        cfg.movplot.gripapt.axis    = 'yz';
        cfg.movplot.xlim            = [-0.5 3];        
        % define trials
        cfg.trldef.param            = 'log';
        cfg.trldef.beg              = 't_stim';
        cfg.trldef.end              = 't_rtrn';
        cfg.trldef.prepad           = -0.5;
        cfg.trldef.postpad          = 0.5;
        % devide data over trials
        cfg.cuttrials.align         = 'no';
        cfg.cuttrials.timevar       = '^t_';
        % check trials
        cfg.trlcheck.param          = {'bvel','bpos'};
        cfg.trlcheck.bvel.axis      = 'xyz';
        cfg.trlcheck.bvel.marker	= 'meanti';
        cfg.trlcheck.bvel.crit      = [-0.05 0.05];
        cfg.trlcheck.bpos.marker    = 'meanti';
        cfg.trlcheck.bpos.refmarker = 'obj';
        cfg.trlcheck.bpos.crit      = '(y<-0.15) & (z<-0.2)';
        cfg.trlcheck.bpos.type      = 'incl';
        % obtain movement parameters
        cfg.movparam.param          = {'rt','mt','tt','rtt','ttt','rttt','itt','ritt','itdiff','ritdiff','gt','rgt','td','tcd','gcd','mv','mtv','mgv','pv','tpv','rtpv','pa','tpa','rtpa','pd','tpd','rtpd','nvc','tnvc','gnvc','pgv','npgv','tgv','mggv','gcgd','mga','tmga','rtmga','tga','dtga','ega','dega','tgo','dtgo','adtgo','ego','dego','adego'};	% 'epos'
        %cfg.movparam.param          = {'rt','mt','td','mv','mav','pv','mga','ega','dego'};	% 'epos'
        cfg.movparam.pos.marker     = 'meanti';
        cfg.movparam.pos.axis       = {'x','y','z'};
        %cfg.movparam.dist.marker    ={'meanti','t','i','MCP'};
        cfg.movparam.dist.marker    = 'meanti';
        cfg.movparam.dist.axis      = 'yz';
        %cfg.movparam.vel.marker     = {'meanti','t','i','MCP'};
        cfg.movparam.vel.marker     = 'meanti';
        cfg.movparam.vel.axis       = 'yz';
        cfg.movparam.tvel.marker    = 'MCP';
        cfg.movparam.tvel.axis      = 'yz';
        cfg.movparam.tvel.crit      = 0.3;
        cfg.movparam.tvel.slide     = 0.03;
        cfg.movparam.ttvel.marker   = 't';
        cfg.movparam.ttvel.axis     = 'yz';
        cfg.movparam.ttvel.crit     = 0.3;
        cfg.movparam.ttvel.slide    = 0.03;
        cfg.movparam.itvel.marker   = 'i';
        cfg.movparam.itvel.axis     = 'yz';
        cfg.movparam.itvel.crit     = 0.3;
        cfg.movparam.itvel.slide    = 0.03;
        cfg.movparam.acc.marker     = 'meanti';
        cfg.movparam.acc.axis       = 'yz';
        cfg.movparam.gripapt.marker	= 'gripti';
        cfg.movparam.gripapt.axis	= 'yz';
        cfg.movparam.gripori.marker	= 'gripti';
        cfg.movparam.gripori.axis   = 'yz';
        cfg.movparam.gripaptvel.marker	= 'gripti';
        cfg.movparam.gripaptvel.axis	= 'yz';
        cfg.movparam.refori         = 'gripori';
        cfg.movparam.refapt         = 'gripapt';
        cfg.movparam.appoffset      = -0.1;
        cfg.movparam.endoffset      = -0.05;
        cfg.movparam.refoffset      = 0.1;
            
    %% movplot
    case 'movplot'
        cfg.dataset	= 'MASKpchip_FILTlp15_MOV_ANA';
        cfg.task    = {'movplot'};
        
        % visually inspect the movements
        cfg.movplot.param           = {'vel','gripapt','gripaptvel','pos'};
        cfg.movplot.pos.marker      = 'meanti';
        cfg.movplot.pos.axis        = 'z';
        cfg.movplot.vel.marker      = 'meanti';
        cfg.movplot.vel.axis        = 'yz';
        cfg.movplot.gripapt.marker  = 'gripti';
        cfg.movplot.gripapt.axis    = 'yz';
        cfg.movplot.xlim            = [-0.1 4];
        
    %% trlrej
    case 'trlrej'
        cfg.dataset	= 'MASKpchip_FILTlp15_MOV_ANA';
        cfg.task    = {'trlrej'};
        
        % reject trials
        %cfg.trlrej.fun             = 'auto';
        cfg.trlrej.fun              = 'exp';
        cfg.trlrej.apriori.fun      = cfg.exp;
        cfg.trlrej.param            = {'all'};
        cfg.trlrej.apriori.param    = {'rt','mt','td','mv','pv','mga','ega','dego'};
        cfg.trlrej.reportparam      = {'apriori','rt','mt','td','mv','pv','mga','ega','dego'};
        %cfg.trlrej.apriori.param    = {'rt','mt','td_meanti','mv_meanti','pv_meanti','mga','ega','dego'};
        %cfg.trlrej.reportparam      = {'apriori','rt','mt','td_meanti','mv_meanti','pv_meanti','mga','ega','dego'};
        
    %% CI
    case 'ci'
        cfg.dataset	= 'MASKpchip_FILTlp15_MOV_ANA';
        cfg.task    = {'movparamCI','eposerr'};
        
        % confidence intervals of the movement parameters
        cfg.movparamCI.param        	= {'rt','mt','td','gcd','tgv','mga','tga','dtga','tgo','dtgo','ega','dega','ego','dego'};
        cfg.movparamCI.trlrej.dataset	= 'MASKpchip_FILTlp15_MOV_ANA_IDXR';
        cfg.movparamCI.trlrej.listwise 	= {'apriori','rt','mt','td','mga','ega','dego'};
        cfg.movparamCI.trlrej.pairwise 	= cfg.movparamCI.param;
        cfg.movparamCI.fact          	= {'vision','verthori','vision x verthori'};
        cfg.movparamCI.lvl.verthori    	= [1 2];
        cfg.movparamCI.lvl.vision     	= [1 2];
        % calculate end point error
        cfg.eposerr.trlrej.dataset  = 'MASKpchip_FILTlp15_MOV_ANA_IDXR';
        cfg.eposerr.trlrej.listwise	= {'apriori','rt','mt','td','mga','ega','dego'};
        cfg.eposerr.marker          = {'i','t','MCP','meanti','obj'};
        cfg.eposerr.axis            = 'xyz';
        cfg.eposerr.movname         = 'tmovement';
        % cfg.eposerr.endoffset     = -0.05;    % for movname: movement
        cfg.eposerr.endoffset       = 0;        % for movname: tmovement
        cfg.eposerr.fact            = {'slant','vision x slant'};
        cfg.eposerr.lvl.slant       = 0:15:90;
        cfg.eposerr.lvl.vision      = [1 2];
        cfg.eposerr.expfun          = cfg.exp;
        
end


%% execute task set
%----------------------------------------
% loop over subjects and sessions
for s = subjsel
    cfg.subj        = subjects{s};
    cfg.dir.data    = fullfile(dir_root,cfg.subj,'kinematics','raw');
    cfg.dir.log     = fullfile(dir_root,cfg.subj,'log');
    cfg.dir.preproc = fullfile(dir_root,cfg.subj,'kinematics','preproc');
    cfg.dir.ana     = fullfile(dir_root,cfg.subj,'kinematics','ana');
    for ss = sesssel
        cfg.sess = sessions{ss};
        
        if ~isfalse(cfg.task)
            km_dotask(cfg);
        end
        
    end
end

% end settask-loop
end

%----------------------------------------
% UNDEFINED:the options below mostly restructure the data, but only
%           marginally affect the data itself > the dataset is not changed
% read:     read in the data based on cfg.dataset > dataset is raw
% chunk:    cut data in chunks at sampling breaks
% hemi:     correct hemifield crossings, only applicable to data acquired
%           using a magnetic system (e.g. minibird, polhemus).
% axes:     transform axes (swap, mirror, rotate), especially applicable to
%           data acquired using a magnetic system (e.g. minibird, polhemus)
% calc:     add extra calculated markers based on the real markers, this
%           can be used to create for example a mean marker.
% interp:   interpolate data to a constant sampling rate. Most acquisition
%           software is prone to missing data points, but a constant
%           sampling rate is required for filtering algorithms.
% artfctdef:define artifacts, based on events or automatically
% grip:     add grip aperture and orientation
% deriv:    add derivatives (velocity and/or acceleration)
% movdef:	define movement on- and off-sets
% movcheck: check defined movements and keep only those that pass
% readlog:  read-in the logfile
% trldef:   define trials based on the log
% cuttrials:devide the data according to trial definitions
%
%----------------------------------------
% MASKART:  data-points are cut-out and replaced by NaNs or interpolated
%           data > the dataset is 'maskart'
% maskart:  mask artifacts with NaNs or interpolate
%
%----------------------------------------
% FILTER:   your filter settings can have a profound effect on your data >
%           the dataset is 'filter'
% filter:   filter your data: lowpass, bandpass, highpass, using a
%           butterworth, chebyshev or finite impulse response filter
%
%----------------------------------------
% ANA:      analyzing the data has little effect on the data (only the
%           cutting of the data into trials), but does change it's
%           organization
% movparam: obtain movement parameters
% eposerr:  calculate the end point error (multi-dimensional confidence
%           intervals)