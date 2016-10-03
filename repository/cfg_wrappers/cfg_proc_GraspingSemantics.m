function cfg_proc_GraspingSemantics
%--------------------------------------------------------------------------
% CFG_PROC_GRASPINGSEMANTICS is a wrapper function for kinematic data
% processing and analyses for the Grasping Semantics experiment.
%
% Evaluate the code below to add the correct paths for this experiment:
% startup('km');
% cdata('gs');
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
dir_root	= fullfile(filesep,'home','action','lenver','Data','GraspingSemantics');
subjects	= {'subj01','subj02','subj03','subj04','subj05',...
               'subj06','subj07','subj08','subj09','subj10',...
               'subj11','subj12','subj13','subj14','subj15',...
               'subj16','subj17','subj18','subj19','subj20',...
               'pilot99'};

%subjsel = 1:4;
%subjsel = 5:8;
%subjsel = 9:12;
%subjsel = 13:15;   % subj 16 has no trigger events
subjsel = 5;

settask = {'movplot'};
% settask = {'read','movdef','movparam','movplot','trlrej','CI'};

%% processing
%----------------------------------------
% loop over task sets
if ~iscell(settask),    settask = {settask};    end
for st = 1:length(settask)
currsettask = settask{st};

% basic configuration
cfg             = [];
cfg.exp        	= 'GraspingSemantics';
cfg.save        = 'no';
       
switch lower(currsettask)
    
    %% read
    case 'read'
        cfg.dataset	= '';
        cfg.task	= {'read','chunk','hemi','axes','calc','interp','readlog','alignlog','artfctdef','maskart','filter','grip','deriv'};
        
        % read in data
        cfg.read.type           = 'polhemus';
        cfg.read.expfun        	= cfg.exp;
        % cut data in chunks at breaks
        cfg.cutbreaks.maxbreak	= '6*samp';     % maximum of 6 missing samples
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
        cfg.alignlog.logvar     = 't_ITI';
        % artifact detection
        cfg.artfctdef.type  	= 'outlier';
        % artifact masking
        cfg.maskart           	= 'pchip';
        % filter
        cfg.filter.param        = 'pos';
        cfg.filter.mode         = 'lowpass';
        cfg.filter.type         = 'but';
        cfg.filter.freq         = 15;
        cfg.filter.ord          = 6;
        % grip
        cfg.grip.marker         = {{'t','i'},{'t','j'},{'i','j'}};
        cfg.grip.label          = {'ti','tj','ij'};
        cfg.grip.apt.axis       = {'xyz','yz'};
        cfg.grip.ori.axis       = {'yz','z'};
        cfg.grip.ori.dir        = 'cw';
        cfg.grip.ori.offset     = 90;
        % derivatives: velocity and acceleration
        cfg.deriv.method        = 'gradient';
        cfg.deriv.axis         	= {'xyz','yz','x'};
        
    %% movdef
    case 'movdef'
        cfg.dataset	= 'MASKpchip_FILTlp15';
        %cfg.task    = {'movdef','movcheck','trldef','cuttrials','movplot'};
        %cfg.task    = {'trldef','cuttrials','movdef','movcheck'};
        cfg.task    = {'movdef','trldef','cuttrials','movcheck'};
                
        % define movement on- and off-sets
        cfg.movdef.msi                  = 'yes';
        % movement onset
        cfg.movdef(1).twin              = 1;
        cfg.movdef(1).param             = {'pos','vel','vel','gripvel','logfile'};
        cfg.movdef(1).pos.marker        = {'t','i'};
        cfg.movdef(1).pos.crit          = '(x>0.15 & x<0.25) & (y<-0.3) & (z<0.2)';
        cfg.movdef(1).pos.fun           = 'mov';
        cfg.movdef(1).pos.twin          = [-1 0];
        cfg.movdef(1).pos.comb          = 'all';
        cfg.movdef(1).vel{1}.marker     = {'t','i'};
        cfg.movdef(1).vel{1}.axis       = 'yz';
        cfg.movdef(1).vel{1}.crit       = 0.02;
        cfg.movdef(1).vel{1}.slide      = 0.005;
        cfg.movdef(1).vel{1}.conseqflg  = 'merge';
        cfg.movdef(1).vel{1}.fun        = 'mov';
        cfg.movdef(1).vel{1}.comb       = 'any';
        cfg.movdef(1).vel{2}.marker     = {'t','i'};
        cfg.movdef(1).vel{2}.axis       = 'yz';
        cfg.movdef(1).vel{2}.crit       = -Inf;
        cfg.movdef(1).vel{2}.slide      = 'no';
        cfg.movdef(1).vel{2}.fun        = '(1 - (dat./(maxdat/2))) .* mov';
        cfg.movdef(1).vel{2}.comb       = 'mean';
        cfg.movdef(1).gripvel.marker    = 'ti';
        cfg.movdef(1).gripvel.axis      = 'yz';
        cfg.movdef(1).gripvel.abs       = 'yes';
        cfg.movdef(1).gripvel.crit      = 0.01;
        cfg.movdef(1).gripvel.slide     = 'yes';
        cfg.movdef(1).gripvel.conseqflg = 'merge';
        cfg.movdef(1).gripvel.fun       = '(1 - (dat./maxdat))  .* mov';
        cfg.movdef(1).gripvel.comb      = 'mean';
        cfg.movdef(1).logfile.var       = '^t_grasp';
        cfg.movdef(1).logfile.twin      = 5;
        cfg.movdef(1).logfile.conv      = 'betapdfmaxcentered';
        cfg.movdef(1).logfile.betapdf.a = 7;
        cfg.movdef(1).logfile.betapdf.b = 3;
        cfg.movdef(1).logfile.fun       = 'mov/maxmov';
        % movement offset
        cfg.movdef(2).twin              = 1;
        cfg.movdef(2).param             = {'pos','pos','pos','pos','vel','vel','gripapt','gripvel','gripvel','probwin'};
        cfg.movdef(2).pos{1}.marker     = {'t','i'};
        cfg.movdef(2).pos{1}.crit       = '(x>0.12 & x<0.25) & (y>-0.25 & y<0) & (z>0.03 & z<0.12)';
        cfg.movdef(2).pos{1}.fun        = 'mov';
        cfg.movdef(2).pos{1}.comb       = 'all';
        cfg.movdef(2).pos{2}.marker     = {'t','i'};
        cfg.movdef(2).pos{2}.crit       = 'x<0.25';
        cfg.movdef(2).pos{2}.twin       = [-0.2 0];
        cfg.movdef(2).pos{2}.fun        = 'mov';
        cfg.movdef(2).pos{2}.comb       = 'any';
        cfg.movdef(2).pos{3}.marker     = {'t','i'};
        cfg.movdef(2).pos{3}.axis       = 'x';
        cfg.movdef(2).pos{3}.crit       = -Inf;
        cfg.movdef(2).pos{3}.fun        = '1 - (dat-0.2)./0.05';
        cfg.movdef(2).pos{3}.comb       = 'mean';
        cfg.movdef(2).pos{4}.marker     = {'t','i'};
        cfg.movdef(2).pos{4}.axis       = 'z';
        cfg.movdef(2).pos{4}.crit       = -Inf;
        cfg.movdef(2).pos{4}.fun        = '1 - dat.*2';
        cfg.movdef(2).pos{4}.comb       = 'mean';
        %cfg.movdef(2).pos{3}.fun        = '1 - (dat-meddat)./0.05';
        cfg.movdef(2).pos{3}.comb       = 'mean';
        cfg.movdef(2).vel{1}.marker     = {'t','i','j','w'};
        cfg.movdef(2).vel{1}.axis       = 'yz';
        cfg.movdef(2).vel{1}.crit       = -Inf;
        cfg.movdef(2).vel{1}.slide      = 'no';
        cfg.movdef(2).vel{1}.fun        = '(1 - (dat./(maxdat/2))) .* mov';
        cfg.movdef(2).vel{1}.comb       = 'mean';
        cfg.movdef(2).vel{2}.marker     = {'t','i'};
        cfg.movdef(2).vel{2}.axis       = 'x';
        cfg.movdef(2).vel{2}.crit       = -Inf;
        cfg.movdef(2).vel{2}.slide      = 'no';
        cfg.movdef(2).vel{2}.fun        = '(1 - (dat./(maxdat/2))) .* mov';
        cfg.movdef(2).vel{2}.comb       = 'mean';
        cfg.movdef(2).gripapt.marker	= 'ti';
        cfg.movdef(2).gripapt.axis      = 'xyz';
        cfg.movdef(2).gripapt.crit      = 'dat>0.03 & dat<0.11';
        cfg.movdef(2).gripapt.fun       = '(1.25 - (dat./0.1)) .* mov';
        cfg.movdef(2).gripapt.comb      = 'max';
        cfg.movdef(2).gripvel{1}.marker = 'ti';
        cfg.movdef(2).gripvel{1}.axis   = 'xyz';
        cfg.movdef(2).gripvel{1}.abs    = 'yes';
        cfg.movdef(2).gripvel{1}.crit   = 'dat<0.05';
        cfg.movdef(2).gripvel{1}.slide  = 'no';
        cfg.movdef(2).gripvel{1}.fun    = 'mov';
        cfg.movdef(2).gripvel{1}.twin   = [0 0.2];
        cfg.movdef(2).gripvel{1}.comb   = 'any';
        cfg.movdef(2).gripvel{2}.marker = 'ti';
        cfg.movdef(2).gripvel{2}.axis   = 'xyz';
        cfg.movdef(2).gripvel{2}.abs    = 'no';
        cfg.movdef(2).gripvel{2}.crit   = 'dat>-0.1 & dat<0';
        cfg.movdef(2).gripvel{2}.slide  = 'no';
        cfg.movdef(2).gripvel{2}.fun    = '(1 - (dat./(mindat/2))) .* mov';
        cfg.movdef(2).gripvel{2}.comb   = 'max';
        cfg.movdef(2).probwin.winidx    = [1 2];
        cfg.movdef(2).probwin.fun       = 'betapdf(x,2,2)/10+0.85';
        %cfg.movdef(2).probwin.fun      = 'betapdf((x/2)+0.03,2,5)/10+0.75';
        %cfg.movdef(2).probwin.fun      = 'betapdf((x+1)/3,2,3)/1.5';
        %cfg.movdef(2).probwin.fun      = '1-0.2*x';
        % define trials
        cfg.trldef.param            = 'log';
        cfg.trldef.beg              = 't_stim';
        cfg.trldef.end              = 't_rtrn';
        cfg.trldef.prepad           = -0.5;
        cfg.trldef.postpad          = 0.5;
        % devide data over trials
        cfg.cuttrials.align         = 'no';
        cfg.cuttrials.timevar       = '^t_';
        % check movements and keep only that pass
        cfg.movcheck.param          = {'dur','peakvel','ega','bpos','epos'};
        cfg.movcheck.dur.on         = [0.4 3];
        cfg.movcheck.dur.off        = 0;
        cfg.movcheck.peakvel.axis	= 'xyz';
        cfg.movcheck.peakvel.marker	= 'w';
        cfg.movcheck.peakvel.crit   = 0.3;
        cfg.movcheck.ega.axis       = 'yz';
        cfg.movcheck.ega.marker     = 'ti';
        cfg.movcheck.ega.crit       = [0.03 0.11];
        cfg.movcheck.bpos.marker    = 'meanti';
        cfg.movcheck.bpos.crit      = '(x>0.15) & (x<0.25) & (y<-0.3) & (z<0.15)';
        cfg.movcheck.bpos.type      = 'incl';
        cfg.movcheck.epos.marker    = 'meanti';
        cfg.movcheck.epos.crit      = '(x>0.12) & (x<0.25) & (y>-0.25) & (y<0) & (z>0.03) & (z<0.15)';
        cfg.movcheck.epos.type      = 'incl';
        % visually inspect the movements
        cfg.movplot.param           = {'vel','gripapt','gripaptvel','pos'};
        cfg.movplot.vel.marker      = 'meanti';
        cfg.movplot.vel.axis        = 'xyz';
        cfg.movplot.gripapt.marker  = 'ti';
        cfg.movplot.gripapt.axis    = 'xyz';
        cfg.movplot.pos.marker      = 'meanti';
        cfg.movplot.pos.axis        = 'x';
        cfg.movplot.xlim            = [-0.1 4];
        
    %% movparam
    case 'movparam'
        cfg.dataset	= 'MASKpchip_FILTlp15_MOV';
        cfg.task    = {'movparam'};
        
        % obtain movement parameters
        cfg.movparam.param          = {'rt','mt','td','mv','pv','tpv','rtpv','pa','tpa','rtpa','pd','tpd','rtpd','mga','tmga','rtmga','ega','dega'};	% 'epos'
        cfg.movparam.pos.marker     = 'meanti';
        cfg.movparam.pos.axis       = {'x','y','z'};
        cfg.movparam.dist.marker    = 'meanti';
        cfg.movparam.dist.axis      = 'yz';
        cfg.movparam.vel.marker     = 'meanti';
        cfg.movparam.vel.axis       = 'xyz';
        cfg.movparam.tvel.marker    = 'j';
        cfg.movparam.tvel.axis      = 'yz';
        cfg.movparam.tvel.crit      = 0.3;
        cfg.movparam.tvel.slide     = 0.03;
        cfg.movparam.acc.marker     = 'meanti';
        cfg.movparam.acc.axis       = 'xyz';
        cfg.movparam.gripapt.marker	= 'ti';
        cfg.movparam.gripapt.axis	= 'xyz';
        cfg.movparam.refapt         = 'gripapt';
        cfg.movparam.endoffset      = -0.05;
        cfg.movparam.refoffset      = 0.2;
    
    %% movplot
    case 'movplot'
        cfg.dataset	= 'MASKpchip_FILTlp15_MOV_ANA';
        %cfg.dataset	= 'MASKpchip_FILTlp15_MOV';
        cfg.task    = {'movplot'};
        
        % visually inspect the movements
        cfg.movplot.param           = {'vel','gripapt','gripaptvel','pos'};
        cfg.movplot.vel.marker      = 'meanti';
        cfg.movplot.vel.axis        = 'xyz';
        cfg.movplot.gripapt.marker  = 'ti';
        cfg.movplot.gripapt.axis    = 'xyz';
        cfg.movplot.pos.marker      = 'meanti';
        cfg.movplot.pos.axis        = 'x';
        cfg.movplot.xlim            = [-0.1 4];
        
    %% trlrej
    case 'trlrej'
        cfg.dataset	= 'MASKpchip_FILTlp15_MOV_ANA';
        cfg.task    = {'trlrej'};
        
        % reject trials
        cfg.trlrej.fun              = 'auto';
        cfg.trlrej.apriorifun       = cfg.exp;
        cfg.trlrej.mode             = {'apriori','rt','mt','td','mv','pv','mga','ega'};
        cfg.trlrej.rt               = [0.3 1.5];
        cfg.trlrej.mt               = [0.4 2];
        cfg.trlrej.td               = [0.2 0.7];
        cfg.trlrej.ega              = [0.03 0.11];
        
    %% CI
    case 'ci'
        cfg.dataset	= 'MASKpchip_FILTlp15_MOV_ANA';
        cfg.task    = {'movparamCI','eposerr'};
        
        % confidence intervals of the movement parameters
        cfg.movparamCI.param        	= {'rt','mt','td','pv','mga'};
        cfg.movparamCI.trlrej.dataset	= 'MASKpchip_FILTlp15_MOV_ANA_IDXR';
        cfg.movparamCI.fact          	= {'vision','size','color','word','vision x size','size x word','color x word'};
        cfg.movparamCI.lvl.vision     	= [1 2];
        cfg.movparamCI.lvl.size     	= [1 2];
        cfg.movparamCI.lvl.color     	= [1 2];
        cfg.movparamCI.lvl.word     	= [1 2 3 4];
        % calculate end point error
        cfg.eposerr.trlrej.dataset  = 'MASKpchip_FILTlp15_MOV_ANA_IDXR';
        cfg.eposerr.marker          = {'t','i','j','meanti'};
        cfg.eposerr.axis            = 'xyz';
        cfg.eposerr.endoffset       = -0.05;
        cfg.eposerr.fact            = {'vision','size','color','word','vision x size','size x word','color x word'};
        cfg.eposerr.lvl.vision     	= [1 2];
        cfg.eposerr.lvl.size     	= [1 2];
        cfg.eposerr.lvl.color     	= [1 2];
        cfg.eposerr.lvl.word     	= [1 2 3 4];
        cfg.eposerr.expfun          = 'no';
        
end


%% execute task set
%----------------------------------------
% loop over subjects and sessions
for s = subjsel
    cfg.subj        = subjects{s};
    cfg.dir.data    = fullfile(dir_root,cfg.subj,'raw');
    cfg.dir.log     = fullfile(dir_root,cfg.subj,'log');
    cfg.dir.preproc = fullfile(dir_root,cfg.subj,'preproc');
    cfg.dir.ana     = fullfile(dir_root,cfg.subj,'ana');
    
    if ~isfalse(cfg.task)
        km_dotask(cfg);
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