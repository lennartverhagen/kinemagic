function CP_analysis
%--------------------------------------------------------------------------
% This is the main function of CP toolbox, 
% all parts of the analysis are started from here.
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------


%% things to choose
%settask	= {'read','ana','movplot','trlrej','eposerr','traject','plotparam','paramreport','plottraject','plotCI','plotepe'};
settask     = {'plotparam'};
subjsel         = [1:12 14];

%% initialize
dir_CP = CP_initialize;

%% basic parameters
% subject and session selection
subjects        = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15'};
sessions        = {''};

% initialize configuration
cfg             = [];
cfg.save        = 'yes';
cfg.exp        	= 'CommPoint';
cfg.dataset     = '';
cfg.subj        = 'placeholder';
cfg.sessions    = sessions;
cfg.dir.root    = dir_CP;
cfg.dir.data    = fullfile(dir_CP,'raw');
cfg.dir.log     = fullfile(dir_CP,'log');
cfg.dir.preproc = fullfile(dir_CP,'preproc');
cfg.dir.ana     = fullfile(dir_CP,'ana');
cfg.dir.report  = fullfile(dir_CP,'group','report');

% default task requirements
flg_exec	= 'cluster';    % 'local', 'cluster'
memreq      = 1024^3;       % 1 GB
timreq      = 20 * 60;      % 20 minutes

%% configure task sets
if ~iscell(settask),    settask = {settask};    end
for st = 1:length(settask)
    currsettask = settask{st};
    
    switch lower(currsettask)
        
        %% read
        case 'read'
            % initialize
            cfg.dataset	= '';
            cfg.task	= {'read','chunk','hemi','axes','calc','interp','readlog','alignlog','filter','grip','deriv','trldef','cuttrials'};
            memreq      = 1024^3;	% 1 GB
            timreq      = 5 * 60;	% 5 minutes
            
            % read in data
            cfg.read.type           = 'polhemus';  	% 'polhemus', 'minibird'
            cfg.read.expfun        	= 'CP_readfun';
            % cut data in chunks at breaks
            cfg.cutbreaks.maxbreak	= '7*samp';     % maximum of 7 missing samples
            cfg.cutbreaks.minchunk	= 4;            % minimum of 4 second chunks after cutting
            % hemifield crossing correction
            cfg.hemifieldcross      = 'loose';      % 'no','loose','strict'
            % axes transformation
            cfg.axes.order          = 'yxz';        % 'no','xyz','xzy', etc
            cfg.axes.mirror         = 'xy';         % 'no','x','y','z';
            cfg.axes.rotate         = 'no';         % 'no' or rotation in degrees over [x y z] axes
            % calculate extra markers
            cfg.calcmarker.type     = 'mean';
            cfg.calcmarker.marker	= {'hi','hp','w'};
            cfg.calcmarker.name     = 'h';
            % interpolate
            cfg.interp.method       = 'pchip';
            cfg.interp.freq         = 'max';
            % read in the logfile
            cfg.readlog.expfun      = 'CP_readlogfun';
            cfg.readlog.prepad      = 'logfile_';
            cfg.readlog.postpad     = '*.txt';
            % align the logfile with the data
            %cfg.alignlog            = 'no';
            cfg.alignlog.evtcode    = 1;
            cfg.alignlog.logvar     = 'tim_trl_instruct';
            cfg.alignlog.initidx    = 1:6;
            cfg.alignlog.err        = 0.1;
            cfg.alignlog.adjust     = 'log';
            cfg.alignlog.timevar    = '^tim_';
            % filter
            cfg.filter.mode         = 'lowpass';
            cfg.filter.type         = 'but';
            cfg.filter.freq         = 15;
            cfg.filter.ord          = 6;
            % grip
            cfg.grip.marker         = {{'hi','i'},{'hi','hp'}};
            cfg.grip.label          = {'index','hand'};
            cfg.grip.apt.axis       = {'xyz','yz'};
            cfg.grip.ori.axis       = {'yz','z'};
            cfg.grip.ori.dir        = 'cw';
            cfg.grip.ori.offset     = 90;
            % derivatives: velocity and acceleration
            cfg.deriv.method        = 'gradient';
            cfg.deriv.axis         	= {'xyz','yz','x'};
            % define trials
            cfg.trldef.catruns      = 'no';
            cfg.trldef.param        = 'log';
            cfg.trldef.beg          = 'tim_trl_instruct';
            cfg.trldef.stim         = 'tim_trl_instruct';
            cfg.trldef.end          = 'tim_trl_end';
            cfg.trldef.prepad     	= -1;
            cfg.trldef.postpad     	= 0.5;
            % devide data over trials
            cfg.cuttrials.align     = 'no';
            cfg.cuttrials.timevar   = '^tim_';
            
            %% movdef
        case 'ana'
            cfg.dataset	= 'FILTlp15';
            cfg.task    = {'movdef','movcheck','trlcheck','movselect'};
            %cfg.task    = {'movdef','movcheck','trlcheck','movplot'};
            memreq      = 500 * 1024^2;	% 500 MB
            timreq      = 5 * 60;       % 5 minutes
            
            % define movement on- and off-sets
            cfg.movdef.param            = 'vel';
            cfg.movdef.vel.marker       = 'i';
            cfg.movdef.vel.axis         = 'yz';
            cfg.movdef.vel.crit         = 0.1;
            cfg.movdef.vel.slide        = 'no';
            cfg.movdef.vel.consecflg    = 'both';
            % visually inspect the movements
            cfg.movplot.param           = {'pos','vel','gripapt','gripaptvel'};
            cfg.movplot.pos.marker      = {'i','hi','hp','w','h'};
            cfg.movplot.pos.axis        = {'y','x','z'};
            cfg.movplot.vel.marker      = {'i','hi','hp','w','h'};
            cfg.movplot.vel.axis        = {'yz','xyz'};
            cfg.movplot.timeref         = 'rel';
            cfg.movplot.xlim            = [-0.5 2.5];
            % check movements and keep only that pass
            cfg.movcheck(1).param      	= {'dur','peakvel','bpos','epos'};
            cfg.movcheck(1).dur.on     	= [0.2 2];
            cfg.movcheck(1).dur.off   	= 0;
            cfg.movcheck(1).peakvel.marker	= 'i';
            cfg.movcheck(1).peakvel.axis	= 'yz';
            cfg.movcheck(1).peakvel.crit	= 0.1;
            cfg.movcheck(1).bpos.marker	= 'i';
            cfg.movcheck(1).bpos.crit	= 'y<-0.5 & z<0.2';
            cfg.movcheck(1).bpos.type	= 'incl';
            cfg.movcheck(1).epos.marker	= 'i';
            cfg.movcheck(1).epos.crit	= 'y<-0.30 & y>-0.7 & x<0.3 & x>-0.1 & z<0.5 & z>0.25';
            cfg.movcheck(1).epos.type	= 'incl';
            cfg.movcheck(2)             = cfg.movcheck(1);
            cfg.movcheck(2).bpos.crit	= cfg.movcheck(1).epos.crit;
            cfg.movcheck(2).epos.crit	= cfg.movcheck(1).bpos.crit;
            % check trials
            cfg.trlcheck.param          = {'bvel','bpos'};
            cfg.trlcheck.bvel.axis      = 'xyz';
            cfg.trlcheck.bvel.marker	= {'i','h'};
            cfg.trlcheck.bvel.crit      = [-0.2 0.4];
            cfg.trlcheck.bpos.marker    = 'i';
            cfg.trlcheck.bpos.crit      = 'y<-0.5 & z<0.2 & x<0.25 & x>0';
            cfg.trlcheck.bpos.type      = 'incl';
            % select movements
            cfg.movselect.nmov          = 2;
            cfg.movselect.param         = 'index';
         
            
            %% movparam
        case 'movparam'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task    = {'movparam'};
            %cfg.task    = {'movdef','movcheck','trlcheck','movplot'};
            memreq      = 500 * 1024^2;	% 500 MB
            timreq      = 5 * 60;       % 5 minutes
            
            % obtain movement parameters
            cfg.movparam.expfun         = 'CP_movparamfun';
            cfg.movparam.movnr          = [1 2 3];
            cfg.movparam.param          = {'rt','mt','tl','mv','pv','tpv','rtpv','nvc','pos_E','mpos','maxpos'};
            cfg.movparam.pos.marker     = {'i','hi','hp','h','w'};
            cfg.movparam.pos.axis       = {'x','y','z'};
            cfg.movparam.dist.marker    = {'i','hi','w'};
            cfg.movparam.dist.axis      = 'xyz';
            cfg.movparam.vel.marker     = {'i','hi','w'};
            cfg.movparam.vel.axis       = 'xyz';
            cfg.movparam.acc.marker     = {'i','hi','w'};
            cfg.movparam.acc.axis       = 'xyz';
            cfg.movparam.endoffset      = 0;
            
        case 'movplot'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task    = {'movplot'};
            
            % visually inspect the movements
            cfg.movplot.param           = {'pos','vel','gripapt','gripaptvel'};
            cfg.movplot.pos.marker      = {'i','hi','hp','w','h'};
            cfg.movplot.pos.axis        = {'y','x','z'};
            cfg.movplot.vel.marker      = {'i','hi','hp','w','h'};
            cfg.movplot.vel.axis        = {'yz','xyz'};
            cfg.movplot.timeref         = 'rel';
            cfg.movplot.xlim            = [-0.5 2.5];
            
        case 'trlrej'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task   	= {'trlrej'};
            %flg_exec    = 'local';
            
            % reject trials
            cfg.trlrej.fun              = 'auto';
            cfg.trlrej.apriori.fun      = 'CP_aprioritrlrejfun';
            cfg.trlrej.apriori.param	= {};
            cfg.trlrej.param            = {'corr_P1','corr_P2','rt_1','mt_1','ht','mt_2'};
            cfg.trlrej.reportparam      = {'apriori','outlier','toolate','nomov','corr_P1','corr_P2','rt_1','mt_1','ht','mt_2'};
            cfg.trlrej.nomov            = [0 1];
            cfg.trlrej.corr_P1          = [1 1];
            cfg.trlrej.corr_P2          = [1 1];
            cfg.trlrej.rt_1             = [0 8];
            cfg.trlrej.mt_1             = [0 2];
            cfg.trlrej.ht               = [0 2];
            cfg.trlrej.mt_2             = [0 2];
            
        case 'eposerr'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task    = {'eposerr'};
            memreq      = 500 * 1024^2;	% 500 MB
            timreq      = 5 * 60;     	% 5 minutes
            
            % calculate end point error
            cfg.eposerr.trlrej.dataset  = 'FILTlp15_MOV_ANA_IDXR';
            %cfg.eposerr.trlrej.listwise	= {'corr_P1','corr_P2','first2trials','toolate','nomov','rt_1'};
            cfg.eposerr.trlrej.listwise	= {'apriori','rt_1'};
            cfg.eposerr.marker          = {'i','hi','hp','w','h'};
            cfg.eposerr.axis            = 'xyz';
            cfg.eposerr.movnr           = 1;
            cfg.eposerr.movname         = 'movement';
            cfg.eposerr.endoffset       = 0;
            cfg.eposerr.fact            = {'sel_P1','action x sel_P1','action x addressee x sel_P1'};
            cfg.eposerr.lvl.goal_P1     = [4 5 6];
            cfg.eposerr.lvl.sel_P1      = [4 5 6];
            cfg.eposerr.lvl.addressee   = [1 2];
            cfg.eposerr.lvl.action      = [1 2];
            
            %% traject
        case 'traject'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task    = {'traject'};
            memreq      = 700 * 1024^2;	% 700 MB
            timreq      = 30 * 60;      % 30 minutes
            
            cfg.traject.type            = {'project'};
            cfg.traject.movname         = 'movement';
            cfg.traject.movnr           = [1]; % 1 = approach movement, 2 = backmovement
            cfg.traject.marker          = {'i','hi','hp','w','h'};
            cfg.traject.trlrej.dataset  = 'FILTlp15_MOV_ANA_IDXR';
            cfg.traject.trlrej.listwise	= {'apriori','rt_1'};
            %cfg.traject.fact            = {'sel_P1','action x sel_P1','action x addressee x sel_P1'};
            cfg.traject.fact            = {'action x sel_P1','action x time_half x sel_P1','action x addressee x sel_P1'};
            cfg.traject.lvl.goal_P1     = [4 5 6];
            cfg.traject.lvl.sel_P1      = [4 5 6];
            cfg.traject.lvl.addressee   = [1 2];
            cfg.traject.lvl.action      = [1 2];
            cfg.traject.lvl.time_half   = [1 2];
            cfg.traject.expfun          = 'no';
            cfg.traject.ntrl.sel_P1                     = [12 Inf];
            cfg.traject.ntrl.actionxsel_P1              = [12 Inf];
            cfg.traject.ntrl.actionxaddresseexsel_P1	= [8 Inf];
            cfg.traject.ntrl.actionxtime_halfxsel_P1	= [8 Inf];
            
        case 'paramreport'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task   	= {'plot'};
            flg_exec    = 'local';
            cfg.save    = 'no';
            
            % fillers
            %cfg.plot.param  = {'rt_1','mt_1','ht','mt_2','mt_P2','pos_E_i_x_1','pv_i_1','tpv_i_1','rtpv_i_1','tl_i_1','tl_i_2','pos_E_i_y_1','pos_E_i_z_1','pos_E_hi_x_1','pos_E_hi_y_1','pos_E_hi_z_1'};
            cfg.plot.param  = {'rt_1','mt_1','ht','mt_2','mt_P2'};
            
            % factor
            cfg.plot.fact   = {'action'};
            cfg.plot.lvl.action     = [1 2];    % 1:communicative, 2:instrumental
            cfg.plot.lvl.coding     = [1 2];    % 1:color, 2:spatial
            cfg.plot.lvl.addressee  = [1 2];    % 1:left, 2:right
            cfg.plot.lvl.sel_P1     = [4 5 6];  % 4:left, 5:middle, 6:right
            cfg.plot.lvl.level      = [1 2];    % 1:easy, 2:medium, 3:difficult
            cfg.plot.lvl.posterr    = [1 2];    % 1:good, 2:posterr trial
            cfg.plot.lvltick.action	= {'comm','instr'};
            cfg.plot.lvltick.coding	= {'color','spatial'};
            cfg.plot.lvltick.addressee	= {'left','right'};
            cfg.plot.lvltick.sel_P1	= {'left','middle','right'};
            cfg.plot.lvltick.level	= {'easy','medium','diff'};
            cfg.plot.lvltick.posterr= {'good','posterr'};            
            
            % descriptives: mean, median, logmean, meanlog, logmedian, medianlog
            cfg.plot.descrip.group	= 'mean';
            cfg.plot.descrip.subj	= 'mean';
            cfg.plot.descrip.err	= 'sem';
            
            % trial rejection
            cfg.plot.trlrej.dataset	= 'FILTlp15_MOV_ANA_IDXR';
            cfg.plot.trlrej.listwise= {'apriori','rt_1'};
            cfg.plot.trlrej.pairwise= cfg.plot.param;
            cfg.paramreport.trlrej.dataset	= 'FILTlp15_MOV_ANA_IDXR';
            cfg.paramreport.trlrej.listwise= {'apriori','rt_1'};
            cfg.paramreport.trlrej.pairwise= cfg.plot.param;
            
            % report settings
            cfg.paramreport.write   = {'grouprepmeas'};
            %cfg.paramreport.param   = {'all','-trlbeg','-trlend','-trloff','-RT','-MT','-GT','-HRT','-early','-late','-still','-pause','-t_*'};
            cfg.paramreport.param   = cfg.plot.param;
            cfg.paramreport.reject  = cfg.paramreport.trlrej.listwise;
            cfg.paramreport.exclude = 'cases';
            cfg.stats               = 'no';
            
        case 'plotparam'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task   	= {'plot'};
            flg_exec    = 'local';
            cfg.save    = 'no';
            % close all figures
            close all;
            
            % what to plot
            %cfg.plot.param  = {'rt_1','mt_1','ht','mt_2','mt_P2','tl_i_1','tl_i_2','pos_E_i_x_1','pos_E_i_y_1','pos_E_i_z_1','pos_E_hi_x_1','pos_E_hi_y_1','pos_E_hi_z_1'};
            %cfg.plot.param  = {'tl_i_1','pv_i_1','tpv_i_1','rtpv_i_1','pos_E_i_x_1','pos_E_i_y_1','pos_E_i_z_1'};
            cfg.plot.param  = {'pos_E_i_x_1'};
            %cfg.plot.param  = {'nvc_hi_1','tl_hi_1'};
            %cfg.plot.param  = {'pos_E_i_x_1'};
            %cfg.plot.param  = {'rt_1','mt_1','ht','nvc_hi_1','tl_hi_1','pv_hi_1','tpv_hi_1','rtpv_hi_1'};
            %cfg.plot.param  = {'pos_E_i_x_1','pos_E_i_y_1','pos_E_i_z_1','pos_E_hi_x_1','pos_E_hi_y_1','pos_E_hi_z_1'};
            %cfg.plot.param  = {'mpos_i_x_1','mpos_i_y_1','mpos_i_z_1','mpos_i_x_3','mpos_i_y_3','mpos_i_z_3'};
            % factor
            cfg.plot.fact   = {'action x sel_P1'};
            %cfg.plot.fact   = {'addresseexsel_P1'};
            %cfg.plot.fact   = {'time_half x addressee x action'};
            %cfg.plot.fact   = {'time_half x action'};
            %cfg.plot.fact   = {'posterr x addressee x action'};
            cfg.plot.lvl.action     = [1 2];    % 1:communicative, 2:instrumental
            cfg.plot.lvl.coding     = [1 2];    % 1:color, 2:spatial
            cfg.plot.lvl.addressee  = [1 2];    % 1:left, 2:right
            cfg.plot.lvl.sel_P1     = [4 5 6];  % 1:left, 2:middle, 3:right
            cfg.plot.lvl.level      = [1 2];    % 1:easy, 2:medium, 3:difficult
            cfg.plot.lvl.posterr    = [1 2];    % 1:good, 2:posterr trial
            cfg.plot.lvl.time_half  = [1 2];    % 1:first half, 2:second half
            cfg.plot.lvl.time_block = [1 2 3 4 5 6];    % block counter (4x10 trials per block of every condition)
            cfg.plot.lvltick.action	= {'comm','instr'};
            cfg.plot.lvltick.coding	= {'color','spatial'};
            cfg.plot.lvltick.addressee	= {'left','right'};
            cfg.plot.lvltick.sel_P1 = {'left','middle','right'};
            cfg.plot.lvltick.level	= {'easy','medium','diff'};
            cfg.plot.lvltick.posterr= {'good','posterr'};
            cfg.plot.lvltick.time_half   = {'firsthalf','secondhalf'};
            cfg.plot.lvltick.time_block   = {'firstblock','secondblock','thridblock','fourthblock','fifthblock','sixthblock'};
            % exclude some conditions
            %cfg.plot.factexcl       = {'sel_P1'};
            %cfg.plot.lvlexcl.sel_P1 = [5 6];
            cfg.plot.lvlexcl.sel_P1 = [];
            % how many trials/subjects should you minimally have?
            cfg.plot.reptcontrib	= 3;
            cfg.plot.subjcontrib	= 1;
            % plot types
            cfg.paramplot.type      = 'bar';
            % descriptives: mean, median, logmean, meanlog, logmedian, medianlog
            cfg.plot.descrip.group	= 'mean';
            %cfg.plot.descrip.subj	= 'mean';
            %cfg.plot.descrip.subj	= 'var';
            cfg.plot.descrip.subj	= 'median';
            %cfg.plot.descrip.subj	= 'meanlog';
            %cfg.plot.descrip.subj	= 'logmean';
            %cfg.plot.descrip.subj	= 'logmedian';
            cfg.plot.descrip.err	= 'sem';
            % trial rejection
            cfg.plot.trlrej.dataset	= 'FILTlp15_MOV_ANA_IDXR';
            cfg.plot.trlrej.listwise= {'apriori','rt_1'};
            cfg.plot.trlrej.pairwise= cfg.plot.param;
            % colors
            cfg.plot.flcl	= str2rgb('rbrbocmy');
            cfg.plot.lncl   = str2rgb('kkkkkkk');
            
        case 'plottseries'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task   	= {'plot'};
            flg_exec    = 'cluster';
            cfg.save    = 'no';
            % close all figures
            close all;
            
            % what to plot
            cfg.plot.tseries = {'pos'};
            cfg.plot.param  = {};
            % factor
            %cfg.plot.fact   = {'action'};
            %cfg.plot.fact   = {'action x sel_P1'};
            cfg.plot.fact   = {'action x addressee x sel_P1'};
            %cfg.plot.fact   = {'addressee x action'};
            %cfg.plot.fact   = {'posterr x action'};
            cfg.plot.lvl.action     = [1 2];    % 1:communicative, 2:instrumental
            cfg.plot.lvl.coding     = [1 2];    % 1:color, 2:spatial
            cfg.plot.lvl.addressee  = [1 2];    % 1:left, 2:right
            cfg.plot.lvl.sel_P1     = [4 5 6];  % 1:left, 2:middle, 3:right
            cfg.plot.lvl.level      = [1 2];    % 1:easy, 2:medium, 3:difficult
            cfg.plot.lvl.posterr    = [1 2];    % 1:good, 2:posterr trial
            cfg.plot.lvl.time_half  = [1 2];    % 1:first half, 2:second half
            cfg.plot.lvl.time_block = [1 2 3 4 5 6];    % block counter (4x10 trials per block of every condition)
            cfg.plot.lvltick.action	= {'comm','instr'};
            cfg.plot.lvltick.coding	= {'color','spatial'};
            cfg.plot.lvltick.addressee	= {'left','right'};
            cfg.plot.lvltick.sel_P1 = {'left','middle','right'};
            cfg.plot.lvltick.level	= {'easy','medium','diff'};
            cfg.plot.lvltick.posterr= {'good','posterr'};
            cfg.plot.lvltick.time_half   = {'firsthalf','secondhalf'};
            cfg.plot.lvltick.time_block   = {'firstblock','secondblock','thridblock','fourthblock','fifthblock','sixthblock'};
            % exclude some conditions
            %cfg.plot.factexcl       = {'sel_P1'};
            %cfg.plot.lvlexcl.sel_P1 = [5];
            cfg.plot.contr  = {};
            %cfg.plot.contr  = {'action'};
            
            % time series settings
            %----------------------------------------
            trlrejlst	= {'apriori','rt_1'};
            movname     = 'movement';
            cfg.tseries.vel.marker      = {'i'};
            cfg.tseries.vel.axis        = {'x'};
            cfg.tseries.vel.trlrej      = trlrejlst;
            cfg.tseries.vel.movname     = movname;
            cfg.tseries.vel.movnr       = 1;
            %cfg.tseries.pos.marker      = {'i','hi'};
            cfg.tseries.pos.marker      = {'i'};
            %cfg.tseries.pos.axis        = {'x','y','z'};
            %cfg.tseries.pos.axis        = {'x','z'};
            cfg.tseries.pos.axis        = {'x'};
            cfg.tseries.pos.trlrej      = trlrejlst;
            cfg.tseries.pos.movname     = movname;
            cfg.tseries.pos.movnr       = 1;
                        
            % how many trials/subjects should you minimally have?
            cfg.plot.reptcontrib	= 3;
            cfg.plot.subjcontrib	= 1;
            % plot types
            cfg.paramplot.type      = 'bar';
            % descriptives: mean, median, logmean, meanlog, logmedian, medianlog
            cfg.plot.descrip.group	= 'mean';
            %cfg.plot.descrip.subj	= 'mean';
            %cfg.plot.descrip.subj	= 'var';
            cfg.plot.descrip.subj	= 'median';
            %cfg.plot.descrip.subj	= 'meanlog';
            %cfg.plot.descrip.subj	= 'logmean';
            %cfg.plot.descrip.subj	= 'logmedian';
            cfg.plot.descrip.err	= 'sem';
            % trial rejection
            cfg.plot.trlrej.dataset	= 'FILTlp15_MOV_ANA_IDXR';
            cfg.plot.trlrej.listwise= {'apriori','rt_1'};
            cfg.plot.trlrej.pairwise= cfg.plot.param;
            % colors
            %cfg.plot.lnstyle	= '------';
            cfg.plot.lnstyle	= {'-','--','-','--','-','--','-','--','-','--','-','--'};
            %cfg.plot.flcl	= str2rgb('rgbocmyrgbocmy');
            %cfg.plot.flcl	= str2rgb('rbrbrb');
            cfg.plot.flcl	= str2rgb('ggbbrrbbrrbb');
            cfg.plot.lncl   = str2rgb('kkkkkkkkkkkk');
            
        case 'plottraject'
            cfg.dataset	= 'FILTlp15_MOV_ANA';
            cfg.task   	= {'plottraject'};
            flg_exec    = 'local';
            cfg.save    = 'no';
            % close all figures
            close all;
            
            % what to plot
            cfg.plot.param  = {'rt_1'};
            cfg.plot.fetchdata  = {'traject'};
            % factor --> indicate this also in CP_plot_trajectory
            %cfg.plot.fact   = {'action x addressee'};
            cfg.plot.fact   = {'action'};
            %cfg.plot.fact   = {'action x sel_P1'};
            cfg.plot.lvl.action     = [1 2];    % 1:communicative, 2:instrumental
            cfg.plot.lvl.coding     = [1 2];    % 1:color, 2:spatial
            cfg.plot.lvl.addressee  = [1 2];    % 1: left, 2:right
            cfg.plot.lvl.time_half  = [1 2];    % 1:good, 2:posterr trial
            cfg.plot.lvl.level      = [1 2];    % 1:easy, 2:medium, 3:difficult
            cfg.plot.lvl.posterr    = [1 2];    % 1:good, 2:posterr trial
            cfg.plot.lvl.sel_P1     = [4 5 6];  % 1:left, 2:middle, 3:right
            cfg.plot.lvl.goal_P1    = [4 5 6];  % 1:left, 2:middle, 3:right
            cfg.plot.lvltick.action	= {'comm','instr'};
            cfg.plot.lvltick.coding	= {'color','spatial'};
            cfg.plot.lvltick.addressee	= {'left','right'};
            cfg.plot.lvltick.time_half	= {'1st half','2nd half'};
            cfg.plot.lvltick.level	= {'easy','medium','diff'};
            cfg.plot.lvltick.posterr= {'good','posterr'};
            cfg.plot.lvltick.goal_P1= {'left','middle','right'};
            cfg.plot.lvltick.sel_P1= {'left','middle','right'};
            % how many trials/subjects should you minimally have?
            cfg.plot.reptcontrib	= 3;
            cfg.plot.subjcontrib	= 1;
            % plot types
            cfg.paramplot.type      = 'bar';
            % descriptives: mean, median, logmean, meanlog, logmedian, medianlog
            cfg.plot.descrip.group	= 'mean';
            cfg.plot.descrip.subj	= 'mean';
            %cfg.plot.descrip.subj	= 'var';
            %cfg.plot.descrip.subj	= 'median';
            %cfg.plot.descrip.subj	= 'meanlog';
            %cfg.plot.descrip.subj	= 'logmean';
            %cfg.plot.descrip.subj	= 'logmedian';
            cfg.plot.descrip.err	= 'sem';
            % trial rejection
            cfg.plot.trlrej.dataset	= 'FILTlp15_MOV_ANA_IDXR';
            cfg.plot.trlrej.listwise= {'apriori','rt_1'};
            cfg.plot.trlrej.pairwise= cfg.plot.param;
            % colors
            cfg.plot.lnstyle	= '------';
            cfg.plot.flcl	= str2rgb('rbrbocmy');
            cfg.plot.lncl   = str2rgb('kkkkkkk');
            
            
        case 'plotci'
            
            %flg = 1; % 'target';
            %flg = 2; % 'target x action';
            flg = 3; % 'target x action x addressee';
            flg_plot_subj = false;
            
            cfg.subj        = subjects(subjsel);
            cfg.sess        = sessions;
            cfg.dataset     = '_FILTlp15_MOV_ANA';
            close all;
            sel_xy = {'pos_E_i_x_1','pos_E_i_z_1'};
            %trlrejparam = {'apriori','rt_1','corr_P1','corr_P2'};
            trlrejparam = {'apriori','first2trials','toolate','nomov','corr_P1','corr_P2','rt_1'};
            code_target = 1:3; nr_target = length(code_target);
            code_action = 1;
            code_addressee = 1;
            if flg > 1
                code_action = 1:2;
            end
            if flg > 2
                code_addressee = 1:2;
            end
            nr_action = length(code_action);
            nr_addressee = length(code_addressee);
            colorlist = reshape(str2rgb('rbok')',[3 2 2]);
            colorlist_light = colorlist + 3*(1-colorlist)/4;
            line_addressee = {'-','-'};
            xy_agg = {nan(0,2)};
            xy_agg = repmat(xy_agg,[nr_target nr_action nr_addressee]);
            CI = cell(nr_target,nr_action,nr_addressee,length(subjsel));
            mu = nan(2,nr_target,nr_action,nr_addressee,length(subjsel));
            area = nan(nr_target,nr_action,nr_addressee,length(subjsel));
            d = 50;
            for s = 1:length(subjsel)
                
                tdcfg = load(fullfile(cfg.dir.ana,sprintf('%s%s%s_CFG',cfg.subj{s},cfg.sess{1},cfg.dataset)));
                dcfg = tdcfg.cfg;
                
                tdidxr = load(fullfile(cfg.dir.ana,sprintf('%s%s%s_IDXR',cfg.subj{s},cfg.sess{1},cfg.dataset)));
                didxr = tdidxr.cfg;
                
                % fetch data
                target = dcfg.trl(:,ismember(dcfg.vars,'sel_P1')) - 3;
                action = dcfg.trl(:,ismember(dcfg.vars,'action'));
                addressee = dcfg.trl(:,ismember(dcfg.vars,'addressee'));
                if flg < 2
                    action = ones(size(action));
                end
                if flg < 3
                    addressee = ones(size(addressee));
                end
                %coding = dcfg.trl(:,ismember(dcfg.vars,'coding'));
                %level  = dcfg.trl(:,ismember(dcfg.vars,'level'));
                dat = dcfg.trl(:,ismember(dcfg.vars,sel_xy));
                
                % replace rejection trials by NaNs
                didxr = km_combidxr(didxr,trlrejparam);
                dat(didxr.reject,:) = NaN;
                
                % initialize figure
                if flg_plot_subj
                    figure(s); clf; hold on; hl = nan(1,nr_action); hl2 = nan(nr_action,nr_addressee);
                end
                
                % loop over targets and actions
                for t = code_target
                    for a = code_action
                        for ad = code_addressee
                            % get data
                            idx = target==t & action==a & addressee==ad;
                            xy = dat(idx,:);
                            xy = xy(~any(isnan(xy),2),:);
                            % aggregate data over subjects
                            xy_agg{t,a,ad} = [xy_agg{t,a,ad}; xy];
                            % plot touch points
                            if flg_plot_subj
                                scatter(xy(:,1),xy(:,2),d,colorlist_light(:,a,ad)','filled');
                            end
                            % get confidence ellipse
                            CI{t,a,ad,s} = km_mdCI(xy);
                            mu(:,t,a,ad,s) = CI{t,a,ad,s}.mu;
                            area(t,a,ad,s) = CI{t,a,ad,s}.xy.area;
                        end
                    end
                    % plot CI on top of dots
                    if flg_plot_subj
                        for a = code_action
                            for ad = code_addressee
                                hl(a) = scatter(mu(1,t,a,ad,s),mu(2,t,a,ad,s),d,colorlist(:,a,ad)','filled');
                                hl2(a,ad) = plot(CI{t,a,ad,s}.xy.xyz(:,1),CI{t,a,ad,s}.xy.xyz(:,2),'marker','none','linewidth',2,'color',colorlist(:,a,ad)','linestyle',line_addressee{ad});
                            end
                        end
                    end
                end
                if flg_plot_subj
                    %axis([-0.15 0.05 0.3 0.6]);
                    axis([-0.05 0.2 -0.75 -0.45]);
                    axis equal;
                    if flg == 2
                        legend(hl(:),'comm','instr','location','NorthOutside');
                    elseif flg == 3
                        legend(hl2(:),'comm-left','instr-left','comm-right','instr-right','location','NorthOutside');
                    end
                    %legend(hl2(:),'comm-left','comm-right','location','NorthOutside');
                end
            end
            
            % average over subject
            mu_avg = nanmean(mu,5);
            area_avg = nanmean(area,4);
            
            % plot mu
            figure(s+1); clf; hold on; 
            if flg == 2
                hl = nan(2,1);
            elseif flg == 3 
                hl = nan(2,2);
            end
            for a = code_action
                for ad = code_addressee
                    hl(a,ad) = scatter(mu_avg(1,:,a,ad),mu_avg(2,:,a,ad),d,colorlist(:,a,ad)','filled');
                end
            end
            axis([-0.15 0.05 0.5 0.55]);
            axis equal;
            if flg == 2
                legend(hl(:),'comm','instr','location','NorthOutside');
            elseif flg == 3
                legend(hl(:),'comm-left','instr-left','comm-right','instr-right','location','NorthOutside');
            end
            %legend(hl2(:),'comm-left','comm-right','location','NorthOutside');
            
                    % plot area
                    figure(s+2); clf; hold on;
                    for t = code_target
                        subplot(1,3,t);
                        hb = bar(area_avg(t,:));
                        % set individual bar face colors
                        hp = get(hb,'Children');
                        set(hp,'CDataMapping','direct','CData',1:2)
                        %colormap(color_action);
                        axis([0.5 2.5 0 0.0015]);
                    end
            
            % now analyse the data aggregated over all subjects
            CI_agg = cell(9,2,2);
            mu_agg = nan(2,9,2,2);
            area_agg = nan(9,2,2);
            % initialize figure
            figure(s+3); clf; hold on; hl = nan(1,2);
            for t = code_target
                for a = code_action
                    for ad = code_addressee
                        % get xy
                        xy = xy_agg{t,a,ad};
                        % plot touch points
                        scatter(xy(:,1),xy(:,2),d,colorlist_light(:,a,ad)','filled');
                        % get confidence ellipse
                        CI_agg{t,a,ad} = km_mdCI(xy);
                        mu_agg(:,t,a,ad) = CI_agg{t,a,ad}.mu;
                        area_agg(t,a,ad) = CI_agg{t,a,ad}.xy.area;
                    end
                end
                % plot CI on top of dots
                for a = code_action
                    for ad = code_addressee
                        hl(a) = scatter(mu_agg(1,t,a,ad),mu_agg(2,t,a,ad),d,colorlist(:,a,ad)','filled');
                        hl2(a,ad) = plot(CI_agg{t,a,ad}.xy.xyz(:,1),CI_agg{t,a,ad}.xy.xyz(:,2),'marker','none','linewidth',2,'color',colorlist(:,a,ad)','linestyle',line_addressee{ad});
                    end
                end
                
            end
            axis([-0.15 0.05 0.25 0.55]);
            axis equal;
            if flg == 2
                %legend(hl(:),'comm','instr','location','NorthOutside');
            elseif flg == 3
                legend(hl2(:),'comm-left','instr-left','comm-right','instr-right','location','NorthOutside');
            end
            %legend(hl2(:),'comm-left','comm-right','location','NorthOutside');
            
%             %save
%             dir_report = fullfile(filesep,'home','action','ankmur','CommPointing','tempplots_kin');
%             if ~exist(dir_report,'dir'), mkdir(dir_report); end
%             figname = fullfile(dir_report,'CI_allsubj');
%             export_fig(figname,'-pdf','openGL');
%             %export_fig(figname,'-png','-r600');
%             
            
            return
            
        otherwise
            error('option not supported');
            
    end
    
    % execute immediately if requested
    if isequal(cfg.task,{'plot'})
        % who to plot
        cfg.subj        = subjects(subjsel);
        cfg.sess        = sessions;
        % execute plotting wrapper
        close all;
        % call plotting wrapper
        km_plot(cfg);
        % return
        return
    elseif isequal(cfg.task,{'plottraject'})
        % who to plot
        cfg.subj        = subjects(subjsel);
        cfg.sess        = sessions;
        % execute plotting wrapper
        close all;
        % call plotting wrapper
        CP_plot_trajectory(cfg);
        % return
        return
    end
    
    
    %% execute tasks
    % clear warnings and errors
    lastwarn(''); %lasterr('');
    % loop over subjects
    tcfg = cell(1,length(subjsel));
    for s = 1:length(subjsel)
        
        % initialize tcfg with the default cfg fields
        tcfg{s}             = cfg;
        
        % setup cfg for each subject individually
        tcfg{s}.subj        = subjects{subjsel(s)};
    end
    
    % execute the task either locally or distribute to the cluster
    if length(subjsel) == 1 || strcmpi(flg_exec,'local')
        % execute single subject tasks locally
        cellfun(@km_distribute_tasks,tcfg,'UniformOutput',false);
    else
        % or distribute each subject separately to the cluster
        qsubcellfun(@km_distribute_tasks,tcfg,'memreq',memreq,'timreq',timreq,'StopOnError',false);
    end
    
    % end settask-loop
end