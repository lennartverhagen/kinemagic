function varargout = km_movplot(cfg,data)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% retrieve input if not supplied
if nargin < 1
	[fname,pname] = uigetfile('*CFG.mat','Select stored CFG file');
    load(fullfile(pname,fname));
end
if nargin < 2
	[fname,pname] = uigetfile('*DATA.mat','Select stored DATA file');
    load(fullfile(pname,fname));
    if nargin < 1 || ~isfield(cfg,'movplot')
        basicparam = km_labelselection({'pos','ori'},fieldnames(data));
        for p = 1:length(basicparam)
            cfg.movplot.(basicparam{p}).axis = data.axis;
            cfg.movplot.(basicparam{p}).marker = data.label;
        end
        derivparam = km_labelselection({'vel','acc','orivel','oriacc'},fieldnames(data));
        for p = 1:length(derivparam)
            cfg.movplot.(derivparam{p}).axis = data.derivaxis;
            cfg.movplot.(derivparam{p}).marker = data.label;
        end
        if isfield(data,'grip')
            gripaptparam = km_labelselection({'apt','aptvel','aptacc'},fieldnames(data.grip));
            gripaptparam = cellfun(@(str) ['grip' str],gripaptparam,'UniformOutput',false);
            for p = 1:length(gripaptparam)
                cfg.movplot.(gripaptparam{p}).axis = data.grip.aptaxis;
                cfg.movplot.(gripaptparam{p}).marker = data.grip.label;
            end
            griporiparam = km_labelselection({'ori','orivel','oriacc'},fieldnames(data.grip));
            griporiparam = cellfun(@(str) ['grip' str],griporiparam,'UniformOutput',false);
            for p = 1:length(gripaptparam)
                cfg.movplot.(griporiparam{p}).axis = data.grip.oriaxis;
                cfg.movplot.(griporiparam{p}).marker = data.grip.label;
            end
        else
            gripaptparam = {}; griporiparam = {};
        end
        cfg.movplot.param = [basicparam; derivparam; gripaptparam; griporiparam];
        cfg.movplot.xlim = [-0.5 3];
    end
end

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end

% FIXME: legacy reasons for app > apt
if isfield(data,'grip')
    if isfield(data.grip,'appaxis')
        data.grip.aptaxis = data.grip.appaxis;
        data.grip = rmfield(data.grip,'appaxis');
    end
    if isfield(data.grip,'app')
        data.grip.apt = data.grip.app;
        data.grip = rmfield(data.grip,'app');
    end
    if isfield(data.grip,'appvel')
        data.grip.aptvel = data.grip.appvel;
        data.grip = rmfield(data.grip,'appvel');
    end
    if isfield(data.grip,'appacc')
        data.grip.aptacc = data.grip.appacc;
        data.grip = rmfield(data.grip,'appacc');
    end
end


% build info structure on which the figure is based
info = build_info(cfg,data);

% build graphical user interface
hf = build_gui(info);

% plot the data
draw_dat(hf);
draw_mov(hf);

% run loop until quit
while ishandle(hf)
    redraw(hf);
    info = guidata(hf);
    if ~info.quit
        uiwait;
    else
        delete(hf);
        break
    end
end

% update dataset
if isempty(regexpi(cfg.dataset,'_mov','once'))
    if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
    cfg.dataset = [cfg.dataset '_MOV'];
end

% update processing directory
cfg = km_setcfg(cfg,'dirana');
cfg.dir.proc = cfg.dir.ana;

% return new movement selection or store in data and configuration structure
if nargout == 1
    varargout = {info.selmov};
else
    if isfield(cfg,'movdef')
        for i = 1:length(cfg.movdef)
            cfg.movdef(i).movement = info.selmov;
        end
    end
    cfg.movement	= info.selmov;
    data.movement	= info.selmov;
    varargout       = {cfg,data};
end
%--------------------------------------------------------------------------


function info = build_info(cfg,data)
%----------------------------------------

% get movement
if isfield(cfg,'movement') && ~isstruct(cfg.movement) && ~isempty(cfg.movement)
    movement = cfg.movement;
elseif isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
    movement = cfg.movdef.movement;
elseif isfield(data,'movement')
    movement = data.movement;
elseif isfield(cfg,'movdef')
    movement = km_movdef(cfg,data);
    if isfield(cfg,'movcheck')
        cfg.movement = movement;
        movement = km_movcheck(cfg,data);
    end
    if isfield(cfg,'trlcheck')
        cfg.movement = movement;
        movement = km_trlcheck(cfg,data);
    end
else
    movement = repmat({zeros(0,2)},1,length(data.time));
end

% if the data has been analyzed already
allmov = movement;
if ~isempty(regexpi(cfg.dataset,'_ana$','once'))
    selmov = movement;
else
    selmov = repmat({zeros(0,2)},1,length(movement));
end

% maybe some movements have been pre-selected by km_movcheck
if isfield(cfg,'movcheck') && isfield(cfg.movcheck,'allmov');
    allmov = cfg.movcheck(1).allmov;
    selmov = movement;
end

% initialize info structure for figure
info = struct;
info.quit = false;

% data is continuous or not
nallmov = cellfun(@(x) size(x,1),allmov);
nselmov = cellfun(@(x) size(x,1),selmov);
if isfield(cfg.movplot,'continuous')
    info.continuous = istrue(cfg.movplot.continuous,'free');
elseif isfield(cfg,'continuous')
    info.continuous = istrue(cfg.continuous,'free');
else
    mintim = cellfun(@min,data.time);
    maxtim = cellfun(@max,data.time);
    if (all(nallmov<10) && all(nselmov<3)) || (length(data.time) > 1 && all(mintim <= min(maxtim)))
        info.continuous = false;
    else
        info.continuous = true;
    end
end

% relative or absolute timing
if isfield(cfg.movplot,'timeref') && ischar(cfg.movplot.timeref)
    info.absreltim = cfg.movplot.timeref;
elseif isfield(cfg,'timeref') && ischar(cfg.timeref)
    info.absreltim = cfg.timeref;
end

% get trial selection
if ~info.continuous && isfield(cfg,'vars') && isfield(cfg,'trl')
    if isfield(cfg,'cuttrials') && isfield(cfg.cuttrials,'timevar')
        timevar = cfg.cuttrials.timevar;
    else
        timevar = '^t_';
    end
    idx = ~cellfun(@isempty,regexp(cfg.vars,timevar));
    info.vars = regexprep(cfg.vars(idx),timevar,'');
    info.trl = cfg.trl(:,idx);
else
    info.vars = {};
    info.trl = [];
end

% basic selections
info.r = 1; % run/trial
info.i = 1; % mov
info.m = 1; % marker
info.a = 1;	% axis
info.p = 1; % parameter
info.currr = [];	% currently plotted run
info.curri = []; 	% currently plotted mov
info.previ = []; 	% previously plotted mov

% movement onsets and offsets
info.checkparam = {};
info.check      = repmat({[]},1,length(allmov));
info.showchecks = false;
%if isempty(regexpi(cfg.dataset,'_ana$','once'))
if isfield(cfg,'movcheck')
    info.checkparam = cfg.movcheck(1).param;
    info.check      = cfg.movcheck(1).check;
    np = length(info.checkparam);
    nc = length(cfg.movcheck);
    if nc > 1
        ic = repmat(1:nc,[np 1]);
        info.checkparam = cellfun(@(x,y) sprintf('%s_%i',x,y),repmat(info.checkparam',[1 nc]),num2cell(ic),'UniformOutput',false);
        info.checkparam = info.checkparam(:)';
        d = np*nc;
        info.check = cellfun(@(x) reshape(x,[size(x,1) d]),info.check,'UniformOutput',false);
    end
    info.showchecks = true;
end
if isfield(cfg,'trlcheck')
    trlcheckparam   = cellfun(@(x) ['trl_' x],cfg.trlcheck.param,'UniformOutput',false);
    info.checkparam = [info.checkparam trlcheckparam];
    trlcheck        = mat2cell(cfg.trlcheck.check,ones(length(info.check),1))';
    info.check      = cellfun(@(x,y) [x repmat(y,[size(x,1) 1])],info.check,trlcheck,'UniformOutput',false);
    info.showchecks = true;
end
%end
info.allmov = allmov;
info.selmov = selmov;
info.nallmov = nallmov;
info.nselmov = nselmov;

% store time and time limits
info.nrun	= length(data.time);
info.time   = data.time;
info.xlim   = cfg.movplot.xlim;

% loop over parameters and find markers, axes and data
info.param = cfg.movplot.param;
info.nparam = length(info.param);
for p = 1:info.nparam
    
    % parameter name
    pname = info.param{p};
    
    % get marker list
    if ~isempty(regexp(pname,'^grip','once'))
        pmarker = data.grip.label;
    else
        pmarker = data.label;
    end
    % get axis list
    if ~isempty(regexp(pname,'(?<!^grip.*)(vel$|acc$)','once'))
        paxis = data.derivaxis;
    elseif ~isempty(regexp(pname,'^gripori','once'))
        paxis = data.grip.oriaxis;
    elseif ~isempty(regexp(pname,'^grip','once'))
        paxis = data.grip.aptaxis;
    else
        paxis = data.axis;
    end
    % get data
    if ~isempty(regexp(pname,'^gripori','once'))
        dat      = data.grip.(pname(5:end));
    elseif ~isempty(regexp(pname,'^grip','once'))
        dat      = data.grip.(pname(5:end));
    else
        dat      = data.(pname);
    end
    
    % store data
    info.(pname).dat    = dat;
    
    % store markers and axis
    info.(pname).marker = pmarker;
    info.(pname).midx   = 1:length(pmarker);
    info.(pname).axis   = paxis;
    info.(pname).aidx   = 1:length(paxis);
    
    % store axes limits
    if isfield(cfg.movplot,pname) && isfield(cfg.movplot.(pname),'ylim')
        info.(pname).ylim = cfg.movplot.(pname).ylim;
    else
        info.(pname).ylim = cfg.movplot.ylim;
    end
    
    % try to find a slide criterium
    if isfield(cfg.movplot,pname) && isfield(cfg.movplot.(pname),'slide')
        info.(pname).slide = cfg.movplot.(pname).slide;
    elseif isfield(cfg.movplot,'slide')
        info.(pname).slide = cfg.movplot.slide;
    elseif isfield(cfg,'movdef') && length(cfg.movdef)==1 && isfield(cfg.movdef,pname) && isfield(cfg.movdef.(pname),'slide')
        info.(pname).slide = cfg.movdef.(pname).slide;
    else
        info.(pname).slide = 0;
    end
    
    if p > 2, continue; end
    
    % select marker from list
    if isfield(cfg.movplot,pname) && isfield(cfg.movplot.(pname),'marker')
        markersel = cfg.movplot.(pname).marker;
    elseif isfield(cfg.movplot,'marker')
        markersel = cfg.movplot.marker;
    elseif isfield(cfg.movdef,'pname') && isfield(cfg.movdef.(pname),'marker')
        markersel = cfg.movdef.(pname).marker;
    else
        markersel = pmarker{1};
    end
    markersel = find(ismember(pmarker,markersel),1,'first');
    if isempty(markersel),  markersel = 1;  end
    % select axis from list
    if isfield(cfg.movplot,pname) && isfield(cfg.movplot.(pname),'axis')
        axissel = cfg.movplot.(pname).axis;
    elseif isfield(cfg.movplot,'axis')
        axissel = cfg.movplot.axis;
    elseif isfield(cfg.movdef,'pname') && isfield(cfg.movdef.(pname),'axis')
        axissel = cfg.movdef.(pname).axis;
    else
        axissel = paxis{1};
    end
    axissel = find(ismember(paxis,axissel),1,'first');
    if isempty(axissel),  axissel = 1;  end
    
    % store selections
    info.m(p) = markersel;
    info.a(p) = axissel;
    info.p(p) = p;
end

% if km_movplot is called by km_movdef, capture movement windows
if strcmpi(cfg.dataset,'movdef');
    info.plotmovdef = true;
    if ~info.continuous
        warning('KM:km_movplot:continuous','It seems like your data is not (expected to be) continuous. This might give problems when plotting movement definitions (as you are doing now).');
    end
    
    % prevent xlim from being negative
    info.xlim(info.xlim<0) = 0;
    % if diff(info.xlim) < 1, info.xlim(2) = info.xlim(1) + 1;    end
    
    info.movdef.param = cfg.movplot.movdef.param;
    info.movdef.win = data.movwin;
    info.movdef.color = repmat('gcmyo',1,ceil(length(info.movdef.param)/5));
    info.movdef.color = num2cell(info.movdef.color(1:length(info.movdef.param)));
    info.movdef.idx.indiv = 1:0;
    info.movdef.idx.comb = 1:length(info.movdef.param);
else
    info.plotmovdef = false;
    info.movdef.idx.indiv = nan(1,0);
    info.movdef.idx.comb = nan(1,0);
end

% set default movwin parameter indices (only used for movdef)
function hf = build_gui(info)
%----------------------------------------
% create new figure
% default position/size: [x y 560 420];
hf = figure(...
    'Visible','off',...
    'WindowStyle','normal',...
    'DockControls','off',...
    'Position',[200 200 600 480],...
    'Color','w',...
    'NumberTitle','off',...
    'Name','Movement Selection');
if info.plotmovdef,	set(hf,'Name','Movement Definition');   end

%% custom pointer
%----------------------------------------
set(hf,'PointerShapeCData',pointer_lrdrag('Cdata'),'PointerShapeHotSpot',pointer_lrdrag('HotSpot'))

%% axes
%----------------------------------------
info.hf = hf;
info.ha(1) = axes('units','normalized',...
    'position',[0.1 0.2 0.6 0.7],...
    'LineWidth',2,...
    'Color','none');
if info.nparam > 1
    info.ha(2) = axes('units','normalized',...
        'position',get(info.ha(1),'Position'),...
        'XAxisLocation','bottom',...
        'YAxisLocation','right',...
        'Color','none',...
        'YColor','b',...
        'LineWidth',2,...
        'XTick',[]);
end

%% trial markers
%----------------------------------------
info.htrl.line	= -1;
info.htrl.text	= -1;

%% movement parameters (lines) and selections (patches)
%----------------------------------------
info.hl         = -ones(size(info.ha));
info.hm.all     = -1;
info.hm.sel     = -1;
info.hm.text    = -1;
info.hm.check   = -1;

%% movement definitions (lines)
%----------------------------------------
if info.plotmovdef
    info.hmdef	= -ones(1,length(info.movdef.param) + 1);
else
    info.hmdef  = -1;
end

%% buttons and text panels
%----------------------------------------
% define strings
if info.continuous, runstr = 'run';
else                runstr = 'trial';   end
rstr        = num2str(info.r);
istr        = num2str(info.i);
absrelstr   = {'abs','rel'};
param       = info.param;
xlimstr     = get_limstr(info.xlim);
ylim1str    = get_limstr(info.(param{1}).ylim);
if info.nparam > 1
    ylim2str= get_limstr(info.(param{2}).ylim);
end
if info.plotmovdef
    movdefparamC= info.movdef.param;
    movdefparamI= cellfun(@(a,b) [a '-' b],info.movdef.color,movdefparamC,'UniformOutput',false);
    set(info.ha,'Units','pixels');
    v = get(info.ha(1),'Position');
    set(info.ha,'Units','normalized');
    tlim = [min(info.time{info.r}) max(info.time{info.r})];
    xwin = diff(info.xlim); twin = diff(tlim);
end

% done/quit button
x = 5;
info.hb.quit        = uicontrol(hf,                     'position',[x      5  40  18],'String','done',    'Backgroundcolor','w','Callback',@b_quit);
x = x + 40 + 10;
% run/trial selection
if info.nrun > 10
info.hb.tjprev      = uicontrol(hf,                     'position',[x      5  25  18],'String','<<',      'Backgroundcolor','w','Callback',@b_tjprev);
info.hb.tjnext      = uicontrol(hf,                     'position',[x+110  5  25  18],'String','>>',      'Backgroundcolor','w','Callback',@b_tjnext);
info.ht.trial       = uicontrol(hf,'Style','text',      'position',[x     23 135  18],'String',runstr,    'Backgroundcolor','w');
x = x + 25;
else info.ht.trial  = uicontrol(hf,'Style','text',      'position',[x     23  85  18],'String',runstr,    'Backgroundcolor','w');
end
info.hb.tprev       = uicontrol(hf,                     'position',[x      5  25  18],'String','<',       'Backgroundcolor','w','Callback',@b_tprev);
info.hb.tedit       = uicontrol(hf,'Style','edit',      'position',[x+25   5  35  18],'String',rstr,      'Backgroundcolor','w','Callback',@e_tedit);
info.hb.tnext       = uicontrol(hf,                     'position',[x+60   5  25  18],'String','>',       'Backgroundcolor','w','Callback',@b_tnext);
if info.nrun > 10,  x = x + 110 + 15;
else                x = x + 85 + 15;  end
% movement selection
if info.continuous && ~info.plotmovdef
if any(info.nallmov > 10)
info.hb.mjprev      = uicontrol(hf,                     'position',[x      5  25  18],'String','<<',      'Backgroundcolor','w','Callback',@b_mjprev);
info.hb.mjnext      = uicontrol(hf,                     'position',[x+110  5  25  18],'String','>>',      'Backgroundcolor','w','Callback',@b_mjnext);
info.ht.mov         = uicontrol(hf,'Style','text',      'position',[x     23 135  18],'String','mov',     'Backgroundcolor','w');
x = x + 25;
else info.ht.mov    = uicontrol(hf,'Style','text',      'position',[x     23  85  18],'String','mov',     'Backgroundcolor','w');
end
info.hb.mprev       = uicontrol(hf,                     'position',[x      5  25  18],'String','<',       'Backgroundcolor','w','Callback',@b_mprev);
info.hb.medit       = uicontrol(hf,'Style','edit',      'position',[x+25   5  35  18],'String',istr,      'Backgroundcolor','w','Callback',@e_medit);
info.hb.mnext       = uicontrol(hf,                     'position',[x+60   5  25  18],'String','>',       'Backgroundcolor','w','Callback',@b_mnext);
if any(info.nallmov > 10),  x = x + 110 + 15;
else                        x = x + 85 + 15;  end
end
% axes limits
if info.plotmovdef
info.ht.time        = uicontrol(hf,'Style','text',      'position',[x     23  80  18],'String','abs time','Backgroundcolor','w');
info.hb.xslider     = uicontrol(hf,'Style','slider',    'position',[v(1)  50 v(3) 20],                    'Backgroundcolor','w','Callback',@s_xslider);
set(info.hb.xslider,'Min',tlim(1),'Max',tlim(2),'Value',xwin/2,'SliderStep',[xwin/(5*twin) xwin/twin]);
else
info.hb.absreltim   = uicontrol(hf,'Style','popupmenu', 'position',[x     26  40  20],'String',absrelstr, 'Backgroundcolor','w','Callback',@p_absreltim);
if ~isfield(info,'absreltim')
    info.absreltim = 'abs';
    if info.continuous || ~all(cellfun(@(x) min(x)<=0 && max(x)>=0,info.time))
        info.absreltim = 'rel';
    end
end
if strcmpi(info.absreltim,'rel')
    set(info.hb.absreltim,'Value',2);
else
    set(info.hb.absreltim,'Value',2);
end
info.ht.time        = uicontrol(hf,'Style','text',      'position',[x+40  23  35  18],'String','time',    'Backgroundcolor','w');
end
info.hb.xprev       = uicontrol(hf,                     'position',[x      5  25  18],'String','<',       'Backgroundcolor','w','Callback',@b_xprev);
info.hb.xlim        = uicontrol(hf,'Style','edit',      'position',[x+25   5  80  18],'String',xlimstr,   'Backgroundcolor','w','Callback',@b_xlim);
info.hb.xnext       = uicontrol(hf,                     'position',[x+105  5  25  18],'String','>',       'Backgroundcolor','w','Callback',@b_xnext);
x = x + 130 + 15;
% show checks
if ~info.plotmovdef
info.hb.checks      = uicontrol(hf,'Style','togglebutton','position',[x    5  80  18],'String','no checks','Backgroundcolor','w','Callback',@b_checks);
set(info.hb.checks,'Value',get(info.hb.checks,'Max'));
if isempty(info.checkparam) || all(cellfun(@isempty,info.check))
    set(info.hb.checks,'String','checks');
    set(info.hb.checks,'Enable','off');
end
x = x + 80 + 15;
end

%% popupmenu parameter selection
%----------------------------------------
% first parameter panel
y = 330;
m1str = info.(param{1}).marker;
a1str = info.(param{1}).axis;
info.hp.param(1)    = uipanel(hf,  'units','pixels',    'position',[490 y    100 115],'Title','Param 1',                    'BackgroundColor','w');
info.hb.param(1)    = uicontrol(hf,'Style','popupmenu', 'position',[500 y+80  80  20],'String',param,   'Value',1,          'Backgroundcolor','w','Callback',{@p_param,1});
info.hb.marker(1)   = uicontrol(hf,'Style','popupmenu', 'position',[500 y+55  80  20],'String',m1str,   'Value',info.m(1),  'Backgroundcolor','w','Callback',{@p_marker,1});
info.hb.axis(1)     = uicontrol(hf,'Style','popupmenu', 'position',[500 y+30  80  20],'String',a1str,   'Value',info.a(1),  'Backgroundcolor','w','Callback',{@p_axis,1});
info.hb.ylim(1)     = uicontrol(hf,'Style','edit',      'position',[525 y+5   55  20],'String',ylim1str,                    'Backgroundcolor','w','Callback',{@b_ylim,1});
info.ht.ylim(1)     = uicontrol(hf,'Style','text',      'position',[500 y+5   25  20],'String','lim',                       'Backgroundcolor','w');

% second parameter panel
if info.nparam > 1
y = 190;
m2str = info.(param{2}).marker;
a2str = info.(param{2}).axis;
info.hp.param(2)    = uipanel(hf,  'units','pixels',    'position',[490 y    100 115],'Title','Param 2',                    'BackgroundColor','w','ForegroundColor','b');
info.hb.param(2)    = uicontrol(hf,'Style','popupmenu', 'position',[500 y+80  80  20],'String',param,   'Value',2,          'Backgroundcolor','w','Callback',{@p_param,2});
info.hb.marker(2)   = uicontrol(hf,'Style','popupmenu', 'position',[500 y+55  80  20],'String',m2str,   'Value',info.m(2),  'Backgroundcolor','w','Callback',{@p_marker,2});
info.hb.axis(2)     = uicontrol(hf,'Style','popupmenu', 'position',[500 y+30  80  20],'String',a2str,   'Value',info.a(2),  'Backgroundcolor','w','Callback',{@p_axis,2});
info.hb.ylim(2)     = uicontrol(hf,'Style','edit',      'position',[525 y+5   55  20],'String',ylim2str,                    'Backgroundcolor','w','Callback',{@b_ylim,2});
info.ht.ylim(2)     = uicontrol(hf,'Style','text',      'position',[500 y+5   25  20],'String','lim',                       'Backgroundcolor','w');
end

if info.plotmovdef
%% movdef control
%----------------------------------------
y = 10;
info.hp.movdef      = uipanel(hf,  'units','pixels',    'position',[490 y    100 170],'Title','MovDef',                     'BackgroundColor','w','ForegroundColor','r');
info.ht.movdef.indiv= uicontrol(hf,'Style','text',      'position',[500 y+140 80  15],'String','plot indiv',                'Backgroundcolor','w');
info.hb.movdef.indiv= uicontrol(hf,'Style','listbox',   'position',[500 y+95  80  45],'String',movdefparamI,                'Backgroundcolor','w');
info.ht.movdef.comb = uicontrol(hf,'Style','text',      'position',[500 y+75  80  15],'String','combined',                  'Backgroundcolor','w');
info.hb.movdef.comb = uicontrol(hf,'Style','listbox',   'position',[500 y+30  80  45],'String',movdefparamC,                'Backgroundcolor','w');
info.hb.movdef.plot = uicontrol(hf,                     'position',[500 y+5   80  20],'String','plot',                      'Backgroundcolor','w','Callback',@b_movdef);
if length(movdefparamC) > 1
set(info.hb.movdef.indiv,'Min',0,'Max',length(movdefparamC),'Value',info.movdef.idx.indiv);
set(info.hb.movdef.comb ,'Min',0,'Max',length(movdefparamC),'Value',info.movdef.idx.comb);
end
else
%% slide control
%----------------------------------------
y = 30;
slidestr = sprintf('%g',info.(param{1}).slide);
info.hp.slide       = uipanel(hf,  'units','pixels',    'position',[490 y    100 135],'Title','Slide',                      'BackgroundColor','w','ForegroundColor','r');
info.hb.slide.slide = uicontrol(hf,                     'position',[500 y+100 80  20],'String','slide',                     'Backgroundcolor','w','Callback',@b_slide);
info.hb.slide.param = uicontrol(hf,'Style','popupmenu', 'position',[500 y+80  80  20],'String',param,   'Value',1,          'Backgroundcolor','w','Callback',@p_slideparam);
info.hb.slide.marker= uicontrol(hf,'Style','popupmenu', 'position',[500 y+55  80  20],'String',m1str,   'Value',info.m(1),  'Backgroundcolor','w');
info.hb.slide.axis  = uicontrol(hf,'Style','popupmenu', 'position',[500 y+30  80  20],'String',a1str,	'Value',info.a(1),  'Backgroundcolor','w');
info.hb.slide.crit  = uicontrol(hf,'Style','edit',      'position',[530 y+5   50  20],'String',slidestr,                    'Backgroundcolor','w','Callback',@b_slidecrit);
info.ht.slide       = uicontrol(hf,'Style','text',      'position',[500 y+5   30  20],'String','crit',                      'Backgroundcolor','w');
set(info.hb.slide.slide,'Enable','off');
end

%% cursor control
%----------------------------------------
set(hf, 'WindowButtonDownFcn',  {@m_select, 'event', 'ButtonDown',  'callback', {@m_select_cb, hf}});
set(hf, 'WindowButtonUpFcn',    {@m_select, 'event', 'ButtonUp',    'callback', {@m_select_cb, hf}});
set(hf, 'WindowButtonMotionFcn',{@m_select, 'event', 'Motion',      'callback', {@m_select_cb, hf}});

%% wrap up
%----------------------------------------
% store info with the figure for callback functions
guidata(hf,info);

% Allow resizing
set(info.hf,'Units','normalized');
set(info.ha,'Units','normalized');
flds = fieldnames(info.hb);
for i = 1:length(flds)
    if isstruct(info.hb.(flds{i}))
        subflds = fieldnames(info.hb.(flds{i}));
        for j = 1:length(subflds)
            set(info.hb.(flds{i}).(subflds{j}),'Units','normalized');
        end
    else
        set(info.hb.(flds{i}),'Units','normalized');
    end
end
flds = fieldnames(info.ht);
for i = 1:length(flds)
    if isstruct(info.ht.(flds{i}))
        subflds = fieldnames(info.ht.(flds{i}));
        for j = 1:length(subflds)
            set(info.ht.(flds{i}).(subflds{j}),'Units','normalized');
        end
    else
        set(info.ht.(flds{i}),'Units','normalized');
    end
end
flds = fieldnames(info.hp);
for i = 1:length(flds)
    set(info.hp.(flds{i}),'Units','normalized');
end
% move the GUI to the centre
movegui(info.hf,'center')
% make the GUI visible
set(info.hf,'Visible','on');


%% function draw_dat
%----------------------------------------
function draw_dat(h,ai)
% make sure the handle is a figure handle
if ~ishandle(h),    return; end
while ~isequal(get(h,'parent'),0)
    h = get(h,'parent');
end
info = guidata(h);
if nargin < 2,  ai = 1:length(info.p);    end

% get patch handles
htl = info.htrl.line;
htt = info.htrl.text;
   
% erase old patches
for jj = 1:length(htl)
    if ishandle(htl(jj)),	delete(htl(jj)); end
end
for jj = 1:length(htt)
    if ishandle(htt(jj)),	delete(htt(jj)); end
end

% get variables from info
p = info.p(ai);
r = info.r;

% select time
tim = info.time{r};

% get xlim
if strcmpi(info.absreltim,'abs') || ...
    (strcmpi(info.absreltim,'rel') && info.i > 0 && ~isempty(info.allmov{info.r}))
    tmpxlim = set_xlim(info);
    idx = (tim >= tmpxlim(1)-1) & (tim <= tmpxlim(2)+1);
    tim = tim(idx);
    if isempty(tim)
        % delete line plots
        for i = 1:length(p)
            if ishandle(info.hl(ai(i))),    delete(info.hl(ai(i))); end
        end
        % remove old movdef lines
        for i = 1:length(info.hmdef)
            if ishandle(info.hmdef(i)),	delete(info.hmdef(i)); end
        end
        % store info in figure
        guidata(h,info);
        return;
    end
else
    idx = true(size(tim));
end

% select data
C = [0 0 0; 0 0 1]; % 1:black and 2:blue
for i = 1:length(p)
    pname = info.param{p(i)};
    m = info.(pname).midx(info.m(ai(i)));
    a = info.(pname).aidx(info.a(ai(i)));
    dat = squeeze(info.(pname).dat{r}(m,a,idx));
    
    % plot data
    hold(info.ha(ai(i)),'on');
    if ishandle(info.hl(ai(i))),    delete(info.hl(ai(i))); end
    info.hl(ai(i)) = plot(info.ha(ai(i)),tim,dat,'LineWidth',2,'Color',C(ai(i),:));
    
    % adjust y limits
    set_ylim(info.ha(ai(i)),info.(pname),info.m(ai(i)),info.a(ai(i)));
    info.(pname).v{ai(i)} = axis;
    
    % set color and YAxisLocation
    set(info.ha(ai(i)),'Color','none','YColor',C(ai(i),:),'LineWidth',2);
    if ai(i) == 1
        set(info.ha(ai(i)),'YAxisLocation','left');
    else
        set(info.ha(ai(i)),'YAxisLocation','right','XTick',[]);
    end
    
end

% movement definitions
if info.plotmovdef
    
    % remove old movdef lines
    for i = 1:length(info.hmdef)
        if ishandle(info.hmdef(i)),	delete(info.hmdef(i)); end
    end
    
    % get ylim to rescale
    ylim = get(info.ha(1),'Ylim');
    ywin = diff(ylim);
    
    % plot in first axis
    hold(info.ha(1),'on');
   
    % individual parameters
    indiv = info.movdef.win{info.r}(idx,info.movdef.idx.indiv);
    % give a 1% offset from bottom and top axes
    indiv = indiv*ywin*0.99 + ywin*0.005 + ylim(1);
    % plot
    for i = 1:size(indiv,2)
        info.hmdef(i) = plot(info.ha(1),tim,indiv(:,i),'LineWidth',1,'LineStyle','--','Color',str2rgb(info.movdef.color{i}));
    end
    
    % combine parameters
    comb = prod(info.movdef.win{info.r}(idx,info.movdef.idx.comb),2);
    if ~isempty(info.movdef.idx.comb)
        if isempty(i),  i = 0;  end
        comb = comb*ywin*0.99 + ywin*0.005 + ylim(1);
        info.hmdef(i+1) = plot(info.ha(1),tim,comb,'LineWidth',2,'LineStyle','--','Color',str2rgb('r'));
    end
    
end

% plot trial markers
nvars = length(info.vars);
htl = nan(1,nvars);
htt = nan(1,nvars);
if nvars > 0
    v = axis;
    yline   = [v(3), v(3) + (v(4)-v(3))/20];
    ytext   = v(3) + (v(4)-v(3))/20;
end
for i = 1:nvars
    x = info.trl(r,i);
    htl(i) = line([x x],yline,'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2);
    htt(i) = text(x,ytext,info.vars{i},'HorizontalAlignment','center','Fontsize',8,'FontWeight','normal','Color',[0.5 0.5 0.5],'VerticalAlignment','bottom','interpreter','none');
end
info.htrl.line = htl;
info.htrl.text = htt;

% store counter of this plot
info.currr = r;
info.previ = [];

% store info in figure
guidata(h,info);


%% function draw_mov
%----------------------------------------
function draw_mov(h)
% make sure the handle is a figure handle
if ~ishandle(h),    return; end
while ~isequal(get(h,'parent'),0)
    h = get(h,'parent');
end
info = guidata(h);

% return quickly if possible
if info.plotmovdef, return; end

% make first axis the current axis
%axes(info.ha(1));

% get movement time on- and off-sets
allmov = info.allmov{info.r};
selmov = info.selmov{info.r};
check  = info.check{info.r};

% combine overlapping selections
allmov = time2logic(allmov,info.time{info.r});
allmov = logic2time(allmov,info.time{info.r});
selmov = time2logic(selmov,info.time{info.r});
selmov = logic2time(selmov,info.time{info.r});
info.allmov{info.r} = allmov;
info.selmov{info.r} = selmov;

% get movement counter
i = info.i;
if i > size(allmov,1); i = size(allmov,1);  end
if i < 1 && size(allmov,1)>0,   i = 1;      end
info.i  = i;
if i < 1,   i = 1;  end

% get patch handles
hma     = info.hm.all;
hms     = info.hm.sel;
hta     = info.hm.text;
htc     = info.hm.check;
   
% erase old patches
for jj = 1:length(hma)
    if ishandle(hma(jj)),	delete(hma(jj)); end
end
for jj = 1:length(hms)
    if ishandle(hms(jj)),	delete(hms(jj)); end
end
for jj = 1:length(hta)
    if ishandle(hta(jj)),	delete(hta(jj)); end
end
for jj = 1:length(htc)
    if ishandle(htc(jj)),	delete(htc(jj)); end
end
    
% adjust x limits
set_xlim(info);
v = axis;
yfill	= v([3 3 4 4]);
ytext   = v(4) - (v(4)-v(3))/10;
ycheck  = v(4) - (v(4)-v(3))/8;

% find movements in time window
if info.continuous && ~info.plotmovdef
    n = find(allmov(:,2)<v(2),1,'last');
    ii = find(selmov(:,1)>v(1),1,'first');
    nn = find(selmov(:,2)<v(2),1,'last');
else
    n = size(allmov,1);
    ii = 1;
    nn = size(selmov,1);
end
if isempty(n),  n  = 0;     end
if isempty(nn), nn = 0;    	end
if isempty(ii), ii = nn+1;	end
    
% plot all movements in the time window
%hold(info.ha(1),'on');
jj = 1;
hma = nan(1,n-i+1);
hta = nan(1,n-i+1);
for j = i:n
    hma(jj) = fill(allmov(j,[1 2 2 1]),yfill,'g',...
        'EdgeColor','g',...
        'LineStyle','--',...
        'LineWidth',2,...
        'FaceAlpha',0.2);
    
    xtext = mean(allmov(j,[1 2]));
    hta(jj) = text(xtext,ytext,num2str(j),'HorizontalAlignment','center','Fontsize',16,'FontWeight','normal','Color',[0.5 0.5 0.5],'interpreter','none');
    if info.showchecks
        checkstr = sprintf('%s\n',info.checkparam{logical(check(j,:))});
        htc(jj) = text(xtext,ycheck,checkstr,'HorizontalAlignment','center','Fontsize',10,'FontWeight','normal','Color',[0.5 0.5 0.5],'VerticalAlignment','top','interpreter','none');
    end
    jj = jj + 1;
end

% retrieve xpositions of the movements
if length(hta)>1
    xtext = cellfun(@(x) x(:,1),get(hta,'Position'));
end

% plot selected movements in the time window
jj = 1;
hms = nan(1,nn-ii+1);
for j = ii:nn
    hms(jj) = fill(selmov(j,[1 2 2 1]),yfill,'r',...
        'EdgeColor','r',...
        'LineStyle','--',...
        'LineWidth',2,...
        'FaceAlpha',0.2);
    jj = jj + 1;
    % find movements from allmov that overlap and emphasize them
    xmov = selmov(j,[1 2]);
    if ~isempty(hta)
        idxtext = xtext>=xmov(1) & xtext<=xmov(2);
        if any(idxtext), set(hta(idxtext),'Fontsize',18,'FontWeight','bold','Color','k');   end
    end
end
%hold(info.ha(1),'off');
    
% store counter of this plot
info.curri = i:n;
if isempty(info.curri)
    info.curri = 0;
end

% set correct linewidth
if ~isempty(hma) || ~isempty(hms)
    set(info.ha,'LineWidth',1.5);
    set(info.hl,'LineWidth',1.5);
else
    set(info.ha,'LineWidth',2);
    set(info.hl,'LineWidth',2);
end

% store patch handles
info.hm.all  = hma;
info.hm.sel	 = hms;
info.hm.text = hta;
info.hm.check= htc;

% store info in figure
guidata(h,info);

% define area boundaries
selinfo.tol = 0.005*(v(2)-v(1));    % accuracy: 0.5%
customarea = selmov(ii:nn,:)';
customarea = customarea(:);
customarea = [customarea-selinfo.tol customarea+selinfo.tol];
selinfo.pointer = {customarea,'custom'; [allmov(i:n,:); selmov(ii:nn,:)],'hand'; v(1:2),'crosshair'; [-Inf Inf],'arrow'};
% FIXME: implement y-restrictions in a more subtle way
selinfo.ylim = v(3:4);
selinfo.yfill = yfill;

% define selections
selinfo.sel = selmov(ii:nn,:);
selinfo.desel = allmov(i:n,:);
for j = 1:size(selinfo.sel,1);
    idx = selinfo.desel(:,1) >= selinfo.sel(j,1) & selinfo.desel(:,2) <= selinfo.sel(j,2);
    selinfo.desel = selinfo.desel(~idx,:);
end
selinfo.sel = [selinfo.sel zeros(size(selinfo.sel,1),1)];
selinfo.desel = [selinfo.desel zeros(size(selinfo.desel,1),1)];
selinfo.slidesel = zeros(0,2);
selinfo.hb.slide.slide = info.hb.slide.slide;
selinfo.idx = [ii nn];
selinfo.h = hms;

% store selection info in figure
setappdata(h,'select_m',selinfo);


%% function redraw
%----------------------------------------
function redraw(h)
% make sure the handle is a figure handle
if ~ishandle(h),    return; end
while ~isequal(get(h,'parent'),0)
    h = get(h,'parent');
end
info = guidata(h);

% enable or disable trial buttons
if info.r <= 1
    set(info.hb.tprev,'Enable','off');
else
    set(info.hb.tprev,'Enable','on');
end
if info.r >= info.nrun
    set(info.hb.tnext,'Enable','off');
else
    set(info.hb.tnext,'Enable','on');
end
if info.nrun > 10
    if info.r <= 1
        set(info.hb.tjprev,'Enable','off');
    else
        set(info.hb.tjprev,'Enable','on');
    end
    if info.r >= info.nrun
        set(info.hb.tjnext,'Enable','off');
    else
        set(info.hb.tjnext,'Enable','on');
    end
end

% enable or disable movement buttons
if info.continuous && ~info.plotmovdef
    if min(info.curri) <= 1 && info.r <= 1
        set(info.hb.mprev,'Enable','off');
    else
        set(info.hb.mprev,'Enable','on');
    end
    if max(info.curri) >= info.nallmov(info.r) && info.r >= info.nrun
        set(info.hb.mnext,'Enable','off');
    else
        set(info.hb.mnext,'Enable','on');
    end
    if any(info.nallmov > 10)
        if min(info.curri) <= 1 && info.r <= 1
            set(info.hb.mjprev,'Enable','off');
        else
            set(info.hb.mjprev,'Enable','on');
        end
        if max(info.curri) >= info.nallmov(info.r) && info.r >= info.nrun
            set(info.hb.mjnext,'Enable','off');
        else
            set(info.hb.mjnext,'Enable','on');
        end
    end
end

% update checks button and x-slider
if info.plotmovdef
    set_xslider(info);
elseif info.showchecks
    set(info.hb.checks,'String','no checks');
else
    set(info.hb.checks,'String','checks');
end

% enable or disable drop down lists
for i = 1:length(info.p)
    pname = info.param{info.p(i)};
    if length(info.(pname).marker) > 1
        set(info.hb.marker(i),'Enable','on');
    else
        set(info.hb.marker(i),'Enable','off');
    end
    if length(info.(pname).axis) > 1
        set(info.hb.axis(i),'Enable','on');
    else
        set(info.hb.axis(i),'Enable','off');
    end
end
if info.nparam > 1
    set(info.hb.param(1),'Enable','on');
else
    set(info.hb.param(1),'Enable','off');
end

if ~info.plotmovdef
    % enable or disable the slide button
    selinfo = getappdata(h,'select_m');
    if isfield(selinfo,'slidesel') && ~isempty(selinfo.slidesel)
        set(info.hb.slide.slide,'Enable','on');
    else
        set(info.hb.slide.slide,'Enable','off');
    end
    % enable or disable slide drop down lists
    list = get(info.hb.slide.param,'String');
    idx = get(info.hb.slide.param,'val');
    pname = list{idx};
    if length(info.(pname).marker) > 1
        set(info.hb.slide.marker,'Enable','on');
    else
        set(info.hb.slide.marker,'Enable','off');
    end
    if length(info.(pname).axis) > 1
        set(info.hb.slide.axis,'Enable','on');
    else
        set(info.hb.slide.axis,'Enable','off');
    end
end

% update trial and movement counters in the edit windows
% run/trial number
set(info.hb.tedit,'String',num2str(info.r));
% movement number
if info.continuous && ~info.plotmovdef
    if ~isempty(info.i) && info.i>0
        set(info.hb.medit,'String',num2str(info.i));
    else
        set(info.hb.medit,'String','n/a');
    end
end
% xlim
set(info.hb.xlim,'String',get_limstr(info.xlim));
% ylim
for i = 1:length(info.p)
    pname = info.param{info.p(i)};
    %ylimstr = get_limstr(info.(pname).ylim);
    v = axis(info.ha(i));
    ylimstr = sprintf('[%g %g]',v(3:4));
    set(info.hb.ylim(i),'String',ylimstr);
end
%--------------------------------------------------------------------------


%% function pointer_lrdrag
%----------------------------------------
function out = pointer_lrdrag(flg)
if nargin < 1,  flg = 'CData';  end
switch lower(flg)
    case 'cdata'
        o=NaN; w=2; k=1;
        out = [...
            o o o o o w w w w w o o o o o o
            o o o o o w k w k w o o o o o o
            o o o o o w k w k w o o o o o o
            o o o o w w k w k w w o o o o o
            o o o w k w k w k w k w o o o o
            o o w k k w k w k w k k w o o o
            o w k k k k k w k k k k k w o o
            w k k k k k k w k k k k k k w o
            o w k k k k k w k k k k k w o o
            o o w k k w k w k w k k w o o o
            o o o w k w k w k w k w o o o o
            o o o o w w k w k w w o o o o o
            o o o o o w k w k w o o o o o o
            o o o o o w k w k w o o o o o o
            o o o o o w w w w w o o o o o o
            o o o o o o o o o o o o o o o o];
    case 'hotspot'
        out = [8 8];
    otherwise
end

%% function m_select
%----------------------------------------
function varargout = m_select(h, eventdata, varargin)

% get the optional arguments
event    = keyval('event',    varargin);
callback = keyval('callback', varargin);

% make sure the handle is a figure handle
while ~isequal(get(h,'parent'),0)
    h = get(h,'parent');
end

% get user data from figure
if ishandle(h)
    info = getappdata(h,'select_m');
else
    info = [];
end

% initialize user specified info
if isempty(info)
	% maybe I could do something here...
    set(h, 'Pointer','arrow');
    return
end

% get pointer location
p = get(gca, 'CurrentPoint');
p = p(1,1:2);
px = p(1); py = p(2);

% return if outside of y limits
if py < min(info.ylim) || py > max(info.ylim)
    set(h, 'Pointer','arrow');
    return
end

% loop over pointer areas to change cursor
if all(info.sel(:,3)==0) && all(info.desel(:,3)==0)
    for i = 1:size(info.pointer,1)
        if any(px >= info.pointer{i,1}(:,1) & px <= info.pointer{i,1}(:,2))
            set(h, 'Pointer',info.pointer{i,2});
            break
        end
    end
end

% control button clicks
switch event
    case 'ButtonDown'
        idxleft = px >= info.sel(:,1)-info.tol & px <= info.sel(:,1)+info.tol;
        idxright = px >= info.sel(:,2)-info.tol & px <= info.sel(:,2)+info.tol;
        idxsel = px >= info.sel(:,1) & px <= info.sel(:,2);
        idxdesel = px >= info.desel(:,1) & px <= info.desel(:,2);
        if any(idxleft)
            info.sel(:,3)	= 0;
            info.sel(find(idxleft,1,'first'),3) = 10;
        elseif any(idxright)
            info.sel(:,3)	= 0;
            info.sel(find(idxright,1,'first'),3) = 20;
        elseif any(idxsel)
            info.desel(:,3) = 0;
            info.sel(:,3)	= 0;
            info.sel(find(idxsel,1,'first'),3) = 1;
        elseif any(idxdesel)
            info.desel(:,3) = 0;
            info.sel(:,3)	= 0;
            info.desel(find(idxdesel,1,'first'),3) = 1;
        else
            sel = [px px 20];
            hms = fill(sel([1 2 2 1]),info.yfill,[0.5 0.5 0.5],...
                'FaceColor','none',...
                'EdgeColor',[0.5 0.5 0.5],...
                'LineStyle','--',...
                'LineWidth',1.5,...
                'FaceAlpha',0.2);
            info.h(end+1) = hms;
            info.sel(end+1,:) = sel;
            % update pointer areas
            idx = ismember(info.pointer(:,2),'hand');
            info.pointer{idx,1} = [info.desel(:,1:2); info.sel(:,1:2)];
        end
        setappdata(h,'select_m',info);
    case 'ButtonUp'
        idxleft = info.sel(:,3)==10;
        idxright = info.sel(:,3)==20;
        idxsel = px >= info.sel(:,1) & px <= info.sel(:,2) & info.sel(:,3)==1;
        idxdesel = px >= info.desel(:,1) & px <= info.desel(:,2) & info.desel(:,3)==1;
        % invoke callback function
        if any(idxleft) || any(idxright)
            evalCallback(callback,info.sel);
            if any(idxleft)
                slidesel = [find(idxleft,1,'first') 1];
            elseif any(idxright)
                slidesel = [find(idxright,1,'first') 2];
            end
            % enable the slide button
            set(info.hb.slide.slide,'Enable','on');
        elseif any(idxsel)
            sel = info.sel(~idxsel,:);
            evalCallback(callback,sel);
            slidesel = zeros(0,2);
            % disable the slide button
            set(info.hb.slide.slide,'Enable','off');
        elseif any(idxdesel)
            sel = [info.sel; info.desel(idxdesel,:)];
            sel = sortrows(sel,1);
            evalCallback(callback,sel);
            slidesel = zeros(0,2);
            % disable the slide button
            set(info.hb.slide.slide,'Enable','off');
        else
            slidesel = zeros(0,2);
        end
        info = getappdata(h,'select_m');
        info.slidesel = slidesel;
        setappdata(h,'select_m',info);
    case 'Motion'
        idxleft = info.sel(:,3)==10;
        idxright = info.sel(:,3)==20;
        if any(idxleft) && ishandle(info.h(idxleft))
            xdat = get(info.h(idxleft),'XData');
            xdat([1 4]) = px;
            set(info.h(idxleft),'XData',xdat);
            %info.sel(idxleft,1:2) = sort([px info.sel(idxleft,2)]);
            info.sel(idxleft,1:2) = sort([px xdat(2)]);
        elseif any(idxright) && ishandle(info.h(idxright))
            xdat = get(info.h(idxright),'XData');
            xdat([2 3]) = px;
            set(info.h(idxright),'XData',xdat);
            %info.sel(idxright,1:2) = sort([info.sel(idxright,1) px]);
            info.sel(idxright,1:2) = sort([xdat(1) px]);
        end
        if any(idxleft) || any(idxright)
            % update pointer areas
            customarea = info.sel(:,2)';
            idx = ismember(info.pointer(:,2),'custom');
            info.pointer{idx,1} = [customarea(:)-info.tol customarea(:)+info.tol];
            idx = ismember(info.pointer(:,2),'hand');
            info.pointer{idx,1} = [info.desel(:,1:2); info.sel(:,1:2)];
            setappdata(h,'select_m',info);
        end
    otherwise
end

%% function evalCallback
%----------------------------------------
function evalCallback(callback, val)
if ~isempty(callback)
    if isa(callback, 'cell')
        % the callback specifies a function and additional arguments
        funhandle = callback{1};
        funargs   = callback(2:end);
        if nargin < 2
            feval(funhandle, funargs{:});
        else
            feval(funhandle, val, funargs{:});
        end
    else
        % the callback only specifies a function
        funhandle = callback;
        if nargin < 2
            feval(funhandle);
        else
            feval(funhandle, val);
        end
    end
end

%% function m_select_cb
%----------------------------------------
function varargout = m_select_cb(sel,h)
info = guidata(h);
selinfo = getappdata(h,'select_m');
selmov = info.selmov{info.r};
idx = [selinfo.idx(1)-1 selinfo.idx(2)+1];
selmovnew = [selmov(1:idx(1),:); sel(:,1:2); selmov(idx(2):end,:)];
info.selmov{info.r} = selmovnew;
info.hm.sel = selinfo.h;
guidata(h,info);
draw_mov(h);

%% function b_quit
%----------------------------------------
function varargout = b_quit(h, eventdata)
info = guidata(h);
info.quit = true;
guidata(h,info);
uiresume;

%% function b_tprev: TRIAL JUMP
%----------------------------------------
function varargout = b_tjprev(h, eventdata)
info = guidata(h);
if info.r > 1
    info.r = max(info.r - 10,1);
    info.i = 1; % reset movement counter
    info.xlim = reset_xlim(info);  % reset xlim
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function b_tprev: TRIAL
%----------------------------------------
function varargout = b_tprev(h, eventdata)
info = guidata(h);
if info.r > 1
    info.r = info.r - 1;
    info.i = 1; % reset movement counter
    info.xlim = reset_xlim(info);  % reset xlim
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function e_tedit: TRIAL EDIT
%----------------------------------------
function varargout = e_tedit(h, eventdata)
info = guidata(h);
r = get(info.hb.tedit,'String');
try
    set(info.hb.tedit,'ForegroundColor','k');
    r = str2num(r);
    if isempty(r),      error(' ');     end
    if r < 1,           r = 1;          end
    if r > info.nrun,	r = info.nrun;   end
    if r ~= info.r
        info.r = r;
        info.i = 1; % reset movement counter
        info.xlim = reset_xlim(info);  % reset xlim
        guidata(h,info);
        draw_dat(h);
        draw_mov(h);
    end
    set(info.hb.tedit,'String',num2str(info.r));
catch
    set(info.hb.tedit,'String','err');
    set(info.hb.tedit,'ForegroundColor','r');
end
uiresume;

%% function b_tnext: TRIAL
%----------------------------------------
function varargout = b_tnext(h, eventdata)
info = guidata(h);
if info.r < info.nrun
    info.r = info.r + 1;
    info.i = 1; % reset movement counter
    info.xlim = reset_xlim(info);  % reset xlim
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function b_tnext: TRIAL JUMP
%----------------------------------------
function varargout = b_tjnext(h, eventdata)
info = guidata(h);
if info.r < info.nrun
    info.r = min(info.r + 10,info.nrun);
    info.i = 1; % reset movement counter
    info.xlim = reset_xlim(info);  % reset xlim
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function b_mjprev: MOV JUMP
%----------------------------------------
function varargout = b_mjprev(h, eventdata)
info = guidata(h);
if info.i > 1
    info.previ = [];
    info.i = max(info.i - 10,1);
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
elseif info.r > 1
    info.r = info.r - 1;
    info.i = info.nallmov(info.r); % reset movement counter
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function b_mprev: MOV
%----------------------------------------
function varargout = b_mprev(h, eventdata)
info = guidata(h);
if info.i > 1
    if ~isempty(info.previ) && info.previ < info.i
        info.i = info.previ;
    else
        info.i = info.i - 1;
    end
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
elseif info.r > 1
    info.r = info.r - 1;
    info.i = info.nallmov(info.r); % reset movement counter
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function e_medit: MOV EDIT
%----------------------------------------
function varargout = e_medit(h, eventdata)
info = guidata(h);
i = get(info.hb.medit,'String');
try
    set(info.hb.medit,'ForegroundColor','k');
    i = str2num(i);
    if isempty(i),                  error(' ');                 end
    if i < 1,                       i = 1;                      end
    if i > info.nallmov(info.r),	i = info.nallmov(info.r);	end
    if i ~= info.i
        info.previ = [];
        info.i = i;
        guidata(h,info);
        draw_dat(h);
        draw_mov(h);
    end
    set(info.hb.medit,'String',num2str(info.i));
catch
    set(info.hb.medit,'String','err');
    set(info.hb.medit,'ForegroundColor','r');
end
uiresume;

%% function b_mnext: MOV
%----------------------------------------
function varargout = b_mnext(h, eventdata)
info = guidata(h);
if info.i < info.nallmov(info.r)
    info.previ = info.i;
    info.i = max(info.curri) + 1;
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
elseif info.r < info.nrun
    info.r = info.r + 1;
    info.i = 1; % reset movement counter
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function b_mnext: MOV JUMP
%----------------------------------------
function varargout = b_mjnext(h, eventdata)
info = guidata(h);
if info.i < info.nallmov(info.r)
    info.previ = [];
    info.i = min(max(info.curri) + 10,info.nallmov(info.r));
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
elseif info.r < info.nrun
    info.r = info.r + 1;
    info.i = 1; % reset movement counter
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
end
uiresume;

%% function p_absreltim: relative or absolute time window
%----------------------------------------
function varargout = p_absreltim(h, eventdata)
info = guidata(h);
str = get(info.hb.absreltim,'String');
idx = get(info.hb.absreltim,'Value');
info.absreltim = lower(str{idx});
set_xlim(info);
guidata(h,info);
%redraw(h);
    
%% function b_xprev: XLIM
%----------------------------------------
function varargout = b_xprev(h, eventdata)
info = guidata(h);
info.xlim = info.xlim - diff(info.xlim);
% update figure
set_xlim(info);
set_xslider(info);
guidata(h,info);
draw_dat(h);
draw_mov(h);
%redraw(h);
uiresume;

%% function b_xlim: change xlim (time window)
%----------------------------------------
function varargout = b_xlim(h, eventdata)
info = guidata(h);
timwin = get(info.hb.xlim,'String');
try
    % reset edit window text color to default
    set(info.hb.xlim,'ForegroundColor','k');
    
    % string to number if requiered
    if ~isequal(timwin,'maxmin')
        try     timwin = str2num(timwin);
        catch,  timwin = eval(timwin);  end
        if ~isnumeric(timwin),  error(' '); end
    end
    
    % return if nothing to do
    if isequal(info.xlim,timwin)
        return
    end
    
    % update configuration in info
    info.xlim = timwin;
    
    % update figure
    set_xlim(info);
    set_xslider(info);
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
    %redraw(h);
catch
    %set(info.hb.xlim,'String','error');
    set(info.hb.xlim,'ForegroundColor','r');
end
uiresume;
    
%% function b_xnext: XLIM
%----------------------------------------
function varargout = b_xnext(h, eventdata)
info = guidata(h);
info.xlim = info.xlim + diff(info.xlim);
% update figure
set_xlim(info);
set_xslider(info);
guidata(h,info);
draw_dat(h);
draw_mov(h);
%redraw(h);
uiresume;

%% function b_ylim: change ylim
%----------------------------------------
function varargout = b_ylim(h, eventdata, ai)
info = guidata(h);
datlim = get(info.hb.ylim(ai),'String');
try
    % reset edit window text color to default
    set(info.hb.ylim(ai),'ForegroundColor','k');
    
    % string to number if requiered
    if ~isequal(datlim,'maxmin') && ~isequal(datlim,'auto')
        try     datlim = str2num(datlim);
        catch,  datlim = eval(datlim);  end
        if ~isnumeric(datlim),  error(' '); end
    end
    
    % return if nothing to do
    pname = info.param{info.p(ai)};
    if isequal(info.(pname).ylim,datlim)
        return
    end
    
    % update configuration in info
    info.(pname).ylim = datlim;
    
    % update figure
    set_ylim(info.ha(ai),info.(pname),info.m(ai),info.a(ai));
    info.(pname).v{ai} = axis;
    guidata(h,info);
    draw_mov(h);
    %redraw(h);
catch
    %set(info.hb.ylim,'String','error');
    set(info.hb.ylim(ai),'ForegroundColor','r');
end
uiresume;

%% function s_xslider: change xlim by xslider
%----------------------------------------
function varargout = s_xslider(h, eventdata)
info = guidata(h);
val = get(info.hb.xslider,'Value');
try
    % reset edit window and slider color to default
    set(info.hb.xlim,'ForegroundColor','k');
    set(info.hb.xslider,'BackgroundColor','w');
    
    % xlim = set_xlim(info);
    xlim = get(info.ha(1),'Xlim');
    xwin = diff(xlim);
    timwin = [val-(xwin/2) val+(xwin/2)];
    
    % return if nothing to do
    if isequal(info.xlim,timwin)
        return
    end
    
    % update configuration in info
    info.xlim = timwin;
    
    % update figure
    set_xlim(info);
    %set_xslider(info);
    guidata(h,info);
    draw_dat(h);
    draw_mov(h);
    %redraw(h);
    
catch
    %set(info.hb.xlim,'String','error');
    set(info.hb.xlim,'ForegroundColor','r');
    set(info.hb.xslider,'BackgroundColor','r');
end
uiresume;

%% function b_param: change param
%----------------------------------------
function varargout = p_param(h, eventdata, ai)
info = guidata(h);
%str = get(info.hb.param(ai),'String');
info.p(ai) = get(info.hb.param(ai),'Value');
pname = info.param{info.p(ai)};
set(info.hb.marker(ai),'String',info.(pname).marker,'Value',1);
set(info.hb.axis(ai),'String',info.(pname).axis,'Value',1);
info.m(ai) = 1;
info.a(ai) = 1;
guidata(h,info);
draw_dat(h,ai);
draw_mov(h);
%redraw(h);

%% function b_marker: change marker
%----------------------------------------
function varargout = p_marker(h, eventdata, ai)
info = guidata(h);
pname = info.param{info.p(ai)};
info.m(ai) = info.(pname).midx(get(info.hb.marker(ai),'Value'));
guidata(h,info);
draw_dat(h,ai);
draw_mov(h);
%redraw(h);

%% function b_param: change axis
%----------------------------------------
function varargout = p_axis(h, eventdata, ai)
info = guidata(h);
pname = info.param{info.p(ai)};
info.a(ai) = info.(pname).aidx(get(info.hb.axis(ai),'Value'));
guidata(h,info);
draw_dat(h,ai);
draw_mov(h);
%redraw(h);

%% function b_checks: show checks
%----------------------------------------
function varargout = b_checks(h, eventdata)
info = guidata(h);
info.showchecks = get(info.hb.checks,'Value')==get(info.hb.checks,'Max');
guidata(h,info);
draw_mov(h);
%redraw(h);

%% function b_slide
%----------------------------------------
function varargout = b_slide(h, eventdata)
while ~isequal(get(h,'parent'),0)
    h = get(h,'parent');
end
selinfo = getappdata(h,'select_m');
slidesel = selinfo.slidesel;
if isempty(slidesel),   return; end;
info = guidata(h);
% get slide data
list = get(info.hb.slide.param,'String');
idx = get(info.hb.slide.param,'Value');
pname = list{idx};
m = get(info.hb.slide.marker,'Value');
a = get(info.hb.slide.axis,'Value');
% slide!
dat = squeeze(info.(pname).dat{info.r}(m,a,:));
selidx = selinfo.idx(1) - 1 + slidesel(1,1);
sel = time2logic(info.selmov{info.r}(selidx,:),info.time{info.r})';
sel = km_slidesel(dat,sel,info.(pname).slide);
sel = logic2time(sel,info.time{info.r});
info.selmov{info.r}(selidx,slidesel(1,2)) = sel(1,slidesel(1,2));
% update figure
guidata(h,info);
draw_mov(h);
redraw(h);
selinfo = getappdata(h,'select_m');
selinfo.slidesel = slidesel;
setappdata(h,'select_m',selinfo);

%% function b_slideparam: change slide parameter
%----------------------------------------
function varargout = p_slideparam(h, eventdata)
info = guidata(h);
list = get(info.hb.slide.param,'String');
idx = get(info.hb.slide.param,'Value');
pname = list{idx};
set(info.hb.slide.marker,'String',info.(pname).marker,'Value',1);
set(info.hb.slide.axis,'String',info.(pname).axis,'Value',1);
set(info.hb.slide.crit,'String',info.(pname).slide);
guidata(h,info);
%redraw(h);

%% function b_slidecrit: slide criterium
%----------------------------------------
function varargout = b_slidecrit(h, eventdata)
info = guidata(h);
crit = get(info.hb.slide.crit,'String');
try
    % reset edit window text color to default
    set(info.hb.slide.crit,'ForegroundColor','k');
    
    % string to number if requiered
    try     crit = str2num(crit);
    catch,  crit = eval(crit);  end
    if ~isnumeric(crit),  error(' '); end
    
    % return if nothing to do
    list = get(info.hb.slide.param,'String');
    idx = get(info.hb.slide.param,'Value');
    pname = list{idx};
    if isequal(info.(pname).slide,crit)
        return
    end
    
    % update configuration in info
    info.(pname).slide = crit;
    
    % update figure
    guidata(h,info);
    %redraw(h);
catch
    %set(info.hb.ylim,'String','error');
    set(info.hb.ylim(ai),'ForegroundColor','r');
end
uiresume;

%% function b_movdef: plot movement definitions
%----------------------------------------
function varargout = b_movdef(h, eventdata)
info = guidata(h);
try
    % reset edit window text color to default
    set(info.hb.movdef.plot,'ForegroundColor','k');
    
    % get parameter indeces
    info.movdef.idx.indiv = get(info.hb.movdef.indiv,'Value');
    info.movdef.idx.comb = get(info.hb.movdef.comb,'Value');
    
    % update figure
    guidata(h,info);
    draw_dat(h);
    %redraw(h);
catch
    set(info.hb.movdef.plot,'ForegroundColor','r');
end
uiresume;


%% function set_xlim
%----------------------------------------
function tmpxlim = set_xlim(info)
allmov = info.allmov{info.r};
i = info.i;
if i > size(allmov,1); i = size(allmov,1);  end

if isequal(info.xlim,'maxmin')
    for ai = 1:length(info.ha)
        xlim(info.ha(ai),'auto');
    end
    v = axis(info.ha(1));
    tmpxlim = v(1:2);
    return
else
    if strcmpi(info.absreltim,'abs')
        tmpxlim = info.xlim;
    elseif strcmpi(info.absreltim,'rel') && i > 0 && ~isempty(info.allmov{info.r})
        tmpxlim = allmov(i,1) + info.xlim;
        tmpxlim = [min([allmov(i,1) tmpxlim(1)]) max([allmov(i,2) tmpxlim(2)])];
        tmpxlim(1) = floor(tmpxlim(1)/0.01)*0.01;
        tmpxlim(2) = ceil(tmpxlim(2)/0.01)*0.01;
    else
        tmpxlim = [-Inf Inf];
    end
    if any(~isfinite(tmpxlim))
        xlim(info.ha(1),'auto');
        v = axis(info.ha(1));
        if ~isfinite(tmpxlim(1)), tmpxlim(1) = v(1);    end
        if ~isfinite(tmpxlim(2)), tmpxlim(2) = v(2);    end
    end
    info.xlim = tmpxlim;
    for ai = 1:length(info.ha)
        xlim(info.ha(ai),tmpxlim);
    end
end

% store info in figure if applicable
if isfield(info,'hf') && ishandle(info.hf)
    guidata(info.hf,info);
end

%% function set_ylim
%----------------------------------------
function tmpylim = set_ylim(ha,info,m,a)
% determine y range
if isequal(info.ylim,'auto')
    nrun = length(info.dat);
    mcell = repmat({info.midx(m)},1,nrun);
    acell = repmat({info.aidx(a)},1,nrun);
    ymin = cellfun(@(x,i,j) min(squeeze(x(i,j,:))),info.dat,mcell,acell);
    ymin = min(ymin(~cutoutlier(ymin,1.5,'iqr','lower')));
    ymax = cellfun(@(x,i,j) max(squeeze(x(i,j,:))),info.dat,mcell,acell);
    ymax = max(ymax(~cutoutlier(ymax,1.5,'iqr','higher')));
    sc = 10.^(floor(log10(ymax-ymin))-1);
    ymin = floor(ymin/sc)*sc;
    ymax = ceil(ymax/sc)*sc;
    tmpylim = [ymin ymax];
else
    tmpylim = info.ylim;
end

% set y limits if applicable
if ishandle(ha)
    if isequal(info.ylim,'maxmin')
        ylim(ha,'auto');
        v = axis;
        tmpylim = v(3:4);
    else
        ylim(ha,tmpylim);
    end
end

%% function set_xslider
%----------------------------------------
function set_xslider(info,xlim)
if ~info.plotmovdef || ~isfield(info.hb,'xslider') || ~ishandle(info.hb.xslider)
    return;
end
if nargin < 2
    %xlim = info.xlim;
    xlim = get(info.ha(1),'Xlim');
end
xwin = diff(xlim);
tlim = [min(info.time{info.r}) max(info.time{info.r})];
twin = diff(tlim);
set(info.hb.xslider,'Min',tlim(1),'Max',tlim(2),'Value',mean(xlim),'SliderStep',[xwin/(5*twin) xwin/twin]);

%% function reset_xlim
%----------------------------------------
function xlim = reset_xlim(info,xlim)
if ~strcmpi(info.absreltim,'abs')
    xlim = info.xlim;
    return
end
if nargin < 2
    %xlim = info.xlim;
    xlim = get(info.ha(1),'Xlim');
end
xwin = diff(xlim);
tlim = [min(info.time{info.r}) max(info.time{info.r})];
xlim = [tlim(1) tlim(1)+xwin];
if xlim(2) > tlim(2),	xlim(2) = tlim(2);	end

%% function get_limstr
%----------------------------------------
function limstr = get_limstr(lim)
% transform a limint specification into a printable string
if isequal(lim,'maxmin') || isequal(lim,'auto')
    limstr = lim;
else
    if ischar(lim)
        try     lim = str2num(lim);
        catch,  lim = eval(lim);	end
    end
    limstr = sprintf('[%g %g]',lim);
end