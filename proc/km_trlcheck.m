function varargout = km_trlcheck(cfg,data)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('KM:Return','%s: Nothing to do...',task);
    return
end

% get movement
if isfield(cfg,'movement') && ~isstruct(cfg.movement) && ~isempty(cfg.movement)
    movement = cfg.movement;
elseif isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
    movement = cfg.movdef.movement;
elseif isfield(data,'movement')
    movement = data.movement;
else
    error('no movement provided');
end
oldmovement = movement;


% loop over trials
ntrl = length(data.time);
nparam = length(cfg.trlcheck.param);
check = false(ntrl,nparam);
for t = 1:ntrl
    
    % get movement time on- and off-sets
    mov = movement{t};
    
    % remove all NaN movement on- and off-sets
    mov = mov(~any(isnan(mov),2),:);
    
    % continue immediately if no movements are present
    if isempty(mov),   continue;   end
    
    % loop over parameters
    for p = 1:nparam
        
        % parameter name
        pname = lower(cfg.trlcheck.param{p});
        tcfg = cfg.trlcheck.(pname);
        
        % switch between check parameters
        switch pname
                
            case 'bvel'
                
                % select data
                markersel	= km_labelselection(tcfg.marker,data.label);
                markersel	= ismember(data.label,markersel);
                axissel     = km_labelselection(tcfg.axis,data.derivaxis);
                axissel     = ismember(data.derivaxis,axissel);
                vel = data.vel{t}(markersel,axissel,:);
                
                % force to absolute values
                if ~isfalse(tcfg.abs),  vel = abs(vel); end
                
                % begin velocity values
                bvel = squeeze(vel(:,:,1));
                
                % remove movements that do not reach the criterium
                if length(tcfg.crit) == 1
                    idx = all(bvel > tcfg.crit);
                else
                    idx = all(bvel > min(tcfg.crit) & bvel < max(tcfg.crit));
                end
                
            case {'bga','ega'}
                
                % select data
                markersel	= km_labelselection(tcfg.marker,data.grip.label);
                markersel	= ismember(data.grip.label,markersel);
                axissel     = km_labelselection(tcfg.axis,data.grip.appaxis);
                axissel     = ismember(data.grip.appaxis,axissel);
                app = data.grip.app{t}(markersel,axissel,:);
                
                % determine begin/end grip apperture
                if strcmpi(pname,'bga')
                    ga = squeeze(app(:,:,1));
                elseif strcmpi(pname,'ega')
                    ga = squeeze(app(:,:,end));
                end
                
                % remove movements that do not reach the criterium
                if length(tcfg.crit) == 1
                    idx = all(ga > tcfg.crit);
                else
                    idx = all(ga > min(tcfg.crit) & ega < max(tcfg.crit));
                end
                
            case {'pos','epos','bpos'}
                
                % select markers
                markersel	= km_labelselection(tcfg.marker,data.label);
                markersel	= ismember(data.label,markersel);
                if ~isfalse(tcfg.refmarker)
                    refsel  = km_labelselection(tcfg.refmarker,data.label);
                    refsel	= ismember(data.label,refsel);
                else
                    refsel = [];
                end
                
                % simple (single axis) criterium
                if ~ischar(tcfg.crit)
                    % get data
                    axissel = km_labelselection(tcfg.axis,data.axis);
                    axissel = ismember(data.axis,axissel);
                    pos = squeeze(data.pos{t}(markersel,axissel,:));
                    if ~isempty(refsel)
                        refpos = squeeze(data.pos{t}(refsel,axissel,:));
                        pos = pos - refpos;
                    end
                    
                    % apply criteria
                    if strcmpi(pname,'bpos')
                        tmp = pos(1);
                    elseif strcmpi(pname,'epos')
                        tmp = pos(end);
                    elseif strcmpi(pname,'pos')
                        tmp = pos;
                    end
                    if length(tcfg.crit) == 1
                        idx = any(tmp > tcfg.crit);
                    else
                        idx = any(tmp > min(tcfg.crit) & tmp < max(tcfg.crit));
                    end
                    
                % complex (or multiple axes) criterium
                else
                    % get data
                    pos = data.pos{t}(markersel,:,:);
                    if ~isempty(refsel)
                        refpos = data.pos{t}(refsel,:,:);
                        pos = pos - refpos;
                    end
                    
                    % get axis data
                    ax = regexp(tcfg.crit,'(?<![a-zA-Z])[a-zA-Z](?![a-zA-Z])','match');
                    ax = unique(ax);
                    for a = 1:length(ax)
                        axissel = km_labelselection(ax{a},data.axis);
                        axissel = ismember(data.axis,axissel);
                        axpos.(ax{a}) = squeeze(pos(:,axissel,:));
                    end
                    
                    % apply criterium
                    for a = 1:length(ax)
                        if strcmpi(pname,'bpos')
                            tmp = axpos.(ax{a})(1);
                        elseif strcmpi(pname,'epos')
                            tmp = axpos.(ax{a})(end);
                        elseif strcmpi(pname,'pos')
                            tmp = axpos.(ax{a});
                        end
                        eval(sprintf('%s = tmp;',ax{a}));
                    end
                    
                    % evaluate cirterium
                    idx = any(eval(tcfg.crit));
                    
                end
                
                % apply selection
                idx = logical(idx);
                if strcmpi(tcfg.type,'excl')
                    idx = ~idx;
                end
                
            otherwise
                error('parameter ''%s'' not recognized',pname);
                
        end
        
        % store selection
        check(t,p) = ~idx;
        
    end
    
    % remove movements definition from selected trials
    if any(check(t,:),2);
        movement{t} = nan(0,2);
    end
        
end
cfg.trlcheck.check = check;
cfg.trlcheck.allmov = oldmovement;

% return movement definition or store in data and configuration structure
if nargout == 1
    varargout = {movement};
else
    if length(cfg.movdef) == 1
        cfg.movdef.movement = movement;
    end
    cfg.movement	= movement;
    data.movement	= movement;
    varargout       = {cfg,data};
end
