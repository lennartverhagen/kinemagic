function varargout = km_movcheck(cfg,data)
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

% loop over runs
check = cell(1,length(movement));
for r = 1:length(data.time)
    
    % get movement time on- and off-sets
    mov = movement{r};
    
    % remove all NaN movement on- and off-sets
    mov = mov(~any(isnan(mov),2),:);
    
    % continue immediately if no movements are present
    if isempty(mov),   continue;   end
    
    % loop over movcheck
    check{r} = false(size(mov,1),length(cfg.movcheck(1).param),length(cfg.movcheck));
    for c = 1:length(cfg.movcheck)
        
        % loop over parameters
        for p = 1:length(cfg.movcheck(c).param)
            
            % parameter name
            pname = lower(cfg.movcheck(c).param{p});
            tcfg = cfg.movcheck(c).(pname);
            
            % switch between check parameters
            switch pname
                
                case 'dur'
                    
                    % define periods without movement
                    nomov = [mov(1:end-1,2) mov(2:end,1)];
                    nomov = diff(nomov,1,2);
                    
                    % identify the movement on- and off-sets to keep
                    idx_on = [inf; nomov] >= tcfg.off;
                    idx_off = [nomov; inf] >= tcfg.off;
                    
                    % combine movements with too short breaks
                    mov = [mov(idx_on,1) mov(idx_off,2)];
                    
                    % now check if the movement that remain are of an
                    % appropriate duration
                    movdur = diff(mov,1,2);
                    if length(tcfg.on) == 1
                        idx = movdur > tcfg.on;
                    else
                        idx = movdur > min(tcfg.on) & movdur < max(tcfg.on);
                    end
                    
                case 'peakvel'
                    
                    % select data
                    markersel	= km_labelselection(tcfg.marker,data.label);
                    markersel	= ismember(data.label,markersel);
                    axissel     = km_labelselection(tcfg.axis,data.derivaxis);
                    axissel     = ismember(data.derivaxis,axissel);
                    vel = data.vel{r}(markersel,axissel,:);
                    
                    % force to absolute values
                    if ~isfalse(tcfg.abs),  vel = abs(vel); end
                    
                    % go from the time to the sample domain
                    idx = time2idx(mov,data.time{r});
                    
                    % determine peak velocity
                    peakvel = nan(size(idx,1),1);
                    for i = 1:size(idx,1)
                        peakvel(i) = max(max(max(vel(:,:,idx(i,1):idx(i,2)))));
                    end
                    
                    % remove movements that do not reach the criterium
                    if length(tcfg.crit) == 1
                        idx = peakvel > tcfg.crit;
                    else
                        idx = peakvel > min(tcfg.crit) & peakvel < max(tcfg.crit);
                    end
                    
                case {'bga','ega'}
                    
                    % select data
                    markersel	= km_labelselection(tcfg.marker,data.grip.label);
                    markersel	= ismember(data.grip.label,markersel);
                    axissel     = km_labelselection(tcfg.axis,data.grip.aptaxis);
                    axissel     = ismember(data.grip.aptaxis,axissel);
                    apt = data.grip.apt{r}(markersel,axissel,:);
                    
                    % go from the time to the sample domain
                    idx = time2idx(mov,data.time{r});
                    
                    % determine end grip apperture
                    if strcmp(pname,'bga')
                        ga = squeeze(apt(:,:,idx(:,1)));
                    elseif strcmp(pname,'ega')
                        ga = squeeze(apt(:,:,idx(:,2)));
                    end
                    
                    % remove movements that do not reach the criterium
                    if length(tcfg.crit) == 1
                        idx = ga > tcfg.crit;
                    else
                        idx = ga > min(tcfg.crit) & ga < max(tcfg.crit);
                    end
                    
                case {'pos','epos','bpos'}
                    
                    % select markers
                    markersel	= km_labelselection(tcfg.marker,data.label);
                    markersel	= ismember(data.label,markersel);
                    if ~isfalse(tcfg.refmarker)
                        refsel   = km_labelselection(tcfg.refmarker,data.label);
                        refsel	= ismember(data.label,refsel);
                    else
                        refsel = [];
                    end
                    
                    % simple (single axis) criterium
                    if ~ischar(tcfg.crit)
                        % get data
                        axissel = km_labelselection(tcfg.axis,data.axis);
                        axissel = ismember(data.axis,axissel);
                        pos = squeeze(data.pos{r}(markersel,axissel,:));
                        if ~isempty(refsel)
                            refpos = squeeze(data.pos{r}(refsel,axissel,:));
                            pos = pos - refpos;
                        end
                        
                        % go from the time to the sample domain
                        posidx = time2idx(mov,data.time{r});
                        
                        % apply criteria
                        idx = nan(size(posidx,1),1);
                        for i = 1:size(posidx,1)
                            if strcmpi(pname,'bpos')
                                tmp = pos(posidx(i,1));
                            elseif strcmpi(pname,'epos')
                                tmp = pos(posidx(i,2));
                            elseif strcmpi(pname,'pos')
                                tmp = pos(posidx(i,1):posidx(i,2));
                            end
                            if length(tcfg.crit) == 1
                                idx(i) = any(tmp > tcfg.crit);
                            else
                                idx(i) = any(tmp > min(tcfg.crit) & tmp < max(tcfg.crit));
                            end
                        end
                        
                        % complex (or multiple axes) criterium
                    else
                        % get data
                        pos = data.pos{r}(markersel,:,:);
                        if ~isempty(refsel)
                            refpos = data.pos{r}(refsel,:,:);
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
                        
                        % go from the time to the sample domain
                        posidx = time2idx(mov,data.time{r});
                        
                        % apply criterium
                        idx = nan(size(posidx,1),1);
                        for i = 1:size(posidx,1)
                            for a = 1:length(ax)
                                if strcmpi(pname,'bpos')
                                    tmp = axpos.(ax{a})(posidx(i,1));
                                elseif strcmpi(pname,'epos')
                                    tmp = axpos.(ax{a})(posidx(i,2));
                                elseif strcmpi(pname,'pos')
                                    tmp = axpos.(ax{a})(posidx(i,1):posidx(i,2));
                                end
                                eval(sprintf('%s = tmp;',ax{a}));
                            end
                            
                            % evaluate cirterium
                            idx(i) = any(eval(tcfg.crit));
                        end
                        
                    end
                    
                    % apply selection
                    idx = logical(idx);
                    if strcmpi(tcfg.type,'excl')
                        idx = ~idx;
                    end
                    
                otherwise
                    error('parameter ''%s'' not recognized',pname);
                    
            end
            
            % store checks
            if ismember('dur',cfg.movcheck(c).param)
                check{r}(idx_on,p,c) = ~idx;
            else
                check{r}(:,p,c) = ~idx;
            end
            
        end
        
        % expand the check to combined movements
        if ismember('dur',cfg.movcheck(c).param) && ~all(idx_on)
            idxcombon   = logic2idx(~idx_off);
            idxcomboff 	= logic2idx(~idx_on);
            for i = 1:size(idxcombon,1)
                if length(cfg.movcheck) > 1
                    error('Right, you wanted to run km_movcheck with multiple options while merging multiple movements. In other words, you have found an option that I have not implemented yet. Congratulations.');
                end
                idx = idxcombon(i,1):idxcomboff(i,2);
                check{r}(idx,:,c) = repmat(check{r}(idxcombon(i,1),:,c),length(idx),1);
            end
        end
        
    end
    
    % store only good movements
    if ismember('dur',cfg.movcheck(c).param)
        idx_select = any(~any(check{r}(idx_on,:,:),2),3);
    else
        idx_select = any(~any(check{r},2),3);
    end
    movement{r} = mov(idx_select,:);
        
end
cfg.movcheck(1).check = check;
cfg.movcheck(1).allmov = oldmovement;

% update dataset
if isempty(regexpi(cfg.dataset,'_mov','once'))
    if strcmpi(cfg.dataset,'_raw'),	cfg.dataset = '';   end
    cfg.dataset = [cfg.dataset '_MOV'];
end

% update processing directory
cfg = km_setcfg(cfg,'dirana');
cfg.dir.proc = cfg.dir.ana;

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
