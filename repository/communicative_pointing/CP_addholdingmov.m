%--------------------------------------------------------------------------
% This script defines and adds holding time to movement params.
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------

% what subject to analyse
subjsel    	= [1:15];
nsubj       = length(subjsel);

% settings
dataset = 'FILTlp15_MOV_ANA';
dir_proc = fullfile(filesep,'home','action','ankmur','CommPointing','Data-analysis','KIN','ana');
subjects        = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15'};

% loop over subjects
for s = 1:nsubj
    
    % select the correct file
    subj = subjects{subjsel(s)};
    fname_cfg = fullfile(dir_proc,sprintf('%s_%s_CFG',subj,dataset));
    
    % load the correct files
    cfg = load(fname_cfg);
    cfg = cfg.cfg;
    
    % get movement
    if isfield(cfg,'movement') && ~isstruct(cfg.movement) && ~isempty(cfg.movement)
        movement = cfg.movement;
    elseif isfield(cfg,'movdef') && isfield(cfg.movdef,'movement')
        movement = cfg.movdef.movement;
    elseif isfield(data,'movement')
        movement = data.movement;
    elseif isfield(cfg,'movdef')
        movement = km_movdef(cfg,data);
    else
        error('no movement provided');
    end
    % make sure that all trials have the same number of movements
    movement = km_setequalnmov(movement);
    
    % check number of movements
    nmov = size(movement{1},1);
    if nmov == 3
        continue
    elseif s ~= 13 && (nmov < 2 || nmov > 3)
        error('something is wrong');
    end
    
    % add holding time as a movement
    cfg.movement = cellfun(@(x) [x(1:2,:); x(1,2) x(2,1)],movement,'UniformOutput',false);

    % save
    cfg.save = 'cfg';
    km_save(cfg);
    
end