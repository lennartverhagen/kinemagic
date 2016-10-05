function [cfg,flg_return] = CP_movparamfun(cfg,~,flg_caller)
%--------------------------------------------------------------------------
% This function fills in experimentally specific fields in
% the cfg and data structures before and after km_movparam
%
% See also KM_MOVPARAM
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------

% if called at the beginning of km_movparam return quickly
flg_return = false;
if ~strcmpi(flg_caller,'end'), return; end

% retrieve fields
vars = cfg.vars;
trl = cfg.trl;

% is holding time already exists, remove it
idx_ht = strcmpi('ht',vars);
if any(idx_ht)
    vars = vars(~idx_ht);
    trl = trl(:,~idx_ht);
end
idx_pname_ht = strcmpi('ht',cfg.movparam.pname);
if any(idx_pname_ht)
    cfg.movparam.pname = cfg.movparam.pname(~idx_pname_ht);
end

% adapt rt_2 to "holding time" (ht)
idx_rt_1 = strcmpi('rt_1',vars);
idx_mt_1 = strcmpi('mt_1',vars);
idx_rt_2 = strcmpi('rt_2',vars);
ht = trl(:,idx_rt_2) - trl(:,idx_rt_1) - trl(:,idx_mt_1);

% store
if isempty(ht)
    error('KM:MovParam:ExpFunKIN','Holding time (ht) could not be correctly calculated');
end
vars{idx_rt_2} = 'ht'; cfg.vars = vars;
trl(:,idx_rt_2) = ht; cfg.trl = trl;
cfg.movparam.pname{strcmpi('rt_2',cfg.movparam.pname)} = 'ht';