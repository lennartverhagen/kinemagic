function idx = km_aprioritrlrejfun_GraspSlant(cfg,data,idx)
%--------------------------------------------------------------------------
% KM_APRIORITRLREJFUN_GRASPSLANT selects trials to be rejected based on
% behavioural parameters
%
% See also KM_TRLREJ
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------


% you need the subject number to know in which category the subject belongs
subjnr = str2double(regexp(cfg.subj,'\d*','match','once'));    

% get trldata
trl     = cfg.trl;

% number of trials
ntrl = size(trl,1);
ndattrl = length(data.time);

% check number of trials
if ntrl ~= ndattrl
    error('number of trials in the data does not match the number of trials in the logfile');
end

% reject trials a priori
idx.apriori = false(ntrl,1);
if subjnr == 999 && strcmpi(cfg.sess,'example')
    idx.apriori(1:12) = true;
end

% exclude trials with rt < ttms
%rt = cfg.trl(:,strcmp(cfg.vars,'rt'));
%ttms = cfg.trl(:,strcmp(cfg.vars,'TTMS'))/1000;
%idx.apriori = idx.apriori | rt<=ttms;