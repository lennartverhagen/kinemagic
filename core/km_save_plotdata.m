function km_save_plotdata(cfg,param,tseries)
%--------------------------------------------------------------------------
%
% See also KM_DOTASK
%
% This file is part of the FieldTripWrapper toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% save configuration and data
if isfalse(cfg.save)
    return
end
    
% check configuration
cfg = km_setcfg(cfg,{'subj','sess','dataset'});
cfg = km_setcfg(cfg,'dirreport');
dir_report = getsubjsubdir(cfg,'Group','report');

% get subj and sess strings
if iscell(cfg.subj)
    if length(cfg.subj) == 1
        str_subj = cfg.subj{1};
    else
        str_subj = 'Group';
    end
else
    str_subj = cfg.subj;
end
str_sess = [cfg.sess{:}];

% save configuration
if istrue(cfg.save,'free') || isequal(cfg.save,'cfg')
	fname = sprintf('%s%s%s',str_subj,str_sess,cfg.dataset);
    save(fullfile(dir_report,fname),'cfg');
end

% save param and tseries data
if istrue(cfg.save,'free') || isequal(cfg.save,'data')
	fname = sprintf('%s%s%s_PARAM',str_subj,str_sess,cfg.dataset);
    save(fullfile(dir_report,fname),'param');
	fname = sprintf('%s%s%s_TSERIES',str_subj,str_sess,cfg.dataset);
    save(fullfile(dir_report,fname),'tseries');
end