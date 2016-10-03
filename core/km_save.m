function km_save(cfg,data)
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
cfg = km_setcfg(cfg,'dirproc');

% if ~isfield(cfg,'subj'),        error('no subject specified: cfg.subj');    end
% if ~isfield(cfg,'sess'),        error('no session specified: cfg.sess');    end
% if ~isfield(cfg,'dataset'),     cfg.dataset     = '';                       end     % 'MASK', 'FILT'
% if ~isempty(cfg.dataset) && ~strcmpi(cfg.dataset(1),'_')
%     cfg.dataset = ['_' cfg.dataset];
% end

% WARNING: Maybe the line below is not the best solution
% check if directory exists, and otherwise, create it
if ~exist(cfg.dir.proc,'dir'), mkdir(cfg.dir.proc); end

% save configuration
if istrue(cfg.save,'free') || isequal(cfg.save,'cfg')
    if ~isempty(regexpi(cfg.dataset,'IDXR$|IDXREJ$|IDXC$|IDXCOND$','once'))
        fname = sprintf('%s%s%s',cfg.subj,cfg.sess,cfg.dataset);
    else
        fname = sprintf('%s%s%s_CFG',cfg.subj,cfg.sess,cfg.dataset);
    end
    save(fullfile(cfg.dir.proc,fname),'cfg');
end

% save data
if istrue(cfg.save,'free') || isequal(cfg.save,'data')
    fname = sprintf('%s%s%s_DATA',cfg.subj,cfg.sess,cfg.dataset);
    save(fullfile(cfg.dir.proc,fname),'data');
end