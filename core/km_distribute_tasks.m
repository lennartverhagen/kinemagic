function varargout = km_distribute_tasks(cfg)
% distribute the tasks (function to be called by qsubcellfun)
%----------------------------------------

% loop over sessions
tcfg = cell(1,length(cfg.sessions));
tdata = cell(1,length(cfg.sessions));
for ss = 1:length(cfg.sessions)
    cfg.sess = cfg.sessions{ss};
    if isequal(cfg.task,{'traject'})
        % this is a special case where trajectory data from indivudual
        % sessions is collapsed over sessions.
        cfg_traject = km_dotask(cfg);
        if isfield(cfg_traject,'tmp')
            cfg.tmp = cfg_traject.tmp;
        end
        if isfield(cfg.traject,'sesscounter')
            cfg.traject.sesscounter = cfg_traject.traject.sesscounter;
        end
        tcfg{ss} = cfg;
    elseif ~isfalse(cfg.task)
        % call km_dotask depending on requested number of outputs
        if nargout == 0
            km_dotask(cfg);
        elseif nargout == 1
            tcfg{ss} = km_dotask(cfg);
        elseif nargout > 1
            [tcfg{ss},tdata{ss}] = km_dotask(cfg);
        end
    end    
end

% sort output
if nargout == 0
    varargout = {};
elseif nargout == 1
    varargout = {tcfg};
elseif nargout > 1
    varargout = {tcfg,tdata};
end