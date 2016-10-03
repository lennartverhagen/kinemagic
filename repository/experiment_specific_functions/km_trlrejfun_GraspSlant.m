function idx = km_trlrejfun_GraspSlant(cfg,~,idx)
%--------------------------------------------------------------------------
% KM_TRLREJFUN_GRASPSLANT selects cases to be rejected based on
% behavioural parameters
%
% See also KM_TRLREJ
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------


% % you need the subject number to know in which category the subject belongs
% subjnr = str2double(regexp(cfg.subj,'\d*','match','once'));    

% get trial data and variable names
vars = cfg.vars;
trl  = cfg.trl;

% loop over automatic trial rejection parameters
for p = 1:length(cfg.trlrej.param)
    
    % param name
    pname = cfg.trlrej.param{p};
    if strcmpi(pname,'apriori'),	continue;	end
    if ~any(strcmp(vars,pname))
        warning('KM:NotSupported','KM %s: Parameter ''%s'' is not supported.',mfilename,pname);
        continue
    end
    
    % get parameter
    dat = trl(:,strcmp(vars,pname));
    
    % get cutoffs
    initcut = get_initcut(pname,cfg.trlrej);
        
    % select trials with too low or high parameter values
    [idx.describ.(pname).low, idxextrlow, cutlow]       = cutoutlier(dat,3,'iqr','lower',initcut,idx.apriori);
    [idx.describ.(pname).high, idxextrhigh, cuthigh]	= cutoutlier(dat,3,'iqr','higher',initcut,idx.apriori);
    idx.describ.(pname).cutlow = cutlow;
    idx.describ.(pname).cuthigh = cuthigh;
    idxextreme = idxextrlow | idxextrhigh;
    
    % combine low-high rejections
    idx.(pname) = idx.describ.(pname).low | idx.describ.(pname).high;
    
    % calculate descriptive metrics
    step = (cuthigh-cutlow)/25;
    edges = (cutlow:step:(cuthigh+step)) - step/2;
    [~,bin] = histc(dat,edges); bin(bin==0) = NaN; modbin = mode(bin);
    if ~isnan(modbin),	modbin = edges(modbin) + step/2;	end
    idx.describ.(pname).all.values	= dat;
    idx.describ.(pname).all.mean    = nanmean(dat);
    idx.describ.(pname).all.median  = nanmedian(dat);
    idx.describ.(pname).all.mode    = modbin;
    [~,bin] = histc(dat(~idxextreme),edges); bin(bin==0) = NaN; modbin = mode(bin);
    if ~isnan(modbin),	modbin = edges(modbin) + step/2;	end
    idx.describ.(pname).sel.values  = dat(~idxextreme);
    idx.describ.(pname).sel.mean    = nanmean(dat(~idxextreme));
    idx.describ.(pname).sel.median  = nanmedian(dat(~idxextreme));
    idx.describ.(pname).sel.mode    = modbin;
    
end
%--------------------------------------------------------------------------

% function get_initcut
%----------------------------------------
function cut = get_initcut(pname,cfg,flgrecursive)
if nargin <3, flgrecursive = false; end

%look for parameter name in list with cutoffs
if isfield(cfg,pname)
    cut = cfg.(pname); return;
end

% test parameter name against experimental default settings
switch lower(pname)
    case 'rt',          cut	= [0.3  1.5 ];
    case 'mt',          cut	= [0.4  2.0 ];
    case {'mt_t','mt_tt','mt_ti'},      cut	= [0.4  1.5];
    case {'rmt_t','rmt_tt','rmt_ti'},	cut	= [0.5  1  ];
    case {'mt_a','mt_at','mt_ai'},     	cut	= [0    0.8];
    case {'rmt_a','rmt_at','rmt_ai'},  	cut	= [0    0.5];
    case 'tl',          cut	= [0.35 0.6 ];
    case 'tl_t',        cut	= [0.3  0.55];
    case 'tl_a',        cut	= [0    0.08];
    case 'mv',          cut	= [0.15 1   ];
    case 'mv_t',        cut	= [0.2  1   ];
    case 'mv_a',        cut	= [0    0.2 ];
    case 'pv',          cut	= [0.5  2   ];
    case 'tpv',         cut	= [0.1  0.75];
    case 'rtpv',        cut	= [0.1  0.5 ];
    case 'pa',          cut	= [0    30 ];
    case 'tpa',         cut	= [0    0.5 ];
    case 'rtpa',        cut	= [0    0.4 ];
    case 'pd',          cut	= [-20  0   ];
    case 'tpd',         cut	= [0.2  2   ];
    case 'rtpd',        cut	= [0.1  1   ];
    case 'nvc',         cut	= [1    35  ];
    case 'nvc_t',       cut	= [1    25  ];
    case 'nvc_a',       cut	= [0    20  ];
    case 'pgv',         cut	= [0    1   ];
    case 'npgv',        cut	= [-2   0   ];
    case 'gv_a',        cut	= [-0.5 0.3 ];
    case 'mgv_a',       cut	= [-0.3 0.1 ];
    case 'zgv_a',       cut	= [-0.1 0.05];
    case 'pga',         cut	= [0.08 0.16];
    case 'tpga',        cut	= [0.2  1.5 ];
    case 'rtpga',       cut	= [0.4  1   ];
    case 'ga_a',        cut	= [0.075 0.15];
    case 'dga_a_re',    cut	= [-0.02 0.06];
    case 'ga_e',        cut	= [0.08 0.11];
    case 'dga_e_re',    cut	= [-0.0025 0.015];
    case 'go_a',        cut	= [-30  120 ];
    case 'dgo_a_re',    cut	= [-25  25  ];
    case 'adgo_a_re',   cut	= [-Inf 20  ];
    case 'go_e',        cut	= [-20  100 ];
    case 'dgo_e_re',    cut	= [-7   7   ];
    case 'adgo_e_re',   cut	= [-Inf 7   ];
    otherwise
        if ~flgrecursive
            % remove trailing sensors and axis indicators from pname
            pname = regexp(pname,'^[a-zA-z]+_?[A-Z]?[a-z]*_?[A-Z]?[a-z]*','match');
            % and run again
            cut = get_initcut(pname{1},cfg,true);
        else
            cut	= [-Inf Inf	];
        end
end