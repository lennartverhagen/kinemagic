function str = km_trlrej_report(cfg)
%--------------------------------------------------------------------------
% provide a textual report of the rejected trials
%
% See also KM_TRLREJ
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% get input
idx = cfg.trlrej.idx;

% number of trials
if isfield(idx,'zeros')
    ntrl = length(idx.zeros);
elseif isfield(idx,'apriori')
    ntrl = length(idx.apriori);
elseif isfield(idx,'nan')
    ntrl = length(idx.nan);
else
    warning('KM:NumberOfTrials','FTW %s: The number of trials could not be determined',mfilename);
    str = '';
    return
end

% automatically assign reportparam if requested
if ~isfalse(cfg.trlrej.reportparam) && strcmpi(cfg.trlrej.reportparam{1},'auto')
    cfg.trlrej.reportparam = fieldnames(idx);
    cfg.trlrej.reportparam = cfg.trlrej.reportparam(~ismember(cfg.trlrej.reportparam,'zeros'));
end

% organize and combine rejection indeces
str = sprintf('Subject %s',cfg.subj);
if ~isempty(cfg.sess)
    sess = cfg.sess; if strcmpi(sess(1),'_'),	sess = sess(2:end);	end
    str = sprintf('%s, Sess %s',str,sess);
end
str = sprintf('%s\n%s\n',str,repmat('-',1,length(str)));
idx.reject = zeros(ntrl,1);
for i = 1:length(cfg.trlrej.reportparam)
    pname = cfg.trlrej.reportparam{i};
    switch lower(pname)
        case 'apriori'
            str = sprintf('%s    apriori:\t\t',str);
        case 'rt'
            str = sprintf('%s    reaction time:\t',str);
        case 'mt'
            str = sprintf('%s    movement time:\t',str);
        case 'tl'
            str = sprintf('%s    trajectory length:\t',str);
        case 'mv'
            str = sprintf('%s    mean velocity:\t',str);
        case 'pv'
            str = sprintf('%s    peak velocity:\t',str);
        case 'pga'
            str = sprintf('%s    max grip aperture:\t',str);
        case 'ga_e'
            str = sprintf('%s    end grip aperture:\t',str);
        case 'dgo_e_re'
            str = sprintf('%s    delta end grip ori:\t',str);
        otherwise
            str = sprintf('%s    %s:\t',str,pname);
            if length(pname) < 12
                str = sprintf('%s\t',str);
            end
            if length(pname) < 4
                str = sprintf('%s\t',str);
            end
            %warning('KM:NotSupported','KM %s: Mode ''%s'' is not supported.',mfilename,pname);
    end
    
    if ~isfield(idx,pname)
        tmpstr = sprintf('not performed\n');
    else
        idx.reject = idx.reject | idx.(pname);
        tmpstr = sprintf('%3d trials',sum(idx.(pname)));
        %if sum(idx.(pname)) > 0
            p = round(percof(sum(idx.(pname)),ntrl));
            tmpstr = sprintf('%s (%d%%)\n',tmpstr,p);
        %else
        %    tmpstr = sprintf('%s\n',tmpstr);
        %end
    end
    str = [str tmpstr];
end

% all criteria combined
str = sprintf('%s    %s\n',str,repmat('-',1,20+length(tmpstr)-1));
str = sprintf('%s    total:\t\t%3d trials\n',str,ntrl);
p = round(percof(sum(idx.reject),ntrl));
str = sprintf('%s    rejections:\t\t%3d trials (%d%%)\n',str,sum(idx.reject),p);
p = round(percof(sum(~idx.reject),ntrl));
str = sprintf('%s    survivors:\t\t%3d trials (%d%%)\n',str,sum(~idx.reject),p);

% display the report in the command window
fprintf(1,'\n%s',str);
