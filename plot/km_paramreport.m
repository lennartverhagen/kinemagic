function km_paramreport(cfg,alldata)
%--------------------------------------------------------------------------
%
% See also KM_PLOT
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task,'strict');

% return if requested
if isfalse(cfg.(task))
    return
end


% select conditions and factors
if ~isfalse(cfg.plot.param)
    fact    = cfg.plot.fact;
    fname	= cfg.plot.fname;
else
    fname   = [];
end
% HACK lenver
%nparam  = length(cfg.plot.param);
% select variables
paramsel = match_str(cfg.vars,km_labelselection(cfg.paramreport.param,cfg.vars));
param = cfg.vars(paramsel);
nparam = length(param);
% END HACK
nfact	= length(fname);
nsubj   = length(cfg.subj);
if nsubj == 1 && ~any(ismember(cfg.paramreport.write,'subj'))
    warning('KM:ParamReport','data from only one subject supplied, changing write settings from ''%s'' to ''subj''.',cfg.paramreport.write{1});
    cfg.paramreport.write = {'subj'};
end

% individual subject report
if any(ismember(cfg.paramreport.write,{'subj','groupuni'}))
    
    % add rejection labels
    rejparam = get_rejparam(cfg);
    if ~any(ismember(rejparam,'reject')),	rejparam = [rejparam {'reject'}]; end
    rejsel = match_str(rejparam,km_labelselection(cfg.paramreport.reject,rejparam));
    rejparam = rejparam(rejsel);
    nrejparam = length(rejparam);
    if nrejparam > 0
        newrejparam = rejparam;
        for i = 1:nrejparam
            if ~strcmpi(rejparam{i},'reject')
                newrejparam{i} = sprintf('rej_%s',rejparam{i});
            end
        end
        newparam = [newrejparam param];
    end
    
    % loop over subjects
    subjhdr = cell(1,nsubj);
    subjdat = cell(1,nsubj);
    for s = 1:nsubj
        trl = cfg.trl{s}(:,paramsel);
        
        % HACK: TODO: do trlreject just like in all other cases
        
        % add rejection selection
        if nrejparam > 0
            % reject trials if requested
            idxr = km_getidxr(cfg,cfg.reportparam);
            
            trlrej = nan(size(trl,1),nrejparam);
            for i = 1:nrejparam
                % transform finite arrays with only ones and zeros to
                % locical arrays
                tmp = idxr.(rejparam{i});
                if isfloat(tmp) && isreal(tmp) && (sum(tmp(:)==1) + sum(tmp(:)==0) == numel(tmp))
                    idxr.(rejparam{i}) = logical(tmp);
                end
                trlrej(:,i) = cfg.idxr{s}.(rejparam{i});
            end
            
            % FIXME: adopt the listwise/pairwise option as used in
            % km_plot_getdata
            % exclude cases or whole trials with NaNs if requested
            if ~istrue(cfg.paramreport.exclude,'free') || any(strcmpi(cfg.paramreport.exclude,{'cases','param','indiv'}))
                for i = 1:nrejparam
                    if ~strcmpi(rejparam{i},'reject')
                        trl(logical(trlrej(:,i)),strcmpi(rejparam{i},param)) = NaN;                        
                    end
                end
            else
                trl(any(trlrej,2),:) = NaN;
            end
            
            trl = [trlrej trl];
        end
            
        % create header string
        hdrstr = sprintf('%s\t',newparam{:});
        % create data string
        expr = repmat('%0.4f\t',1,size(trl,2));
        datstr = sprintf([expr '\r\n'],trl');
        
        % store the report
        subjhdr{s} = newparam;
        subjdat{s} = trl;
        
        if any(ismember(cfg.paramreport.write,'subj'))
            % combine header and data into a report
            report = sprintf('%s\r\n%s',hdrstr,datstr);
            
            % save the report
            if istrue(cfg.save)
                % strings used for figure naming
                sessstr = cfg.plot.sessstr;
                % remove any low dash '_' in sessstr
                sessstr = strrep(sessstr,'_','');
                % append low dash '_' to end of sessstr
                if ~isempty(sessstr) && ~strcmpi(sessstr(end),'_')
                    sessstr = [sessstr '_'];
                end
                
                % define subject specific report directory
                cfg = km_setcfg(cfg,'dirreport');
                dir_report = getsubjsubdir(cfg,s,'report');
                % create file name and save
                reportname = fullfile(dir_report,sprintf('%strials.txt',sessstr));
                % open file
                fid = fopen(reportname,'w');
                % write file
                fprintf(fid,report);
                % close file
                fclose(fid);
            else
                % display the report in the matlab command window
                fprintf(1,'\nSubject: %s\n',cfg.subj{s});
                if ispc
                    disp(report);
                else
                    disp(regexprep(report,'\r',''));
                end
            end
        end
        
    end
    
end

% return if no group report
if ~any(ismember(cfg.paramreport.write,{'groupuni','grouprepmeas'})) || nsubj == 1
    return
end

% univariate group report
if any(ismember(cfg.paramreport.write,'groupuni'))
    % create header string
    hdr = [{'subj'} subjhdr{1}];
    hdrstr = sprintf('%s\t',hdr{:});
    % create data string
    datstr = cell(1,nsubj);
    for s = 1:nsubj
        % you need the subject number to know in which category the subject belongs
        subjnr = str2double(regexp(cfg.subj{s},'\d*','match','once'));
        if isempty(subjnr) || isnan(subjnr),    subjnr = s; end
        dat = [repmat(subjnr,size(subjdat{s},1),1) subjdat{s}];
        expr = repmat('%0.4f\t',1,size(dat,2));
        datstr{s} = sprintf([expr '\r\n'],dat');
        datstr{s} = datstr{s}(1:end-2);
    end
    % combine header and data into a report
    report = sprintf('%s\r\n',hdrstr,datstr{:});
    report = report(1:end-2);
    
    % save the report
    if istrue(cfg.save)
        % strings used for figure naming
        sessstr = cfg.plot.sessstr;
        % remove any low dash '_' in sessstr
        sessstr = strrep(sessstr,'_','');
        % append low dash '_' to end of sessstr
        if ~isempty(sessstr) && ~strcmpi(sessstr(end),'_')
            sessstr = [sessstr '_'];
        end
        
        % define group report directory
        cfg = km_setcfg(cfg,'dirreport');
        dir_report = getsubjsubdir(cfg,'Group','report');
        if exist(dir_report,'dir') ~= 7
            error('KM:DirReport','The group report directory does not exist. You are adviced to create that folder yourself: %',dir_report);
        end
        % create file name and save
        if isempty(sessstr)
            reportname = fullfile(dir_report,'group_univariate.txt');
        else
            reportname = fullfile(dir_report,sprintf('%sunivariate.txt',sessstr));
        end
        % open file
        fid = fopen(reportname,'w');
        % write file
        fprintf(fid,report);
        % close file
        fclose(fid);
    else
        % display the report in the matlab command window
        fprintf(1,'\n Group Univariate Report:\n');
        if ispc
            disp(report);
        else
            disp(regexprep(report,'\r',''));
        end
    end
end

% return if no group repeated measures report
if ~any(ismember(cfg.paramreport.write,'grouprepmeas'))
    return
end

% loop over variables
phdrstr = cell(nfact,nparam);
pdat = cell(nfact,nparam);
for p = 1:nparam
    
    % parameter name
    pname = param{p};
    
    % loop over factors for parameters
    for f = 1:nfact
        
        % HACK
        % get data descriptives
        %subjdescrip = cfg.plot.descrip.subj;
        %dat = alldata.(pname).(fname{f}).subj.(subjdescrip);
        dat = alldata.(pname).(fname{f}).subj;
        
        % get factor levels
        factlvl = cell(1,length(fact{f}));
        nlvl = nan(size(factlvl));
        for ff = 1:length(fact{f})
            factlvl{ff} = cfg.plot.lvl.(fact{f}{ff})(:)';
            nlvl(ff) = length(factlvl{ff});
        end
        
        % combine levels of all factors
        lvl = arraycomb(factlvl{:},'reverse');
        
        % check number of levels
        dims = size(dat);
        dims = dims(2:end);
        if ~isequal(nlvl,dims) && ~isequal([1 nlvl],dims)
            error('the number of levels [%s] does not match the number of data dimensions [%s]',num2str(nlvl),num2str(dims));
        end
        
        % create header string
        n = size(lvl,1);
        hdr = cell(1,n);
        phdr = cell(1,n);
        for i = 1:n
            tmp = [fact{f};num2cell(lvl(i,:))];
            hdr{i} = sprintf('%s%d',tmp{:});
            phdr{i} = sprintf('%s_%s',upper(pname),hdr{i});
        end
        hdrstr = sprintf('%s\t',hdr{:});
        
        % create header string for appending other parameters
        phdrstr{f,p} = sprintf('%s\t',phdr{:});
            
        % permute data to fit header
        dat = permute(dat,[1 ndims(dat):-1:2]);
        
        % reshape data to fit nsubj x n matrix
        dat = reshape(dat,[nsubj n]);
        
        % create data cell array for appending other parameters
        pdat{f,p} = dat;
        
        % create data string
        expr = repmat('%0.4f\t',1,n);
        datstr = sprintf([expr '\r\n'],dat');
        
        % combine header and data into a report
        report = sprintf('%s\r\n%s',hdrstr,datstr);
        
        % display the report in the matlab command window
        fprintf(1,'\n Group Repeated Measures Report:\n');
        %fprintf(report);
        if ispc
            disp(report);
        else
            disp(regexprep(report,'\r',''));
        end
        
        % combine parameter reports
        if nparam > 1 && p == nparam
            % create header string 
            hdrstr = sprintf('%s',phdrstr{f,:});
            % create data string
            dat = [pdat{f,:}];
            expr = repmat('%0.4f\t',1,size(dat,2));
            datstr = sprintf([expr '\r\n'],dat');
            % combine header and data into a report
            preport = sprintf('%s\r\n%s',hdrstr,datstr);
            % display the report in the matlab command window
            if ispc
                disp(preport);
            else
                disp(regexprep(preport,'\r',''));
            end
        end
        
        % save the report
        if istrue(cfg.save)
            % strings used for figure naming
            sessstr = cfg.plot.sessstr;
            % remove any low dash '_' in sessstr
            sessstr = strrep(sessstr,'_','');
            % append low dash '_' to end of sessstr
            if ~isempty(sessstr) && ~strcmpi(sessstr(end),'_')
                sessstr = [sessstr '_'];
            end
            
            % define group report directory
            cfg = km_setcfg(cfg,'dirreport');
            dir_report = getsubjsubdir(cfg,'Group','report');
            if exist(dir_report,'dir') ~= 7
                error('KM:DirReport','The group report directory does not exist. You are adviced to create that folder yourself: %',dir_report);
            end
            % create file name and save
            reportname = fullfile(dir_report,sprintf('%s%s_%s.txt',sessstr,fname{f},pname));
            % open file
            fid = fopen(reportname,'w');
            % write file
            fprintf(fid,report);
            % close file
            fclose(fid);
            
            % save the combined parameter report
            if nparam > 1 && p == nparam
                pname = sprintf('%s_',cfg.plot.param{:});
                pname = pname(1:end-1);
                % create file name and save
                reportname = fullfile(cfg.dir.report,sprintf('%s%s_%s.txt',sessstr,fname{f},pname));
                % open file
                fid = fopen(reportname,'w');
                % write file
                fprintf(fid,preport);
                % close file
                fclose(fid);
            end
        end
        
    end
    
end
%--------------------------------------------------------------------------


% function get_rejparam
%----------------------------------------
function rejparam = get_rejparam(cfg)

if isfield(cfg,'idxr')
    if iscell(cfg.idxr)
        idxr = cfg.idxr{1};
    else
        idxr = cfg.idxr;
    end
    fnames = fieldnames(idxr);
    for f = 1:length(fnames)
        tmp = idxr.(fnames{f});
        % allow only logical arrays
        if ~islogical(tmp)
            % and allow finite arrays with only ones and zeros
            if ~(isfloat(tmp) && isreal(tmp) && (sum(tmp(:)==1) + sum(tmp(:)==0) == numel(tmp)))
                fnames{f} = '';
            end
        end
    end
    rejparam = fnames(~ismember(fnames,{'','zeros','descrip','report','save','sess','dir','dataset'}))';
else
    rejparam = cfg.plot.trlrej.mode;
end
%--------------------------------------------------------------------------
