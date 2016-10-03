function km_paramplot(cfg,alldata)
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
if ~isfield(cfg,task) && ~isempty(alldata)
    cfg.(task) = 'yes';
end
cfg = km_setcfg(cfg,task,'strict');

% return if requested
if isfalse(cfg.(task)) || isempty(alldata)
    return
end


% select conditions and factors
fact    = cfg.plot.fact;
fname	= cfg.plot.fname;
nparam  = length(cfg.plot.param);
nfact	= length(fname);
nsubj   = length(cfg.subj);

% loop over parameters
for p = 1:nparam
    
    % parameter name
    pname = cfg.plot.param{p};
    pdata = alldata.(pname);
    
    % loop over factors for parameters
    for f = 1:nfact
        tmpfact = fact{f}(~ismember(fact{f},cfg.plot.contr{f}));
        
        % get plot type
        plottype = cfg.paramplot.type{min(length(cfg.paramplot.type),f)};
        
        % get plotting specifics
        spec.flcl       = str2rgb(cfg.plot.flcl{f});
        spec.lncl       = str2rgb(cfg.plot.lncl{f});
        spec.lnstyle	= cfg.plot.lnstyle{f};
        spec.lnwidth	= cfg.plot.lnwidth{f};
        spec.marker     = cfg.plot.marker{f};
        
        % get data descriptives
        %subjdescrip = cfg.plot.descrip.subj;
        %groupdescrip = cfg.plot.descrip.group;
        %if nsubj == 1
        %    dat = pdata.(fname{f}).(subjdescrip);
        %else
        %    dat = pdata.(fname{f}).(groupdescrip).(subjdescrip);
        %end
        dat = pdata.(fname{f}).dat;
        
        % reshape dat to fit in one figure
        dims = size(dat);
        dims = [prod(dims(1:end-1)) dims(end)];
        dat = reshape(dat,dims);
        
        % initiate new figure
        h = figure;
        
        % get levels
        lvl = cell(1,length(tmpfact));
        nlvl = nan(size(lvl));
        for ff = 1:length(tmpfact)
            lvl{ff} = cfg.plot.lvl.(tmpfact{ff})(:)';
            nlvl(ff) = length(lvl{ff});
        end
        if isempty(lvl),    lvl = {1};  end
        clvl = arraycomb(lvl{1:end-1});
        
        % bar plot
        if strcmpi(plottype,'bar')
            if numel(dat) == 1
                hb = bar(lvl{end}(:),dat',0.6);
                xlim([lvl{end}-1 lvl{end}+1]);
            elseif numel(dat) == 2
                hb = bar(lvl{end}(:),dat',0.6);
                xlim([min(lvl{end})-abs(diff(lvl{end})/2) max(lvl{end})+abs(diff(lvl{end})/2)]);
            else
                hb = bar(lvl{end}(:),dat');
            end
            xdat = getbarx(hb)';
            if length(hb) > 1
                for i = 1:length(hb);
                    flcl  = spec.flcl(i,:);
                    lncl  = spec.lncl(i,:);
                    set(hb(i),'FaceColor',flcl,'EdgeColor',lncl);
                    if ~all(cellfun(@(x) isequal(spec.lncl(1,:),x),num2cell(spec.lncl,2)))
                        set(hb(i),'LineWidth',2);
                    end
                end
            end
            hold on
        else
            % get x points for line plot (with offset)
            %shiftfact = 1/100;
            shiftfact = 0;
            xdat = lvl{end};
            if length(tmpfact) > 1
                flvl = lvl{end-1};
                d = max(xdat)-min(xdat);
                xdat = repmat(xdat,dims(1),1);
                shift = mc(flvl,[min(flvl) max(flvl)])*floor(length(flvl)/2);
                shift = shiftfact .* d .* shift(:);
                [~,~,tmplvl] = unique(clvl(:,end));
                shift = shift(tmplvl);
                xdat = xdat + repmat(shift,1,size(xdat,2));
            end
        end
        
        % plot with errorbars
        if ~isfalse(cfg.plot.descrip.err)
            % get error
            %errdescrip = cfg.plot.descrip.err;
            if nsubj == 1
                %err = pdata.(fname{f}).(errdescrip);
                err = pdata.(fname{f}).err;
            elseif ~isfalse(cfg.plot.descrip.contrerr) && ...
                    ( (length(nlvl) > 1 && any(nlvl(1:end-1) == 2)) || ...
                    (length(nlvl) == 1 && any(nlvl == 2)) )
                %err = pdata.(fname{f}).diff.(errdescrip).(subjdescrip);
                err = pdata.(fname{f}).diff.err;
                dim = find(nlvl == 2,1,'first');
                tile = ones(1,length(nlvl));
                tile(dim) = 2;
                if length(tile) == 1,   tile = [1 tile];    end
                err = repmat(err,tile);
            else
                %err = pdata.(fname{f}).(errdescrip).(subjdescrip);
                err = pdata.(fname{f}).err;
            end
            
            % reshape err to fit in one figure
            err = reshape(err,dims);
            
            % adapt plottype
            if strcmpi(plottype,'bar')
                % add the errorbars to the bar plot
                plottype = 'errorbar';
                % no lines and markers though
                spec.lnstyle = repmat({'no'},1,length(spec.lnstyle));
                spec.marker = repmat({'no'},1,length(spec.marker));
                % and black errorbars
                spec.flcl = repmat([0 0 0],size(spec.flcl,1),1);
            end
        end
        
        if ~strcmpi(plottype,'bar')
            for i = 1:size(dat,1)
                
                % plotting specifics
                flcl  = spec.flcl(i,:);
                lnstyle     = spec.lnstyle{i};
                lnwidth     = spec.lnwidth(i);
                marker      = spec.marker{i};
                
                switch plottype
                    case 'errorbar'
                        errorbar(xdat(i,:),dat(i,:),err(i,:),'Color',flcl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                    case 'errorcloud'
                        % plot data
                        plot(xdat(i,:),dat(i,:),'Color',flcl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                        % plot error cloud in same color
                        ploterrorcloud(xdat(i,:),dat(i,:),err(i,:),flcl);
                    otherwise
                        plot(xdat(i,:),dat(i,:),'Color',flcl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                end
                hold on
            end
        end
        hold off
        
        % HACK: reset xlim
        if any(ismember(cfg.plot.bfact,'slant'))
            xlim([-7.5 97.5]);
            set(gca,'XTick',0:15:90);
        end
        
        % adjust y limits
        yticks = getaxrange(dat,err);
        ylim([min(yticks) max(yticks)]);
        
        % set baseline value of bar plots
        if exist('hb','var') && all(ishandle(hb))
            if min(yticks) >= 0
                set(hb,'BaseValue',min(yticks));
                set(get(hb(1),'BaseLine'),'LineStyle','no')
            elseif max(yticks) <= 0
                set(hb,'BaseValue',max(yticks));
                set(get(hb(1),'BaseLine'),'LineStyle','no')
            else
                set(hb,'BaseValue',0);
                set(get(hb(1),'BaseLine'),'LineWidth',1.5,'LineStyle','--')
            end
            
            % set bar face and edge colors for ungrouped bars
            if length(hb) == 1
                hp = get(hb,'Children');
                % set individual bar face colors
                % PRINT does not yet support RGB colormode, therefore a
                % colormap is used
                % set(hp,'CData',shiftdim(spec.flcl,-1))
                set(hp,'CDataMapping','direct','CData',1:length(dat))
                colormap(str2rgb(cfg.plot.flcl{f}));
                % check if bars should have individually colored edges
                if ~all(cellfun(@(x) isequal(spec.lncl(1,:),x),num2cell(spec.lncl,2)))
                    % plot new patches with only edges over the original bars
                    hold on;
                    % PRINT does not yet support RGB colormode, therefore
                    % a colormap is used
                    % lncl = repmat(shiftdim(spec.lncl,-1),4,1);
                    % hpedge = patch(get(hp,'XData'),get(hp,'YData'),zeros(size(get(hp,'XData'))),'FaceColor','none','EdgeColor','flat','CData',lncl,'LineWidth',2);
                    hpedge = patch(get(hp,'XData'),get(hp,'YData'),zeros(size(get(hp,'XData'))),'FaceColor','none','EdgeColor','flat','CData',length(dat)+repmat(1:length(dat),4,1),'LineWidth',2);
                    colormap([spec.flcl; spec.lncl]);
                    uistack(hpedge,'down',1);
                else
                    % set single edge color
                    set(hp,'EdgeColor',spec.lncl(1,:));
                end
            end
        end
        
        % set y-ticks
        %set(gca,'YTickMode','auto')
        set(gca,'YTick',yticks);
        
        % pimp plot
        km_pimpplot(h)
                 
        % wrap up the figure
        fig_wrapup(cfg,h,pname,fname{f});
        
    end
    
end
%--------------------------------------------------------------------------


%% function fig_wrapup
%----------------------------------------
function fig_wrapup(cfg,h,paramstr,namestr)
% when done with plotting...

if ishandle(h)
    figure(h);
    
    % strings used for figure naming
    subjstr = cfg.plot.subjstr;
    sessstr = cfg.plot.sessstr;
    % append low dash '_' to end of sessstr
    if ~isempty(sessstr) && ~strcmpi(sessstr(end),'_')
        sessstr = [sessstr '_'];
    end
    
    % make the figure background white and give the figure a name
    figname = sprintf('cond:%s  param:%s  subj:%s  sess:%s',spellout(namestr,1),paramstr,subjstr,sessstr);
    set(h,'Color','w'); set(gca,'Box','off');
    set(h,'NumberTitle','off','Name',figname);
    
    % set the y-axis label
    switch lower(paramstr)
        case 'rt',  fullparamstr = 'reaction time';
        case 'mt',  fullparamstr = 'movement duration';
        case 'tt',  fullparamstr = 'transport duration';
        case 'rtt', fullparamstr = 'rel. transport duration';
        case 'gt',  fullparamstr = 'approach duration';
        case 'rgt', fullparamstr = 'rel. approach duration';
        case 'ttt', fullparamstr = 'thumb transport duration';
        case 'rttt',fullparamstr = 'rel. thumb transport duration';
        case 'itt', fullparamstr = 'index transport duration';
        case 'ritt',fullparamstr = 'rel. index transport duration';
        case 'itdiff',  fullparamstr = 'diff. index-thumb duration';
        case 'ritdiff', fullparamstr = 'rel. diff. index-thumb duration';
        case 'td',	fullparamstr = 'total trajectory';
        case 'tcd',	fullparamstr = 'transport trajectory';
        case 'gcd',	fullparamstr = 'approach trajectory';
        case 'mv',	fullparamstr = 'mean velocity';
        case 'mtv',	fullparamstr = 'mean transport velocity';
        case 'mgv',	fullparamstr = 'mean appraoch velocity';
        case 'pv',	fullparamstr = 'peak veloctiy';
        case 'tpv',	fullparamstr = 'time to peak velocity';
        case 'rtpv',fullparamstr = 'rel. time to peak velocity';
        case 'mga',	fullparamstr = 'maximum grip aperture';
        case 'tmga',fullparamstr = 'time to maximum grip aperture';
        case 'rtmga',fullparamstr = 'rel. time to maximum grip aperture';
        case 'tga',	fullparamstr = 'approach grip aperture';
        case 'dtga',fullparamstr = 'diff. approach grip aperture';
        case 'tgo',	fullparamstr = 'approach grip orientation';
        case 'dtgo',fullparamstr = 'diff. approach grip orientation';
        case 'adtgo',fullparamstr = 'abs. diff. approach grip orientation';
        case 'tgv',	fullparamstr = 'approach grip velocity';
        case 'ega',	fullparamstr = 'end grip aperture';
        case 'dega',fullparamstr = 'diff. end grip aperture';
        case 'ego',	fullparamstr = 'end grip orientation';
        case 'dego',fullparamstr = 'diff. end grip orientation';
        case 'adego',fullparamstr = 'abs. diff. end grip orientation';
        case 'test',	fullparamstr = 'test';
        otherwise,	fullparamstr = paramstr;
    end
    switch lower(paramstr)
        case {'rt','mt','tt','gt','ttt','itt','itdiff','tpv','tpa','tpd','tmga'}
            unitstr = 's';
        case {'td','tcd','gcd','mga','tga','dtga','ega','dega'}
            unitstr = 'log(m)';
        case {'mv','mtv','mgv','pv','tgv'};
            unitstr = 'm/s';
        case {'pa','pd'};
            unitstr = 'm/s²';
        case {'rtt','rgt','rttt','ritt','ritdiff','rtpv','rtpa','rtpd','rtmga'};
            unitstr = '%mt';
        case {'tgo','dtgo','adtgo','ego','dego','adego'};
            unitstr = 'deg';
        otherwise
            unitstr = 'a.u.';
    end
    ystr = sprintf('%s (%s)',fullparamstr,unitstr);
    set(get(gca,'YLabel'),'String',ystr,'FontSize',14,'interpreter','none');
    
    % update the factor based on the contrast
    f = find(strcmpi(cfg.plot.fname,namestr));
    tmpfact = cfg.plot.fact{f}(~ismember(cfg.plot.fact{f},cfg.plot.contr{f}));
    
    % set the x-axis label
    if ~isempty(tmpfact)
        xstr = tmpfact{end};
    else
        xstr = '';
    end
    set(get(gca,'XLabel'),'String',xstr,'FontSize',14,'interpreter','none');
    
    % tick labels and legend
    if isfield(cfg.plot,'lvltick')
        % set the x-tick labels
        allticks = cell(1,length(tmpfact));
        lvl = cell(1,length(tmpfact));
        for ff = 1:length(tmpfact)
            fname = tmpfact{ff};
            allticks{ff} = cfg.plot.lvltick.(fname);
            lvl{ff} = cfg.plot.lvl.(fname)(:)';
        end
        if isempty(lvl)
            lvl = {1};
            allticks = {{'contrast'}};
        end
        lvls = lvl{end};
        if min(lvls) >= 1 && max(lvls) <= length(allticks{end})
            set(gca,'XTickLabel',allticks{end}(lvls));
        end
        
        % create a legend
        if length(tmpfact) == 2
            legend(gca,allticks{1}{lvl{1}},'Location','Best');
        elseif length(tmpfact) > 2
            clvl = arraycomb(lvl{1:end-1});
            legendstr = cell(1,size(clvl,1));
            for i = 1:size(clvl,1)
                legendstr{i} = '';
                for j = 1:size(clvl,2)
                    legendstr{i} = sprintf('%s%s x ',legendstr{i},allticks{j}{clvl(i,j)});
                end
                legendstr{i} = legendstr{i}(1:end-3);
            end
            %legend(gca,legendstr{:},'Location','BestOutside');
            legend(gca,legendstr{:},'Location','Best','interpreter','none');
        end
    end
    
    % set axes font size
    set(gca,'FontSize',12);
    
    % save the figure
    if istrue(cfg.save)
        cfg = km_setcfg(cfg,'dirreport');
        dir_report = getsubjsubdir(cfg,'Group','report');
        % create file name and save
        figname = sprintf('%s%s_%s',sessstr,namestr,paramstr);
        if strcmp(figname(1),'_'), figname = figname(2:end); end
        figname = fullfile(dir_report,figname);
        try
            if ~isfield(cfg.plot,'export'), cfg.plot.export = {'-pdf'}; end
            for j = 1:length(cfg.plot.export)
                if iscell(cfg.plot.export{j})
                    export_fig(figname,cfg.plot.export{j}{:});
                else
                    export_fig(figname,cfg.plot.export{j});
                end
            end
        catch
            saveas(h, figname,'jpg');
        end
    end
end