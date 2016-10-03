function km_tseriesplot(cfg,alldata)
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
if ~isfield(cfg,task) && ~isempty(alldata) && ~isempty(fieldnames(alldata))
    cfg.(task) = 'yes';
end
cfg = km_setcfg(cfg,task,'strict');

% return if requested
if isfalse(cfg.(task)) || isempty(alldata) || isempty(fieldnames(alldata))
    return
end


% select conditions and factors
fact    = cfg.plot.fact;
fname	= cfg.plot.fname;
ntseries= length(cfg.plot.tseries);
nfact	= length(fname);
%ncond	= length(cname);
nsubj   = length(cfg.subj);

% loop over time-series
for ts = 1:ntseries
    
    % the time is normalized
    tim = 0:100;
        
    % time-series name
    tsname = cfg.plot.tseries{ts};
    tsdata = alldata.(tsname);
    
    % draw sensors in plot or not
    if ~isfalse(cfg.plot.drawsens)
        drawsens = cfg.plot.drawsens(ts);
        if ~isfalse(drawsens)
            if iscell(drawsens.marker)
                %[~,drawsens.marker] = ismember(drawsens.marker,tsdata.marker);
                drawsens.marker = km_labelselection(drawsens.marker,tsdata.marker);
                drawsens.marker = match_str(tsdata.marker,drawsens.marker);
            end
            stim = drawsens.tpoint;
            if isempty(drawsens.tpoint)
                stim = linspace(tim(1),tim(end),drawsens.nmark);
            else
                drawsens.nmark = length(stim);
            end
        end
    else
        drawsens = false;
    end
    
        
    % loop over factors for parameters
    for f = 1:nfact
                
        % get plotting specifics
        spec.lncl       = cfg.plot.flcl{f};
        spec.lnstyle	= cfg.plot.lnstyle{f};
        spec.lnwidth	= cfg.plot.lnwidth{f};
        spec.marker     = cfg.plot.marker{f};
        
        % get data descriptives
        dat = tsdata.(fname{f}).mean;
        dimord = tsdata.dimord;
        axs = tsdata.axis;
        mrkr = km_labelselection(cfg.tseries.(tsname).marker,tsdata.marker);
        mrkr = match_str(tsdata.marker,mrkr);
        %[~,mrkr] = ismember(cfg.tseries.(tsname).marker,tsdata.marker);
        
        % get levels
        lvl = cell(1,length(fact{f}));
        nlvl = nan(size(lvl));
        for ff = 1:length(fact{f})
            lvl{ff} = cfg.plot.lvl.(fact{f}{ff})(:)';
            nlvl(ff) = length(lvl{ff});
        end
        clvl = arraycomb(lvl{1:end-1});
        
        % get error
        if nsubj == 1
            err = tsdata.(fname{f}).sem;
        elseif ((length(nlvl) > 1 && any(nlvl(1:end-1) == 2)) || ...
                (length(nlvl) == 1 && any(nlvl == 2)) ) && ...
                isfield(tsdata.(fname{f}),'diff')
            err = tsdata.(fname{f}).diff.sem;
            dim = find(nlvl == 2,1,'first');
            tile = ones(1,length(nlvl));
            tile(dim) = 2;
            if length(tile) == 1,   tile = [1 tile];    end
            err = repmat(err,tile);
        else
            err = tsdata.(fname{f}).sem;
        end
        
        % % reshape dat to n x m
        %dims = size(dat);
        %dims = [prod(dims(1:end-1)) dims(end)];
        %dat = reshape(dat,dims);
        
        % % reshape err to fit in one figure
        % err = reshape(err,dims);
        
        % initiate new figure
        h = figure;
        
        % loop over cells
        for i = 1:numel(dat)
            % plotting specifics
            if ischar(spec.lncl)
                lncl  = spec.lncl(i);
            elseif size(spec.lncl,2)==3
                lncl  = spec.lncl(i,:);
            end
            if ~isfalse(drawsens)
                lnstyle	= 'none';
                lnwidth	= 1;
                if ~isfalse(drawsens.connect)
                    ds_lnstyle	= spec.lnstyle{i};
                    ds_lnwidth	= spec.lnwidth(i);
                else
                    ds_lnstyle	= 'none';
                    ds_lnwidth	= 1;
                end    
            else
                lnstyle	= spec.lnstyle{i};
                lnwidth	= spec.lnwidth(i);
            end
            %marker	= spec.marker{i};   % the plot will be too cramped if
            %you use a marker to specify the level of the third factor
            marker	= 'no';
            
            switch lower(dimord)
                case {'sensor_axis_time','marker_axis_time'}
                    % loop over sensors
                    for m = mrkr
                        if length(axs) > 2
                            x = squeeze(dat{i}(m,1,:));
                            y = squeeze(dat{i}(m,2,:));
                            z = squeeze(dat{i}(m,3,:));
                            plot3(x,y,z,'Color',lncl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                            % I would not know how to plot an 3D error
                            % cloud...
                            % xerr = squeeze(err{i}(m,1,:));
                            % yerr = squeeze(err{i}(m,2,:));
                            % zerr = squeeze(err{i}(m,3,:));
                            % plot3errorcloud(xerr,yerr,zerr,lncl);
                        elseif length(axs) == 2
                            x = squeeze(dat{i}(m,1,:));
                            y = squeeze(dat{i}(m,2,:));
                            plot(x,y,'Color',lncl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                            % I would not know how to plot an 2D error
                            % cloud...
                            % xerr = squeeze(err{i}(m,1,:));
                            % yerr = squeeze(err{i}(m,2,:));
                            % plot2errorcloud(xerr,yerr,lncl);
                        else
                            y = squeeze(dat{i}(m,1,:))';
                            plot(tim,y,'Color',lncl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                            yerr = squeeze(err{i}(m,1,:))';
                            ploterrorcloud(tim,y,yerr,lncl);
                        end
                        hold on
                    end
                    % draw sensors (and connect them)
                    if ~isfalse(drawsens)
                        sdat = interpdat(dat{i}(drawsens.marker,:,:),tim,stim);
                        for n = 1:drawsens.nmark
                            if length(axs) > 2
                                x = squeeze(sdat(:,1,n));
                                y = squeeze(sdat(:,2,n));
                                z = squeeze(sdat(:,3,n));
                                plot3(x,y,z,'Color',lncl,'Marker','o','LineStyle',ds_lnstyle,'LineWidth',ds_lnwidth);
                            elseif length(axs) == 2
                                x = squeeze(sdat(:,1,n));
                                y = squeeze(sdat(:,2,n));
                                plot(x,y,'Color',lncl,'Marker','o','LineStyle',ds_lnstyle,'LineWidth',ds_lnwidth);
                            else
                                y = squeeze(sdat(:,1,n));
                                plot(stim,y,'Color',lncl,'Marker','o','LineStyle',ds_lnstyle,'LineWidth',ds_lnwidth);
                            end
                            hold on
                        end
                    end
                case 'axis_time'
                    if length(axs) > 2
                        plot3(dat{i}(1,:),dat{i}(2,:),dat{i}(3,:),'Color',lncl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                    elseif length(axs) == 2
                        plot(dat{i}(1,:),dat{i}(2,:),'Color',lncl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                    else
                        plot(tim,dat{i}(1,:),'Color',lncl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                        yerr = squeeze(err{i}(1,:));
                        ploterrorcloud(tim,dat{i}(1,:),yerr,lncl);
                    end
                case 'time'
                    plot(tim,dat{i},'Color',lncl,'Marker',marker,'LineStyle',lnstyle,'LineWidth',lnwidth);
                otherwise
                    warning('KM:TSeriesPlot','This dimension order [%s] is not supported',dimord);
            end
            hold on
        end
        hold off
        
        % adjust axis
        if length(axs) > 1,  axis equal; end
        
        % % adjust y limits
        % yticks = getaxrange(dat,err);
        % ylim([min(yticks) max(yticks)]);
        
        % % set y-ticks
        % %set(gca,'YTickMode','auto')
        % set(gca,'YTick',yticks);
        
        % pimp plot
        km_pimpplot(h)
                 
        % wrap up the figure
        fig_wrapup(cfg,h,tsname,fname{f},axs,tsdata.marker(mrkr));
        
    end
    
end
%--------------------------------------------------------------------------

%% function interpdat
%----------------------------------------
function dat = interpdat(dat,x,xi)

% number of samples
dims = size(dat);
% original data
y = permute(dat,[3 2 1]);
y = reshape(y,dims([3 2 1]));
% interpolate
yi = interp1(x,y,xi,'cubic','extrap');
% shape back to original size
dims(3) = size(yi,1);
yi = permute(yi,[3 2 1]);
dat = reshape(yi,dims);


%% function fig_wrapup
%----------------------------------------
function fig_wrapup(cfg,h,tseriesstr,namestr,axs,marker)
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
    figname = sprintf('cond:%s  tseries:%s  subj:%s  sess:%s',spellout(namestr,1),tseriesstr,subjstr,sessstr);
    set(h,'Color','w'); set(gca,'Box','off');
    set(h,'NumberTitle','off','Name',figname);
    
    % set the y-axis label
    switch lower(tseriesstr)
        case {'pos','gripapp'}
            unitstr = 'm';
        case {'vel','gripvel'}
            unitstr = 'm/s';
        case {'acc'};
            unitstr = 'm/s²';
        case {'ori','gripori'};
            unitstr = 'deg';
        otherwise
            unitstr = 'a.u.';
    end
    if length(axs) > 2
        xstr = sprintf('%s-%s (%s)',tseriesstr,axs{1},unitstr);
        ystr = sprintf('%s-%s (%s)',tseriesstr,axs{2},unitstr);
        zstr = sprintf('%s-%s (%s)',tseriesstr,axs{3},unitstr);
        set(get(gca,'XLabel'),'String',xstr);
        set(get(gca,'YLabel'),'String',ystr);
        set(get(gca,'ZLabel'),'String',zstr);
    elseif length(axs) == 2
        xstr = sprintf('%s-%s (%s)',tseriesstr,axs{1},unitstr);
        ystr = sprintf('%s-%s (%s)',tseriesstr,axs{2},unitstr);
        set(get(gca,'XLabel'),'String',xstr);
        set(get(gca,'YLabel'),'String',ystr);
    else
        xstr = sprintf('time (s)');
        ystr = sprintf('%s (%s)',tseriesstr,unitstr);
        set(get(gca,'XLabel'),'String',xstr);
        set(get(gca,'YLabel'),'String',ystr);
    end
    
    % tick labels and legend
    if isfield(cfg.plot,'lvltick')
        % get labels
        f = find(strcmpi(cfg.plot.fname,namestr));
        allticks = cell(1,length(cfg.plot.fact{f}));
        lvl = cell(1,length(cfg.plot.fact{f}));
        for ff = 1:length(cfg.plot.fact{f})
            fname = cfg.plot.fact{f}{ff};
            allticks{ff} = cfg.plot.lvltick.(fname);
            [~, lvl{ff}] = unique(cfg.plot.lvl.(fname)(:)');
        end
        % create a legend
        clvl = arraycomb(lvl{:});
        legendstr = cell(1,length(marker)*size(clvl,1));
        k = 1;
        for i = 1:size(clvl,1)
            for m = 1:length(marker)
                if length(marker) == 1
                    legendstr{k} = '';
                else
                    legendstr{k} = [marker{m} ' - '];
                end
                for j = 1:size(clvl,2)
                    legendstr{k} = sprintf('%s%s x ',legendstr{k},allticks{j}{clvl(i,j)});
                end
                legendstr{k} = legendstr{k}(1:end-3);
                k = k + 1;
            end
        end
        legend(gca,legendstr{:},'Location','Best');
    end
    
    % save the figure
    if istrue(cfg.save)
        cfg = km_setcfg(cfg,'dirreport');
        dir_report = getsubjsubdir(cfg,'Group','report');
        if ~exist(dir_report,'dir'), mkdir(dir_report); end
        % create file name and save
        figname = fullfile(dir_report,sprintf('TS_%s%s_%s',sessstr,namestr,tseriesstr));
        try
            if ~isfield(cfg.plot,'export'), cfg.plot.export = {'-pdf'}; end
            for j = 1:length(cfg.plot.export)
                if iscell(cfg.plot.export{j})
                    export_fig(figname,cfg.plot.export{j}{:});
                else
                    export_fig(figname,cfg.plot.export{j});
                end
            end
        catch errmsg
            % report on errmsg
            % FIXME: implement errmsg report
            saveas(h, figname,'jpg');
        end
    end
end