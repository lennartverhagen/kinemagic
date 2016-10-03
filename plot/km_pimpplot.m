function varargout = km_pimpplot(varargin)
%--------------------------------------------------------------------------
% KM_PIMPPLOT checks if the file ft_pimpplot can be used. If not, it uses a
% copy of the code kept here. But this code copy is not regularly checked
% nor updated...
%
% See also KM_PLOT
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

funname = 'ft_pimpplot';
funhandle = str2func(funname);
if exist(funname,'file')==2
    [varargout{1:nargout}] = funhandle(varargin{:});
    return
else
    fundir = fullfile(filesep,'home','common','matlab','eeg-tms','fieldtripfork');
    funfile = fullfile(fundir,[funname '.m']);
    if exist(funfile,'file')==2
        tmpdir = pwd;
        cd(fundir);
        [varargout{1:nargout}] = funhandle(varargin{:});
        cd(tmpdir);
        return
    end
end


%--------------------------------------------------------------------------
% FT_PIMPPLOT can pimp lineplots, time frequency and topographical
% representations created by FieldTrip. Basically, it replots the plots
% created by imagesc and surface using contourf.
%
% See also FT_PLOT
%
% This file is a FieldTrip fork as part of the FieldTripWrapper toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% set input to ft_pimpplot
if nargin > 0,      hfig = varargin{1}; end
if nargin > 1,      cmap = varargin{2}; end
if nargin > 2,      pimp = varargin{2}; end

% check input
if nargin < 3,      pimp = true;    end
if isempty(pimp),	pimp = true;	end
if nargin < 2,      cmap = @jet;    end
if isempty(cmap),	cmap = @jet;	end
if (islogical(cmap) || ischar(cmap)) && numel(cmap)==1
    pimp = cmap;    cmap = @jet;
end
if nargin < 1,      hfig = gcf;     end

% set output
varargout = {};
if nargout > 0
    varargout = {hfig};
end
    
% return if requested
if isfalse(pimp)
    return
end

% specify number of steps in colormap
stps = 2^9;

% create colormap if needed
if isa(cmap,'function_handle')
    cmap = cmap(stps);
end

% you might want to make patches and image layers non-transparant to be
% able to export them properly to pdf and eps (where you could set the
% transparency again, if you would like)
nontransparantflg = true;

% get plot and axes handles
hs = findobj(hfig,'Type','image','-and','Tag','cip');       % single channel
ht = findobj(hfig,'Type','surface');                        % topoplot
hcb = findobj(hfig,'Tag','Colorbar');                       % colorbar axis
hla = get(findobj(hfig,'Type','line'),'Parent');            % lineplot axis
if iscell(hla), hla = unique([hla{:}]);	end
hta = get(ht,'Parent');                                     % topoplot axis
if iscell(hta), hta = unique([hta{:}]);	end
hla = setdiff(hla,hta);                                     % lineplot but not topoplot axis
hsa = get(hs,'Parent');                                     % single channel axis
if iscell(hsa), hsa = unique([hsa{:}]);	end
hla = setdiff(hla,hsa);                                     % lineplot but not topoplot axis
%ha = findobj(hfig,'Type','axes');                          % all axis

% pimp the single channel time-frequency representation
for i = 1:length(hs)
    
    % get title and data
    %title_str = get(get(get(hs(i),'Parent'),'Title'),'String');
    X = get(hs(i),'XData');
    Y = get(hs(i),'YData');
    Z = get(hs(i),'CData');
    A = get(hs(i),'AlphaData');
    %c = caxis(get(hs(i),'Parent'));
    
    % give warning in case of true color mode
    if ndims(Z) == 3
        % FIXME: export_fig can not handle transparent layers, therefore
        % you might want to lighten your data using true color mode.
        % Unfortunately this is not supported yet by ft_pimpplot
        warning('FT:PimpPlot','True color mode (size(C)==[m n 3]) is not yet supported.')
        return
    end
    
    % expand matrix. This is necessary because the original plot was made
    % using image
    x_step = mode(abs(diff(X)));
    y_step = mode(abs(diff(Y)));
    X = [X(1)-x_step X X(end)+x_step];
    Y = [Y(1)-y_step Y Y(end)+y_step];
    Z = [Z(:,1) Z Z(:,end)];
    Z = [Z(1,:); Z; Z(end,:)];
    A = [A(:,1) A A(:,end)];
    A = [A(1,:); A; A(end,:)];
    
    % create new figure if requested
    if pimp < 0
        if ~exist('hnfig','var')
            hnfig = copyobj(hfig,0);
        end
        figure(hfig);
    end
    
    % re-plot data using contourf
    set(hfig,'CurrentAxes',get(hs(i),'Parent'));
    delete(hs);
    hold on;
    [~,hn] = contourf(X,Y,Z,stps,'LineStyle','none');
    colormap(cmap);
    set(gca,'Box','off');
    axis([min(X)+x_step/2 max(X)-x_step/2 min(Y)+y_step/2 max(Y)-y_step/2]);
    % axis tight;
    
    % add the alpha-map from the original image
    if ~all(A(:))
        whitemat = ones([size(A) 3]);
        A(A>=1) = 0;    % A(A> 0 & A<=0.5) = 0.75;
        image(X,Y,whitemat,'AlphaData',A);
    end
    
    % restack the objects
    hl = findobj(hfig,'Type','line');
    uistack(hl,'top');
    %uistack(hn,'bottom');
    
    % make the axis lines thicker and ticks markers smaller (not the labels)
    set(get(hn,'Parent'),'LineWidth',1.5);
    set(get(hn,'Parent'),'TickLength',[0.003; 0.025]);
    set(get(hn,'Parent'),'TickDir','out');
    
    % % add line at time = 0
    % hold on; line([0 0],[min(Y) max(Y)],'Color','k','LineWidth',2);
    hold off;
    
end

% pimp the topographical plot
for i = 1:length(ht)
    
    % get data
    X = get(ht(i),'XData');
    Y = get(ht(i),'YData');
    Z = get(ht(i),'CData');
    
    % adjust X and Y to have centre at zero. This is necessary because the
    % original plot was made using surface
    X = X + 0.5/size(X,1);
    Y = Y + 0.5/size(Y,1);
    
    % create new figure if requested
    if pimp < 0
        if ~exist('hnfig','var')
            hnfig = copyobj(hfig,0);
        end
        figure(hfig);
    end
        
    % re-plot data using contourf
    set(hfig,'CurrentAxes',get(ht(i),'Parent'));
    delete(ht(i));
    hold on;
    [~,hn] = contourf(X,Y,Z,stps,'LineStyle','none');
    colormap(cmap);
    
    % restack the objects
    uistack(hn,'bottom');
    
end

% reset the tick markers of the colorbar
for i = 1:length(hcb)
    set(hcb(i),'TickLength',[0.003; 0.025]);
    set(hcb(i),'TickDir','out');
end

% in the lineplots make the axis lines thicker and ticks markers smaller (not the labels)
for i = 1:length(hla)
    lnwidth = get(hla(i),'LineWidth');
    set(hla(i),'LineWidth',max(lnwidth,1.5));
    if ~strcmpi(get(get(hla(i),'parent'),'Type'),'figure')
        hparent = get(hla(i),'parent');
    else
        hparent = hla(i);
    end
    axwidth = get(hparent,'LineWidth');
    set(hparent,'LineWidth',max(axwidth,1.5));
    set(hparent,'TickLength',[0.003; 0.025]);
    set(hparent,'TickDir','out');
    if nontransparantflg
        % find the children and make the error clouds non-transparant
        hp = findobj(hla(i),'Type','patch');	% patch objects in lineplot axis
        for j = length(hp):-1:1
            A = get(hp(j),'FaceAlpha');
            C = get(hp(j),'FaceColor');
            C = C + (1-C)*(1-A);
            set(hp(j),'FaceAlpha',1);
            set(hp(j),'FaceColor',C);
            % restack the patches
            uistack(hp(j),'bottom');
        end
        % restack the lines
        hl = findobj(hla(i),'Type','line');     % line objects in lineplot axis
        for j = length(hl):-1:1
            uistack(hl(j),'top');
        end
    end
end


% % continue working on new figure if requested
% if pimp < 0
%     hfig = hnfig;
% end

% remove the title
title('');
    
% make background white and remove box
set(hfig,'Color','w');
set(gca,'Box','off');

% set the renderer to OpenGL
set(gcf,'Renderer','OpenGL')


% set correct printing properties
set(hfig,'PaperType','<custom>');
%set(h,'PaperUnits','centimeters');
if strcmpi(get(hfig,'PaperUnits'),'centimeters')
    set(hfig,'PaperSize',[20 15]);
    set(hfig,'PaperPosition',[0 0 20 15]);
elseif strcmpi(get(hfig,'PaperUnits'),'inches')
    set(hfig,'PaperSize',[8 6]);
    set(hfig,'PaperPosition',[0 0 8 6]);
end

% set output
if nargout > 0
    varargout = {hfig};
end


