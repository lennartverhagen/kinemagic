function km_plot(cfg)
%--------------------------------------------------------------------------
% plot kinematic parameters as line or bar plots. The plotting of kinematic
% time-series is not yet supported.
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% set configuration
task = km_settask(mfilename);
cfg = km_setcfg(cfg,task);

% return if requested
if isfalse(cfg.(task))
    warning('FTW:Return','FTW %s: Nothing to do...',task);
    return
end


%% Get data for plotting
%----------------------------------------
if istrue(cfg.load)
    % Load data ready to plot
    %[cfg,param,tseries] = km_plot_loaddata(cfg);
    
else
    % Load raw dataset(s)
    [cfg,Sparam,SCI,Sepe,Stseries] = km_plot_getdata(cfg);
    
    % Construct contrasts
    %[cfg,Stseries] = km_plot_getcontr(cfg,Stseries);
    
    % Get descriptives
    [cfg,param,tseries] = km_plot_descriptives(cfg,Sparam,SCI,Sepe,Stseries);
    clear Sparam Stseries;
    
    % Save data
    %km_save_plotdata(cfg,param,tseries);
end


%% Reporting for statistics
%----------------------------------------
km_paramreport(cfg,param);


%% Plotting
%----------------------------------------

% parameter plotting
km_paramplot(cfg,param);

% time-series plotting
km_tseriesplot(cfg,tseries);

%--------------------------------------------------------------------------
