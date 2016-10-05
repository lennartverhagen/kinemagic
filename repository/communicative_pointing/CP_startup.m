function CP_startup(varargin)
%--------------------------------------------------------------------------
% This function sets the path and other workspace variable at matlab startup
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------

% store current path to add at end of script
dir_caller = pwd;

% Set path to default
warning('off','MATLAB:dispatcher:pathWarning');
%path(pathdef);
restoredefaultpath;
warning('on','MATLAB:dispatcher:pathWarning');

% clear command line
clc;

% report
disp('Initializing for CommPoint analyses');

% Sort input
kinemagic = true;
reset = false;
for v = 1:length(varargin)
    if any(strcmpi(varargin{v},{'km','kinemagic','kinematic'}))
        kinemagic = true;
    elseif isempty(varargin{v}) || any(strcmpi(varargin{v},{'reset','clean','empty'}))
        reset = true;
    elseif ischar(varargin{v})
        fprintf(1,'Input ''%s'' is not supported.\n',varargin{v});
    else
        disp('Input not recognized.');
    end
end

% reset all flgs if requested
if reset
    kinemagic = false;
end

% find home directory and common-matlab directory
if ispc
    % override default home directory when on PC and use central storage
    % instead
    username = getenv('USERNAME');
    % home_dir = getenv('USERPROFILE');
    % dir_home = 'M:\';
    % dir_common = fullfile('N:','matlab');
    dir_home = [filesep filesep fullfile('fileserver','home','action',username)];
    dir_common = [filesep filesep fullfile('fileserver','home','common','matab')];
else
    dir_home = getenv('HOME');
    dir_common = fullfile(filesep,'home','common','matlab');
end

% Add kinematic analysis program to path
if kinemagic
    disp('Setting up path for Kinemagic');
    kp = genpath(fullfile(dir_common,'kinemagic'));
    % do not inlcude the private folder
    pp = genpath(fullfile(dir_common,'kinemagic','private'));
    kp = strrep(kp,pp,'');
    % add path
    addpath(kp);
end

% Add the qsub function to allow distributed processing jobs
addpath(fullfile(dir_common,'fieldtrip','qsub'));

% Add root paths to top of path
addpath(dir_home);

% Save path and go to root
cd(dir_home);
savepath(dir_home);

% Add matlab path after saving to prevent warning on startup
addpath(fullfile(dir_home,'matlab'));

% Add original calling path and go to directory
addpath(dir_caller);
cd(dir_caller)

% figure color defaults
% This is a reduced version of the command: colordef('white');
whitebg(0,[1 1 1]);
set(0,'DefaultFigureColor',[1 1 1])
%set(0,'DefaultAxesColor',[0 0 0])
set(0,'DefaultAxesColor','none')
set(0,'DefaultAxesColorOrder',[1 0 0;0 .5 0;0 0 1;0 .75 .75;.75 0 .75;.75 .75 0;.25 .25 .25]) % rgbcmyk
set(0,'DefaultFigureColormap',jet(64))
set(0,'DefaultSurfaceEdgeColor',[0 0 0])

% dock figures by default
set(0,'DefaultFigureWindowStyle','docked')