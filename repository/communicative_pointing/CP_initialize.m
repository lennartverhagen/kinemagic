function dir_CP = CP_initialize
%--------------------------------------------------------------------------
% This function inilializes the directories and 
% adds paths (see CP_startup) for CP analysis.
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------


%% initialize
% check if 'CommPoint' is on top of the path (if not: run startup)
path_str = regexp(path,':','split');
if isempty(strfind(pwd,'CommPoint')) || ...
    isempty(strfind(path_str{1},'CommPoint')) || ...
    all(cellfun(@isempty,strfind(path_str,'kinemagic')))

    % get home directory
    if ispc
        username = 'ankethuis';
        dir_home = fullfile(filesep,'CommPointing','Data-analysis');
    else
        username = getenv('USER');
        if strcmpi(username,'lenver')
            dir_home = fullfile(filesep,'home','action','ankmur');
        else
            dir_home = fullfile(filesep,'home','action',username);
        end
    end
    
    % get and go to CommPoint directory
    if strcmpi(username,'lenver') || strcmpi(username,'ankmur')
        %dir_CP = fullfile(dir_home,'Data','CommPoint','KIN');
    %elseif strcmpi(username,'ankmur')
        dir_CP = fullfile(dir_home,'CommPointing','Data-analysis','KIN');
    elseif strcmpi(username,'ankethuis')
        dir_CP = fullfile(dir_home,'KIN'); % please adjust
    else
        dir_CP = fullfile(dir_home,'does not exist');
    end
    cd(fullfile(dir_CP,'mfiles'));
    
    % initialize matlab
    CP_startup('kinemagic');
    cd(dir_CP);
else
    % hard coded
    username = getenv('USER');
    if strcmpi(username,'lenver')
        dir_home = fullfile(filesep,'home','action','ankmur');
    else
        dir_home = fullfile(filesep,'home','action',username);
    end
    dir_CP = fullfile(dir_home,'CommPointing','Data-analysis','KIN');
    % soft coded (but be aware in what folder you currently are)
    %dir_CP = pwd;
    %dir_CP = strrep(dir_CP,'/mfiles','');
end