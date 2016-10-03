% initialize
dir_root	= fullfile(filesep,'home','action','lenver','Data','TMSEEG');
subjects	= {      '',      '',      '',      '',      '',...
                     '','subj07','subj08','subj09','subj10',...
               'subj11','subj12','subj13','subj14','subj15',...
               'subj16','subj17','subj18','subj19','subj20',...
               'subj21','subj22','subj23','subj24','subj25',...
               'subj26','subj27','subj28','subj29','subj30',...
               'subj31','subj32','subj33','subj34','subj35',...
               'subj00'};

% select subjects
subjsel = 8:35;

% define the action to be taken
action = 'find';
%action = 'delete';
%action = 'rename';
%action = 'copy';
%action = 'deletedir';
%action = 'finddir';

% define the filter string to select files
tmpdir = fullfile('kinematics','ana');
filtstr = '*.mat';
newstr = '*FILTlp15_ANA*.mat';

% type determination

if isempty(regexpi(action,'dir','once'))
    diraction = false;
    typestr = 'files';
else
    diraction = true;
    typestr = 'directories';
    action = strrep(action,'dir','');
end

% report initiation
disp(' ');
if strcmpi(action(end),'e')
    str = ['action: ' action(1:end-1) 'ing'];
else
    str = ['action: ' action 'ing'];
end
if diraction
    disp([str ' a directory']);
else
    disp(str);
    disp(['filt string: ' filtstr]);
end


% loop over subjects
A = cell(1,max(subjsel));
for s = subjsel
    
    A{s} = 0;
    
    workdir = fullfile(dir_root,subjects{s},tmpdir);
    
    % work on directories
    if diraction
        if exist(workdir,'dir')
            switch lower(action)
                case 'delete'
                    rmdir(workdir,'s');
                case 'find'
            end
            A{s} = A{s} + 1;
        end
        continue
    end
    
    % work on files
    fname = dir(fullfile(workdir,['*' filtstr '*']));
    for f = 1:length(fname)
        filtfile = fullfile(workdir,fname(f).name);
        switch lower(action)
            case {'rename','move'}
                newfile = fullfile(workdir,strrep(fname(f).name,filtstr,newstr));
                movefile(filtfile,newfile);
            case 'copy'
                newfile = fullfile(workdir,strrep(fname(f).name,filtstr,newstr));
                copyfile(filtfile,newfile);
            case 'delete'
                delete(filtfile);
            case 'find'
                if fname(f).bytes == 0
                    continue
                end
                % do nothing
            otherwise
                error('The action %s is not recognized',action)
        end
        A{s} = A{s} + 1;
    end
end

% report result
B = [A{:}];
if any(strcmpi(action,{'rename','move','copy'}))
    disp(['new string: ' newstr]);
end
if length(unique(B))==1
    disp(['all subjects have the same number of ' typestr]);
else
    disp(['not all subject have the same number of ' typestr]);
end
disp(num2str([subjsel; B]));
if max(B) > 0
    disp([num2str(floor((100*sum(B))/(max(B)*length(B)))) ' percent complete'])
end
