function [cfg,data] = km_readexpfun_GraspingSemantics(cfg,data)
%--------------------------------------------------------------------------
% KM_READEXPFUN_GRASPINGSEMANTICS fills in experimentally specific fields in
% the cfg and data structures after the data is read in
%
% See also KM_READ
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------


% you need the subject number to know in which category the subject belongs
subjnr = str2double(regexp(cfg.subj,'\d*','match','once'));    

% number of markers
d = finddim(data.dimord,'marker');
nmarker = size(data.pos{1},d);
if nmarker ~= 4
    error('the number of markers needs to be exactly 4, for the GraspingSemantics experiment')
end

% marker labels
% w: wrist
% j: metacarpophalangeal joint
% t: thumb
% i: index
if ismember(subjnr,[1:20 99])
    data.label = {'w','j','t','i'};
else
    error('this subject (%d) is not supported',subjnr);
end

% make time relative to start of run
for r = 1:length(data.time)
    data.time{r} = data.time{r} - data.time{r}(1);
end

% swap axis
if ismember(subjnr,[1 2 3 99])
    data.axis = {'y','x','z'};
    
    % flip y axis
    idx = ismember(data.axis,'y');
    for p = 1:length(data.pos)
        data.pos{p}(:,idx,:) = -data.pos{p}(:,idx,:);
    end
else
    data.axis = {'x','y','z'};
end
    