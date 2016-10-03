function [cfg,data] = km_readexpfun_GraspSlant(cfg,data)
%--------------------------------------------------------------------------
% KM_READEXPFUN_GRASPSLANT fills in experimentally specific fields in
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
    error('the number of markers needs to be exactly 4, for the GraspSlant experiment')
end

% marker labels
% w: wrist
% mcp: metacarpophalangeal joint
% t: thumb
% i: index
% obj: object
if subjnr == 7
    data.label = {'w','t','i','j'};
elseif ismember(subjnr,8:10)
    data.label = {'j','t','i','o'};
elseif ismember(subjnr,11:35)
    data.label = {'t','i','j','o'};
else
    error('this subject (%d) is not supported',subjnr);
end

% make time relative to start of run
for r = 1:length(data.time)
    data.time{r} = data.time{r} - data.time{r}(1);
end
