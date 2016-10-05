function [cfg,data] = CP_readfun(cfg,data)
%--------------------------------------------------------------------------
% This function fills in experimentally specific fields in
% the cfg and data structures after the data is read in.
%
% See also KM_READ
%
% This file is part of the CommPoint toolbox,
% an extension of the KineMagic toolbox
% Copyright (C) 2014, Anke Murillo Oosterwijk
% a.murillooosterwijk@donders.ru.nl
% version 1
%--------------------------------------------------------------------------


% you need the subject number to know in which category the subject belongs
subjnr = str2double(regexp(cfg.subj,'\d*','match','once'));    

% number of markers
d = finddim(data.dimord,'marker');
nmarker = size(data.pos{1},d);
if nmarker ~= 4
    error('the number of markers needs to be exactly 4, for the graspingTMSEEG experiment')
end

% ajust logfile association of runs, subj13
if subjnr == 13
    data.run = [1 1];
end

% marker labels
% i: index
% hi: hand at index, first metacarpophalangeal joint
% hp: hand at pink, fourth metacarpophalangeal joint
% w: wrist
data.label = {'i','hi','hp','w'};

% make time relative to start of run
for r = 1:length(data.time)
    data.time{r} = data.time{r} - data.time{r}(1);
end
