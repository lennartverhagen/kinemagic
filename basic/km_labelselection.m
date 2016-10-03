function label = km_labelselection(label,datalabel)

% try to use the standard version supplied by FieldTrip, otherwise use the
% script below.
funname   = 'ft_channelselection';
if exist(funname,'file') == 2
    funhandle = str2func(funname);
    label = funhandle(label, datalabel);
    return
end


%% Rip off of FieldTrip function
%----------------------------------------
if any(size(label) == 0)
  % there is nothing to do if it is empty
  return
end

if isnumeric(label)
  % change index into labelname
  label = datalabel(label);
  return
end

if ~iscell(label)
  % ensure that a single input argument like 'all' also works
  label = {label};
end

if ~iscell(datalabel)
  % ensure that a single input argument like 'all' also works
  datalabel = {datalabel};
end

% ensure that both inputs are column vectors
label     = label(:);
datalabel = datalabel(:);

% remove labels that occur more than once, this sorts the labels alphabetically
[label, indx] = unique(label);
% undo the sorting, make the order identical to that of the data labels
[~, indx] = sort(indx);
label = label(indx);

[~, labelindx] = match_str(datalabel, label);
if length(labelindx)==length(label)
  % there is a perfect match between the labels and the datalabels, only some reordering is needed
  label = label(labelindx);
  % no need to look at label groups
  return
end

% define the known groups with label labels
labelall  = datalabel;

% figure out if there are bad labels or label groups that should be excluded
findbadlabel = strncmp('-', label, length('-'));      % bad labels start with '-'
badlabel = label(findbadlabel);
if ~isempty(badlabel)
  for i=1:length(badlabel)
    badlabel{i} = badlabel{i}(2:end);                 % remove the '-' from the label label
  end
  badlabel = km_labelselection(badlabel, datalabel); % support exclusion of label groups
  label(findbadlabel) = [];                           % remove them from the labels to be processed
end

% use regular expressions to deal with the wildcards
labelreg = false(size(datalabel));
findreg = [];
for i=1:length(label)
  if length(label{i})>1 && label{i}(1)=='*'
    % the wildcard is at the start
    labelreg = labelreg | ~cellfun(@isempty, regexp(datalabel, ['.*' label{i}(2:end) '$'], 'once'));
    findreg  = [findreg i];
  elseif length(label{i})>1 && label{i}(end)=='*'
    % the wildcard is at the end
    labelreg = labelreg | ~cellfun(@isempty, regexp(datalabel, ['^' label{i}(1:end-1) '.*'], 'once'));
    findreg  = [findreg i];
  elseif length(label{i})>1 && any(label{i}=='*')
    % the wildcard is in the middle
    sel  = strfind(label{i}, '*');
    str1 = label{i}(1:(sel-1));
    str2 = label{i}((sel+1):end);
    labelreg = labelreg | ~cellfun(@isempty, regexp(datalabel, ['^' str1 '.*' str2 '$'], 'once'));
    findreg  = [findreg i];
  end
end
labelreg = datalabel(labelreg);

% determine if any of the known groups is mentioned in the label list
findall        = find(strcmp(label, 'all'));

% remove any occurences of groups in the label list
label([findall;findreg]) = [];

% add the full channel labels to the channel list
if findall,        label = [label; labelall]; end
if findreg,        label = [label; labelreg]; end

% remove label labels that have been excluded by the user
badindx = match_str(label, badlabel);
label(badindx) = [];

% remove label labels that are not present in the data
labelindx = match_str(label, datalabel);
label  = label(labelindx);

% remove labels that occur more than once, this sorts the labels alphabetically
label = unique(label);

% undo the sorting, make the order identical to that of the data labels
[~, indx] = match_str(datalabel, label);
label = label(indx);