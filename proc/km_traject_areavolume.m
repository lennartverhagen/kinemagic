function varargout = km_traject_areavolume(val,k,flg)

% sort input
if nargin < 3
    flg = 'area';
end

dims = size(k);
nlvl = dims(1);
nrpt = dims(2);

% number of subjects
if length(dims)==3
    nsubj = dims(3);
else
    nsubj = 1;
end

% initialize
Marea = nan(nlvl,nrpt,nsubj);
if strcmpi(flg,'volume')
    Mvolume = nan(nlvl,nrpt,nsubj);
end
% loop over subjects
for s = 1:nsubj
    % loop over levels
    for i = 1:nlvl
        % loop over repetitions
        for j = 1:nrpt
            % calculate scaled eigenvalues
            abc = k(i,j,s) * sqrt(diag(val(:,:,i,j,s)));
            % store area
            Marea(i,j,s) = pi * prod(abc(2:3)) * 100^2; % bring from m^2 to cm^2
            if strcmpi(flg,'volume')
                % store volume
                Mvolume(i,j,s) = 4/3 * pi * prod(abc) * 100^3; % bring from m^3 to cm^3
            end
        end % nrpt
    end % nlvl
end % nsubj

% sort output
if strcmpi(flg,'volume')
    varargout = {Marea, Mvolume};
else
    varargout = {Marea};
end