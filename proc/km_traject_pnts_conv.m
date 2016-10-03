function pnts_conv = km_traject_pnts_conv(vec,val,avg,k,pnts)
% convolve a circle (pnts) with the principle axis of the ellipse

dims = size(k);
nlvl = dims(1);
nrpt = dims(2);
nax  = size(pnts,2);

% number of subjects
if length(dims)==3
    nsubj = dims(3);
else
    nsubj = 1;
end

% initialize
pnts_conv = nan(size(pnts,1),nax,nlvl,nrpt,nsubj);
% loop over subjects
for s = 1:nsubj
    % loop over levels
    for i = 1:nlvl
        % loop over repetitions
        for j = 1:nrpt
            % convolve circle with Eigen variates
            pnts_ij = k(i,j,s)*pnts*sqrt(val(:,:,i,j,s))*vec(:,:,i,j,s)';
            % translate origin
            pnts_conv(:,:,i,j,s) = bsxfun(@plus,avg(:,i,j,s)',pnts_ij);
        end
    end
end