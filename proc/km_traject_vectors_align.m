function [vec, val] = km_traject_vectors_align(vec,val,gra)
% restructure eigvec and eigval of each ellipsoid to match its neighbour,
% starting with the middle point

dims = size(vec);
% add dimension for levels if not included
if length(dims) == 3
    dims = [dims(1:2) 1 dims(3)];
    vec = reshape(vec,dims);
    val = reshape(val,dims);
end

% number of levels
nlvl = dims(3);
% number of repetitions (and reference point)
nrpt = dims(4);
hnrpt = round(nrpt/2);

% number of subjects
if length(dims)==5
    nsubj = dims(5);
else
    nsubj = 1;
end

% loop over subjects
for s = 1:nsubj
    % loop over all levels
    for i = 1:nlvl
        % use eigenvectors of the middle point of subject 1 as a reference
        tmp0_vec = vec(:,:,i,hnrpt,1);
        
        % determine vector closest to the direction of the average trajectory and
        % place it up front
        [~,m] = max(abs(vectangle(gra(:,i,hnrpt),tmp0_vec)));
        tmp0_vec = tmp0_vec(:,[m setdiff(1:3,m)]);
        
        % loop over all repetitions, starting from the middle
        for j = [hnrpt:nrpt (hnrpt-1):-1:1]
            % when doing the second half use i=hnrpt not i==nrpt as reference
            if j == (hnrpt-1)
                tmp0_vec = vec(:,:,i,hnrpt,s);
            end
            % retrieve current eigenvectors
            tmp_vec = vec(:,:,i,j,s);
            [tmp_vec,ni] = matchvectors(tmp_vec,tmp0_vec,gra(:,i,j,s),3);
            % store for next iteration
            vec(:,:,i,j,s) = tmp_vec;
            tmp0_vec = tmp_vec;
            % retrieve eigenvalues
            tmp_val = val(:,:,i,j,s);
            % and place them in the same order as the eigenvectors
            val(1,1,i,j,s) = tmp_val(ni(1),ni(1));
            val(2,2,i,j,s) = tmp_val(ni(2),ni(2));
            val(3,3,i,j,s) = tmp_val(ni(3),ni(3));
        end
    end
end