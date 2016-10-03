function pnts = km_traject_pnts_align(pnts,avg)
% set the origin of each ellipsoide to the point with highest x value

dims = size(pnts);
npnt = dims(1);
nlvl = dims(3);
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
    % loop over levels
    for i = 1:nlvl
        % use points of previous/first level of first subject as reference
        tmp_ref = pnts(:,:,max(i-1,1),hnrpt,1);
        % ignore last point (it overlaps with the first)
        tmp_ref = tmp_ref(1:end-1,:);
        % find highest x value and circshift
        [~,idx] = max(tmp_ref(:,1),[],1);
        tmp_ref = circshift(tmp_ref,[1-idx 0]);
        % add data point to make circular again
        tmp_ref(end+1,:) = tmp_ref(1,:);
        % make sure that the reference points are oriented clock-wise (when
        % looking along the trajectory)
        if sum(tmp_ref(round(npnt/4),2:3) * diag([-1 1])) > sum(tmp_ref(round(3*npnt/4),2:3) * diag([-1 1]))
            tmp_ref = flipud(tmp_ref);
        end
        % get average trajectory as origin
        orig = squeeze(avg(:,i,:,s));
        for j = [hnrpt:nrpt (hnrpt-1):-1:1]
            % when doing the second half use j=hnrpt not j==nrpt as reference
            if j == (hnrpt-1)
                tmp_ref = pnts(:,:,i,hnrpt,s);
            end
            % ignore last point (it overlaps with the first)
            tmp = pnts(1:end-1,:,i,j,s);
            % find highest x value and circshift
            [~,idx] = max(tmp(:,1),[],1);
            tmp = circshift(tmp,[1-idx 0]);
            % add data point to make circular
            tmp(end+1,:) = tmp(1,:);
            % subtract origin
            tmp_ref0 = bsxfun(@minus,tmp_ref,orig(:,j)');
            tmp0 = bsxfun(@minus,tmp,orig(:,j)');
            % % now calculate distance of all points to all points
            % %C = cellfun(@(x) euclid(ttmp0,x),num2cell(ttmp,2),'UniformOutput',false);
            % %[~,idx] = min(cat(2,C{:}),[],1);
            % compare distance and flipped distance
            d_orig = sum(euclid(tmp_ref0,tmp0).^2);
            d_flip = sum(euclid(tmp_ref0,flipud(tmp0)).^2);
            if d_flip < d_orig
                tmp = flipud(tmp);
            end
            % use as reference for next round
            tmp_ref = tmp;
            % store
            pnts(:,:,i,j,s) = tmp;
        end
    end
end