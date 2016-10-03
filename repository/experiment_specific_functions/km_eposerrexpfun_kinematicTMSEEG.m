function epe = km_eposerrexpfun_kinematicTMSEEG(epe)
%--------------------------------------------------------------------------
% KM_EPOSERREXPFUN_KINEMATICTMSEEG averages the end-point-error values over
% slant to create a vert-hori factor
%
% See also KM_EPOSERR
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% return if no slant factor present
if ~any(ismember(epe.bfact,'slant'))
    return
end

% loop over factors
nfact = length(epe.fact);
for f = 1:nfact
    
    % get factor(s)
    fact = epe.fact{f};
    fname = epe.fname{f};
    
    % continue if no slant present in this factor
    if ~any(ismember(fact,'slant'))
        continue
    end
    
    % verthori factor
    newfact = fact;
    newfact(ismember(fact,'slant')) = {'verthori'};
    newfname = regexprep(fname,'slant','verthori');
    
    % get levels
    factlvl = cell(1,length(fact));
    for ff = 1:length(fact)
        factlvl{ff} = epe.lvl.(fact{ff});
    end
    
    % get indeces
    idx_slant = ismember(fact,'slant');
    if length(idx_slant) == 1
        idx_slant = [idx_slant false];
    end
    idx_vert = epe.lvl.slant < 45;
    idx_hori = epe.lvl.slant > 45;
    
    % determine dimensions
    dims = cellfun(@(x) length(x),factlvl);
    if length(dims) == 1
        dims = [dims 1];
    end
    newdims = dims;
    newdims(idx_slant) = 2;
    
    % calculate the dimension to permute
    alldims = 1:length(dims);
    permdim = [find(idx_slant) alldims(~idx_slant)];
                
    % loop over markers
    tmpnew = struct;
    marker = fieldnames(epe.(fname));
    marker = marker(~ismember(marker,'lvl'));
    for nm = 1:length(marker)
        tmpmark = epe.(fname).(marker{nm});
        
        % loop over axes
        axisname = fieldnames(tmpmark);
        for an = 1:length(axisname)
            
            % loop over fields
            flds = fieldnames(tmpmark.(axisname{an}));
            paramflds = flds(~ismember(flds,{'eigvec','eigval','xyz'}));
            for fn = 1:length(paramflds)
                % get data
                tmp = tmpmark.(axisname{an}).(paramflds{fn});
                % permute and reshape data to bring the slant to the first dimension
                tmp = permute(tmp,permdim);
                tmp = reshape(tmp,dims(permdim));
                % calculate the vert and hori averages
                tmpavg = nan(2,size(tmp,2));
                tmpavg(1,:) = nanmean(tmp(idx_vert,:),1);
                tmpavg(2,:) = nanmean(tmp(idx_hori,:),1);
                % reshape and permute back to fit the new dimensions
                tmpavg = reshape(tmpavg,newdims(permdim));
                tmpavg = ipermute(tmpavg,permdim);
                % store in new field
                tmpnew.(marker{nm}).(axisname{an}).(paramflds{fn}) = tmpavg;
            end
            % FIXME: in principle the other fields could be averaged as
            % well... But personally, I think it would be best to calculate
            % new coordinates based on the eigenvalues and vectors at the
            % time of reporting
            cellflds = flds(ismember(flds,{'eigvec','eigval','xyz'}));
            for fn = 1:length(cellflds)
                % get data
                tmp = tmpmark.(axisname{an}).(cellflds{fn});
                % permute and reshape data to bring the slant to the first dimension
                tmp = permute(tmp,permdim);
                tmp = reshape(tmp,dims(permdim));
                % calculate the vert and hori averages
                tmpdims = size(tmp{1});
                tmp = cellfun(@(x,y) reshape(x,[1 y]),tmp,repmat({tmpdims},size(tmp)),'UniformOutput',false);
                tmpavg = cell(2,size(tmp,2));
                for i = 1:size(tmp,2)
                    tmpavg{1,i} = reshape(nanmean(vertcat(tmp{idx_vert,i}),1),tmpdims);
                    tmpavg{2,i} = reshape(nanmean(vertcat(tmp{idx_hori,i}),1),tmpdims);
                end
                % reshape and permute back to fit the new dimensions
                tmpavg = reshape(tmpavg,newdims(permdim));
                tmpavg = ipermute(tmpavg,permdim);
                % store in new field
                tmpnew.(marker{nm}).(axisname{an}).(cellflds{fn}) = tmpavg;
            end
        end
    end
    
    % add verthori
    newf = length(epe.fact) + 1;
    epe.nfact = newf;
    epe.fact{newf} = newfact;
    epe.fname{newf} = newfname;
    if ~ismember('verthori',epe.bfact)
        epe.bfact = [epe.bfact {'verthori'}];
        epe.nbfact = epe.nbfact + 1;
    end
    if ~isfield(epe.lvl,'verthori')
        epe.lvl.verthori = 1:2;
    end
    
    % store new factor
    epe.(newfname) = tmpnew;
    tmplvl = factlvl;
    tmplvl(ismember(fact,'slant')) = {epe.lvl.verthori};
    epe.(newfname).lvl = arraycomb(tmplvl{:});
        
end