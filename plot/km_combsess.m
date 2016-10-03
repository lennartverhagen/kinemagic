function [data] = km_combsess(sessdata,flg)
%--------------------------------------------------------------------------
% append multiple sessions into one
%
% See also KM_PLOT
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-02-01
%--------------------------------------------------------------------------

% return quickly if only one session is requested
nrsess = length(sessdata);
if nrsess == 1
    data = sessdata{1};
    return
end

switch lower(flg)
    case {'raw','tseries'}
        % append the indeces arrays
        data = structcat(sessdata,1,{},{'dimord','axis','marker'});
        
        % continue if no data is provided
        if isempty(data) || isempty(fieldnames(data))
            return
        end
        
        % set stable fields
        data.dimord = sessdata{1}.dimord;
        data.axis   = sessdata{1}.axis;
        data.marker = sessdata{1}.marker;
        
    case 'param'
        % vertically concatinate the sessions
        data = vertcat(sessdata{:});
        
        
    case 'ci'
        % initialize the data
        data = sessdata{1};
        
        % continue if no epe is provided
        nsess = length(sessdata);
        if isempty(data) || nsess  == 1
            return
        end
        
        % loop over factors
        nfact = length(data.fact);
        for f = 1:nfact
            
            % get factor(s)
            fact = data.fact{f};
            fname = data.fname{f};
            
            % get levels
            factlvl = cell(1,length(fact));
            for ff = 1:length(fact)
                factlvl{ff} = data.lvl.(fact{ff});
            end
            
            % determine dimensions
            dims = cellfun(@(x) length(x),factlvl);
            if length(dims) == 1
                dims = [dims 1];
            end
            
            % add sess
            newf = length(data.fact) + 1;
            data.nfact = newf;
            data.fact{newf} = [{'sess'} fact];
            newfname = ['sessx' fname];
            data.fname{newf} = newfname;
            if ~ismember('sess',data.bfact)
                data.bfact = [data.bfact {'sess'}];
                data.nbfact = data.nbfact + 1;
            end
            if ~isfield(data.lvl,'sess')
                data.lvl.sess = 1:nsess;
            end
            
            % loop over parameters
            tmpnew = struct;
            tmpavg = struct;
            param = fieldnames(data.(fname));
            param = param(~ismember(param,'lvl'));
            for np = 1:length(param)
                tmppar = data.(fname).(param{np});
                
                % loop over fields
                flds = fieldnames(tmppar);
                %flds = flds(~ismember(flds,{'eigvec','eigval','xyz'}));
                for fn = 1:length(flds)
                    % loop over sessions
                    tmp = cell(nsess,1);
                    for ss = 1:nsess
                        tmp{ss} = reshape([sessdata{ss}.(fname).(param{np}).(flds{fn})],[1 dims]);
                    end
                    tmp = vertcat(tmp{:});
                    tmpnew.(param{np}).(flds{fn}) = tmp;
                    tmpavg.(param{np}).(flds{fn}) = reshape(nanmean(tmp,1),dims);
                end
                
            end
            
            % store session data
            tmplvl = [{data.lvl.sess}; factlvl(:)];
            data.(newfname) = tmpnew;
            data.(newfname).lvl = arraycomb(tmplvl{:});
            
            % store averaged data
            %tmpavg.lvl = data.(fname).lvl;
            %tmpavg.lvl = data.lvl.(fname);
            tmpavg.lvl = arraycomb(factlvl{:});
            data.(fname) = tmpavg;
            
        end
        
    case 'epe'
        % initialize the data
        data = sessdata{1};
        
        % continue if no epe is provided
        nsess = length(sessdata);
        if isempty(data) ||nsess  == 1
            return
        end
        
        % loop over factors
        nfact = length(data.fact);
        for f = 1:nfact
            
            % get factor(s)
            fact = data.fact{f};
            fname = data.fname{f};
            
            % get levels
            factlvl = cell(1,length(fact));
            for ff = 1:length(fact)
                factlvl{ff} = data.lvl.(fact{ff});
            end
            
            % determine dimensions
            dims = cellfun(@(x) length(x),factlvl);
            if length(dims) == 1
                dims = [dims 1];
            end
            
            % add sess
            newf = length(data.fact) + 1;
            data.nfact = newf;
            data.fact{newf} = [{'sess'} fact];
            newfname = ['sessx' fname];
            data.fname{newf} = newfname;
            if ~ismember('sess',data.bfact)
                data.bfact = [data.bfact {'sess'}];
                data.nbfact = data.nbfact + 1;
            end
            if ~isfield(data.lvl,'sess')
                data.lvl.sess = 1:nsess;
            end
            
            % loop over markers
            tmpnew = struct;
            tmpavg = struct;
            marker = fieldnames(data.(fname));
            marker = marker(~ismember(marker,'lvl'));
            for nm = 1:length(marker)
                tmpmark = data.(fname).(marker{nm});
                
                % loop over axes
                axisname = fieldnames(tmpmark);
                for an = 1:length(axisname)
                    
                    % loop over fields
                    flds = fieldnames(tmpmark.(axisname{an}));
                    paramflds = flds(~ismember(flds,{'eigvec','eigval','xyz'}));
                    for fn = 1:length(paramflds)
                        % loop over sessions
                        tmp = cell(nsess,1);
                        for ss = 1:nsess
                            tmp{ss} = reshape(sessdata{ss}.(fname).(marker{nm}).(axisname{an}).(paramflds{fn}),[1 dims]);
                        end
                        tmp = vertcat(tmp{:});
                        tmpnew.(marker{nm}).(axisname{an}).(paramflds{fn}) = tmp;
                        tmpavg.(marker{nm}).(axisname{an}).(paramflds{fn}) = reshape(nanmean(tmp,1),dims);
                    end
                    
                    % FIXME: in principle the other fields could be averaged as
                    % well... But personally, I think it would be best to calculate
                    % new coordinates based on the eigenvalues and vectors at the
                    % time of reporting
                    cellflds = flds(ismember(flds,{'eigvec','eigval','xyz'}));
                    for fn = 1:length(cellflds)
                        % loop over sessions
                        tmp = cell(nsess,1);
                        for ss = 1:nsess
                            tmp{ss} = reshape(sessdata{ss}.(fname).(marker{nm}).(axisname{an}).(cellflds{fn}),[1 dims]);
                        end
                        tmp = vertcat(tmp{:});
                        tmpnew.(marker{nm}).(axisname{an}).(cellflds{fn}) = tmp;
                        % reshape tmp to nsess x prod(dims)
                        tmp = reshape(tmp,[nsess prod(dims)]);
                        % reshape fields to nsess x prod(tmpdims)
                        tmpdims = size(tmp{1});
                        tmp = cellfun(@(x,y) reshape(x,[1 prod(y)]),tmp,repmat({tmpdims},size(tmp)),'UniformOutput',false);
                        % calculate the average over the first dimension
                        tmp = mat2cell(nanmean(cell2mat(tmp),1),1,repmat(prod(tmpdims),[1 prod(dims)]));
                        % reshape back to tmpdims
                        tmp = cellfun(@(x,y) reshape(x,y),tmp,repmat({tmpdims},[1 prod(dims)]),'UniformOutput',false);
                        %reshape back to dims
                        tmp = reshape(tmp,dims);
                        tmpavg.(marker{nm}).(axisname{an}).(cellflds{fn}) = tmp;
                    end
                    
                end
            end
            
            % store session data
            tmplvl = [{data.lvl.sess}; factlvl(:)];
            data.(newfname) = tmpnew;
            data.(newfname).lvl = arraycomb(tmplvl{:});
            
            % store averaged data
            tmpavg.lvl = data.(fname).lvl;
            data.(fname) = tmpavg;
            
        end
        
    case 'idxc'
        % check if idxc is a structure or a matrix
        if ~isstruct(sessdata{1})
            data = vertcat(sessdata{:});
            return
        end
        
        % append the indeces arrays
        data = structcat(sessdata,1,{'cell','char'});
        
        % assign trials to session numbers
        t = 1;
        for ss = 1:length(sessdata)
            ntrl = size(sessdata{ss}.all,2);
            sessname = sprintf('sess%d',ss);
            data.(sessname) = false(size(data.all));
            data.(sessname)(:,t:(t+ntrl-1)) = true;
            t = t+ntrl;
        end
    
    case {'idx','idxr'}
        % check if idx is a structure or a matrix
        if ~isstruct(sessdata{1})
            data = vertcat(sessdata{:});
            return
        end
        
        % append the indeces arrays
        data = structcat(sessdata,1,{'cell','char'});
    
    case 'cfg'
        % initialize the configuration structure
        data = sessdata{1};
        
        % append the trl matrix in the configuration structure
        tmp = km_findcfg(sessdata,'trl');
        if ~isempty(tmp)
            tmp = cellfun(@(x) km_findcfg(x,'trl'),sessdata,'UniformOutput',false);
            data.trl = cat(1,tmp{:});
        end
        
        % the selectart array is no longer valid
        if isfield(data,'selectart')
            data = rmfield(data,'selectart');
        end

    otherwise
        error('Unrecognized class/type of data. Nothing to append')
end



