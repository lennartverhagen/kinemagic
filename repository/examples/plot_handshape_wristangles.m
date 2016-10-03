tseries = tmp;
tmp = tseries;
tmp2 = cfg.tseries.gripapt.marker;

tseries.gripapt.marker = {'MCP','i','t','rd'};
cfg.tseries.gripapt.marker = {'MCP','i','t','rd'};
for i = 1:7
    dat = tseries.gripapt.slant.mean{i};
    
    % get sides
    a = squeeze(dat(1,1,:));
    b = squeeze(dat(2,1,:));
    c = squeeze(dat(3,1,:));
    
    % use cosine-rule to obtain angles
    alpha = acos( (b.^2 + c.^2 - a.^2) ./ (2.*b.*c) );
    beta  = acos( (a.^2 + c.^2 - b.^2) ./ (2.*a.*c) );
    gamma = acos( (a.^2 + b.^2 - c.^2) ./ (2.*a.*b) );
    
    % convert to degrees
    dat = reshape([alpha'; beta'; gamma'],[3 1 101]);
    dat = real(dat).*180./pi;
    
    % calculate disection ratio of MCP projected on t-i line
    dat(4,1,:) = (cos(gamma).*b) ./ (cos(beta).*c);
    %dat(5,1,:) = (pi/2 - gamma) ./ ( pi/2 - beta);
    
    % store
    tseries.gripapt.slant.mean{i} = dat;
end


cfg.tseries.gripapt.marker = 't';
km_tseriesplot(cfg,tseries);
cfg.tseries.gripapt.marker = 'i';
km_tseriesplot(cfg,tseries);
cfg.tseries.gripapt.marker = 'MCP';
km_tseriesplot(cfg,tseries);
cfg.tseries.gripapt.marker = 'rd';
km_tseriesplot(cfg,tseries);


% draw hand at 60% of movement, relative to t-i line with t = 0;
figure;
n = 61;
lncl = str2rgb('rgbcmoh');
for i = 1:7
    
    % get angle
    g = tseries.gripapt.slant.mean{i}(3,1,n);
    g = g ./ (180./pi); % in radians
    
    % get distance
    dtMCP = tmp.gripapt.slant.mean{i}(2,1,n);
    dti = tmp.gripapt.slant.mean{i}(1,1,n);    
    
    x = [0 (cos(pi/2 - g).*dtMCP) 0];
    y = [0 (sin(pi/2 - g).*dtMCP) dti];
    
    plot(x,y,'Color',lncl(i,:),'Marker','o','LineStyle','-','LineWidth',2);
    hold on;
end
axis equal


% draw hand at 60% of movement, relative to t-i line with MCP = 0;
figure;
n = 61;
lncl = str2rgb('rgbcmoh');
for i = 1:7
    
    % get angle
    g = tseries.gripapt.slant.mean{i}(3,1,n);
    g = g ./ (180./pi); % in radians
    
    % get distance
    dtMCP = tmp.gripapt.slant.mean{i}(2,1,n);
    dti = tmp.gripapt.slant.mean{i}(1,1,n);    
    
    x = [-(cos(pi/2 - g).*dtMCP) 0 0-(cos(pi/2 - g).*dtMCP)];
    y = [-(sin(pi/2 - g).*dtMCP) 0 dti-(sin(pi/2 - g).*dtMCP)];
    
    plot(x,y,'Color',lncl(i,:),'Marker','o','LineStyle','-','LineWidth',2);
    hold on;
end
axis equal


% draw hand at 60% of movement, relative to t-i line with meanti = 0;
figure;
n = 61;
lncl = str2rgb('rgbcmoh');
for i = 1:7
    
    % get angle
    g = tseries.gripapt.slant.mean{i}(3,1,n);
    g = g ./ (180./pi); % in radians
    
    % get distance
    dtMCP = tmp.gripapt.slant.mean{i}(2,1,n);
    dti = tmp.gripapt.slant.mean{i}(1,1,n);    
    
    x = [0 (cos(pi/2 - g).*dtMCP) 0];
    y = [-dti/2 (sin(pi/2 - g).*dtMCP)-dti/2 dti/2];
    
    plot(x,y,'Color',lncl(i,:),'Marker','o','LineStyle','-','LineWidth',2);
    hold on;
end
axis equal