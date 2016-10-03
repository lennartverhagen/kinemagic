%% fake some data
if true
    n = 200;
    pos = randn(n,3) .* (ones(n,1)*[1 2 3]);
    orig_pos = pos;
else
    load example_pos
end


%% rotate it 30 degrees in the xy plane
%----------------------------------------
r    = pi/4;
rotx = [1     	0     	0
        0       cos(r)	-sin(r)
        0       sin(r)	cos(r)  ];
roty = [cos(r)	0       sin(r)
        0       1      	0
        -sin(r) 0       cos(r) 	 ];
rotz = [cos(r)	-sin(r)	0
        sin(r)	cos(r)	0
        0     	0     	1       ];
pos = orig_pos * rotx * roty * rotz;


%% use my KineMagic script to get the confidence ellipsoid.
%----------------------------------------
CI = km_mdCI(pos,0.95,'equal','chisq');
%CI = km_mdCI(pos,0.95,'equal','f');


%% plot km_mdCI
%----------------------------------------
% plot the datapoints
figure(1); hold off;
idx_out = CI.xyz.idx_out==1;
idx_in = ~idx_out;
plot3(pos(idx_in,1),pos(idx_in,2),pos(idx_in,3),'marker','o','linestyle','none','color','b'); hold on;
plot3(pos(idx_out,1),pos(idx_out,2),pos(idx_out,3),'marker','o','linestyle','none','color','r','markerface','r'); hold on;
title(sprintf('confidence interval: %d%%\nnumber of outliers: %d/%d = %0.1f%%',100*CI.conf,sum(idx_out),size(pos,1),100*sum(idx_out)/size(pos,1)));

% plot the 2D ellipses
xyz = CI.xy.xyz;
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'color','r','linewidth',3);
xyz = CI.zx.xyz;
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'color','b','linewidth',3);
xyz = CI.yz.xyz;
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'color','g','linewidth',3);

% plot the 3D ellipsoid
xyz = CI.xyz.xyz;
% mesh(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),'edgecolor','k','facecolor','k');
surf(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3)); colormap gray
alpha(0.3); camlight; axis tight; axis equal;

% save figure
if ~ispc
    figdir = fullfile(filesep,'home','action','lenver','tmp');
    if ~exist(figdir,'dir'),  mkdir(figdir);  end
    figname = fullfile(figdir,'mdCI_example.pdf');
    export_fig(figname,'-pdf');
end