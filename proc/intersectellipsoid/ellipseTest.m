%define an origin centred ellipsoid
ellipsMatrix    = [6, 4, 2];
ellipseRotation = eye(3);
%define an origin centred plane
normalVector    = [15, -4, 40]';
%find the intersecting ellipse
[newEllipse, newRotation] = planeEllipsoidIntersection(ellipsMatrix, normalVector, ellipseRotation);

%sample the ellipsoid for plotting
[x, y, z] = ellipsoid(0, 0, 0, 1, 1, 1);
xyz = [x(:) y(:) z(:)] * diag(ellipsMatrix) * ellipseRotation;
x(:) = xyz(:, 1); 
y(:) = xyz(:, 2); 
z(:) = xyz(:, 3);
%sample the ellipse for plotting
t = linspace(0, 2 * pi, 40)';
circle = [cos(t) sin(t) zeros(length(t),1)];
xyzEllipse = circle * diag(newEllipse) * newRotation';
%sample the plane for plotting
[xx, yy] = ndgrid(-10:1:10, -10:1:10);
zz = (-normalVector(1) * xx - normalVector(2) * yy) / normalVector(3);

figure;
%plot the ellipsoid surface
surf(x, y, z); colormap('gray'); axis equal; hold on;
%plot the plane
surf(xx, yy, zz); alpha(0.3);
%plot the ellipse
plot3(xyzEllipse(:,1),xyzEllipse(:,2),xyzEllipse(:,3),'color','k', 'linewidth', 3);

% I think it works, though there might be some small issues with rotations
% still
% 
% Happy new Year :)


