function [newEllipse, newRotation] = planeEllipsoidIntersection(ellipsoid, planeNormal, rotation) 
%FINDELLIPSE
%   [EO, RO] = FINDELLIPSE(EI, P, RI)
%   Given an ellipsoid and a plane, the intersecting ellipse will be found
%   EI is the input ellipsoid (3 X 1 matrix) giving the three principle
%   axes of the ellipsoid. The ellipsoid must be origin centred.
%	P is the plane normal for a plane going through the origin.
%   RI is the input rotation of the ellipsoid, given by a 3 X 3 matrix.
%
%   EO is the output ellipse given by the three principle axes. The third
%   coordinate will always be 0. This is done for consistency with the rest
%   of the three-dimensional objects
%   RO is the rotation of the output ellipse
%
%   The method follows the line of reasoning given explained in:
%   http://mathforum.org/library/drmath/view/71275.html
%   The rest of the mathematics has been worked out by Tim van Mourik
%
%   Written by Tim van Mourik, 30 December 2012

if any(ellipsoid <= 0)
    error('TVM:findEllipse:InvalidEllipse', 'All axes of the ellipse should be positive')
end
if nargin < 3
    rotation = eye(3);
end

%normalise the plane normal vector
planeNormal = planeNormal / norm(planeNormal);
%rotate the normal with the inverse rotation of the ellipse
newNormal = rotation \ planeNormal;

%the coordinates are sorted in order to take the best samples possible
[~, index] = min(newNormal .^2);
if index == 1
    indices = [1, 2, 3];
elseif index == 2
    indices = [2, 1, 3];
elseif index == 3
    indices = [3, 2, 1];
end

%the ellipse parameters
a = ellipsoid(indices(1));
b = ellipsoid(indices(2));
c = ellipsoid(indices(3));
%the plane normal paramters
d = newNormal(indices(1));
e = newNormal(indices(2));
f = newNormal(indices(3));

%compute three xyz values. These three values always exist along the
%correct axis.
x = [-min(ellipsoid), 0, min(ellipsoid)] / sqrt(3);
y = (- d * e * x + sqrt((d * e * x) .^ 2 + ((c * f / b) ^ 2 + e * e) * (c * c * f * f - ((c * f / a) ^ 2 + d * d) * x .^ 2))) / ((c * f / b) ^ 2 + e * e);
z = - (d * x + e * y) / f;
%a matrix of samples
samples = [x; y; z];
%reshuffle to the actual xyz-coordinates
samples = [samples(indices(1), :); samples(indices(2), :); samples(indices(3), :)];

%find the rotation of the plane
rotationMatrix = findRotation(newNormal);
%rotate the ellipse with this rotation. After rotation, the z-coordinates
%should be zero
ellipsePoints = rotationMatrix * samples;

%use the three samples to find the three unknowns: rotation and size in x
%and y
coefficients = [ellipsePoints(1, :) .^ 2; ellipsePoints(2, :) .^ 2; ellipsePoints(1, :) .* ellipsePoints(2, :)]';
CDE = coefficients \ ones(3, 1);

%compute the new angle
theta   = atan2(CDE(3), (CDE(2) - CDE(1))) / 2;
A       = cos(theta);
B       = sin(theta);

%the orientation of the 2D ellipse is given by:
orientation = [A, -B, 0;
               B,  A, 0;
               0,  0, 1];
           
%The size of the ellipse is:
sizeX       = sqrt((A ^ 4 - B ^ 4) / (A ^ 2 * CDE(1) - B ^ 2 * CDE(2)));
sizeY       = sqrt((A ^ 4 - B ^ 4) / (A ^ 2 * CDE(2) - B ^ 2 * CDE(1)));
newEllipse  = [sizeX, sizeY, 0];

%The rotation of the new ellipse is the identity matrix...
newRotation = eye(3);
%...transformed back with the right orientation...
newRotation = orientation \ newRotation;
%...transformed back from plane space to ellipse space...
newRotation = rotationMatrix \ newRotation;
%...transformed back from ellipse space to world space.
newRotation = rotation \ newRotation;

% Happy New Year :)

end %end function

function rotationMatrix = findRotation(normal)
%sets the rotation angles
alpha = acos(normal(3));
beta  = 0;
gamma = atan2(normal(1), normal(2));

%computes the sine and cosine of the angle
sinAlpha = sin(alpha);
sinBeta  = sin(beta);
sinGamma = sin(gamma);
cosAlpha = cos(alpha);
cosBeta  = cos(beta);
cosGamma = cos(gamma);

%sets the rotation matrices
rx = [1, 0, 0;
      0, cosAlpha, -sinAlpha;
      0, sinAlpha, cosAlpha];
  
ry = [cosBeta,  0, sinBeta;
      0,        1, 0;
      -sinBeta, 0, cosBeta];
  
rz = [cosGamma, -sinGamma,  0;
      sinGamma, cosGamma,   0;
      0,        0,          1];

rotationMatrix = rx * ry * rz;

end %end function





