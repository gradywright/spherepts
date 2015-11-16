function [th,phi,r] = cart2sphm(x)
%CART2SPHM Transform Cartesian to spherical coordinates.
%   [TH,PHI,R] = cart2sph(X) transforms corresponding elements of
%   data stored in Cartesian coordinates to spherical
%   coordinates (azimuth TH, elevation PHI, and radius R). The matrix X is
%   assumed to contain N rows and 3 columns where each row contains the
%   (x,y,z) Cartesian coordinate of the point.
%
%   Same as matlab function cart2sph, but this accepts a single matrix with
%   the data.

% Author: Grady Wright, 2014

[th,phi,r] = cart2sph(x(:,1),x(:,2),x(:,3));

end
