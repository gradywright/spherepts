function Y = sphHarm(l,m,c1,c2,c3)
%SPHHARM Computes the normalized, real-valued, spherical harmonic of degree
%   L order M at a given set of locations on the sphere.
%
%   Y = sphHarm(L,M,LAM,TH) returns the degree L order M normalized
%   spherical harmonic at the points (LAM,TH) on the sphere expressed in
%   longitude-latitude coordinates (or azimuthal-elevation).  Here
%            -pi <= lam <= pi   is the longitude (azimuthal) coordinate and
%           -pi/2 <= th <= pi/2 is the latitude (elevation) coordinate.
%   
%   Y = SPHHARM(L,M,X,Y,Z) returns the degree L order M normalized
%   spherical harmonic at the points (x,y,z) on the sphere expressed in
%   Cartesian coordinates.
%
%   Example 1: Spherical coordinates
%       [lam,th] = meshgrid(linspace(-pi,pi,81),linspace(-pi/2,pi/2,41));
%       f = sphHarm(6,0,lam,th) + sqrt(14/11)*sphHarm(6,5,lam,th);
%       surf(lam,th,f), shading interp;
%
%   Example 2: Cartesain coordinates
%       [x,y,z] = sphere(101);
%       f = sphHarm(6,0,x,y,z) + sqrt(14/11)*sphHarm(6,5,x,y,z);
%       surf(x,y,z,f), shading interp, axis equal;

% Author: Grady Wright, 2014

% The degree l must be greater than or equal to the magnitude of the order
% m
if l < abs(m)
    error('RBFSPHERE:sphHarm','The degree of the spherical harmonic must be greater than or equal to the magnitude of the order');
end

if nargin == 4      % Spherical coordinates used.
    lam = c1;
    th = c2;
    c1 = []; c2 = [];
elseif nargin == 5  % Cartesian coordinates used.
    % Convert to spherical
    [lam,th] = cart2sph(c1,c2,c3);
    c1 = []; c2 = []; c3 = [];
else
    error('RBFSPHERE:sphHarm','Proper number of input arguments not given');
end

% Flatten and transpose th and lam so they work with the legendre function
sz = size(th); th = th(:)'; lam = lam(:)';

% Normalization
a = sqrt((2*l+1)/2/pi*factorial(l-abs(m))/factorial(l+abs(m))*(2-double(m==0)));
Y = legendre(l,sin(th));
% Get the right associated legendre function
Y = squeeze(Y(abs(m)+1,:,:));
% Determine if the cos or sin term should be added.
pos = abs(max(0,sign(m+1)));
% Compute the spherical harmonic
Y = (pos*cos(m*lam) + (1-pos)*sin(m*lam)).*(a*Y);
% Reshape so it is the same size as the th and lam that were passed in.
Y = reshape(Y,sz);

end