function [lam,th,us,vs] = cartv2sphv(x,y,z,uc,vc,wc)
%CARTV2SPHV Converts vectors tangent to the sphere from Cartesian to
%spherical coordinates.
%
%   [LAM,TH,US,VS] = cartv2sphv(X,Y,Z,UC,VC,WC) converts a vector field
%   tangent to the sphere that is expressed with respect to Cartesian
%   coordinates to a one that is expressed with respect to spherical
%   coordinates. X, Y, Z are the locations where the field is sampled, with
%   UC, VC, and WC being the values of the field in the (X,Y,Z) directions.
%   LAM and TH are the respective azimuthal and elevation (longitude and
%   latitude) values of X, Y, and Z, where -pi < LAM <= PI and 
%   -pi/2 <= TH <= pi/2.  US and VS are the values of the field in the LAM
%   and TH directions, respectively.
%
%   Example:
%       [x,y,z] = sphere(20);
%       u = -y;
%       v = x;
%       w = 0*x;
%       [lam,th,us,vs] = cartv2sphv(x,y,z,u,v,w);
%       quiver(lam,th,us,vs)

sz = size(x);

% Flatten all the inputs:
x = x(:); y = y(:); z = z(:); uc = uc(:); vc = vc(:); wc = wc(:);

[lam,th] = cart2sph(x,y,z);

% Vectors for converting the field in Cartesian coordinates to a field
% in spherical coordinates.
c2s_u = [-sin(lam) -sin(th).*cos(lam)];
c2s_v = [cos(lam) -sin(th).*sin(lam)];
c2s_w = [zeros(size(lam)) cos(th)];

% Convert the vectors.
us = c2s_u(:,1).*uc + c2s_v(:,1).*vc + c2s_w(:,1).*wc;
vs = c2s_u(:,2).*uc + c2s_v(:,2).*vc + c2s_w(:,2).*wc;

us = reshape(us,sz);
vs = reshape(vs,sz);

end

