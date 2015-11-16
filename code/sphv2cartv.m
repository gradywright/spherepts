function [x,y,z,uc,vc,wc] = sphv2cartv(lam,th,us,vs)
%SPHV2CARTV Converts vectors tangent to the sphere from spherical coordinates to
%Cartesian coordinates.
%
%   [X,Y,Z,UC,VC,WC] = cartv2sphv(X,Y,Z,UC,VC,WC) converts a vector field
%   tangent to the sphere that is expressed with respect to spherical
%   coordinates to a one that is expressed with respect to Cartesian
%   coordinates. LAM and TH are the respective azimuthal and elevation
%   (longitude and latitude) locations where the field is sampled, where
%   -pi < LAM <= PI and -pi/2 <= TH <= pi/2.  US and VS being the values of
%   the field in those directions.  X, Y, Z are the respective (X,Y,Z)
%   Cartesian coordinates of LAM and TH. UC, VC, WC are the values of the
%   field in the X, Y, and Z directions, respectively.
%
%   Example:
%       [lam,th] = meshgrid(-pi:pi/16:pi,-pi/2:pi/8:pi/2);
%       u = sin(th).*cos(lam);
%       v = -sin(lam);
%       [x,y,z,uc,vc,wc] = sphv2cartv(lam,th,u,v);
%       sphere(100); 
%       colormap(jet(1)); shading flat; hold on;
%       quiver3(x,y,z,uc,vc,wc,'k-');
%       axis([-1 1 -1 1 -1 1]);

sz = size(lam);

[x,y,z] = sph2cart(lam,th,1+0*lam);

% Flatten all the inputs:
th = th(:); lam = lam(:); us = us(:); vs = vs(:);

% Transformation for converting the latitudinal velocity to Cartesian velocity.
s2c_u = [-sin(lam) cos(lam) zeros(size(lam))];
% Transformation for converting the longitudinal velocity to Cartesian velocity.
s2c_v = [-cos(lam).*sin(th) -sin(lam).*sin(th) cos(th)];

uc = s2c_u(:,1).*us + s2c_v(:,1).*vs;
vc = s2c_u(:,2).*us + s2c_v(:,2).*vs;
wc = s2c_u(:,3).*us + s2c_v(:,3).*vs;

uc = reshape(uc,sz);
vc = reshape(vc,sz);
wc = reshape(wc,sz);

end


