function [xh,yh] = sph2hammer(lam,th)
%SPH2HAMMER Maps points in spherical coordinates to points in Hammer
%coordinates.  
%
% [XH,YH] = sph2hammer(LAM,TH) maps points on the surface of the sphere
%   represented in spherical coordiantes to points in the 2D Hammer
%   coordinates.  This is a popular projection for plotting on the sphere.
%   LAM is the azimuthal angle (in radians) of the point and TH is the 
%   elevation angle (in radians) measured from the equator 
%   (i.e. -pi/2 <= TH <= pi/2).
%
%   Example: Visualize a function on the sphere using Hammer projection
%     [LAM,TH] = meshgrid(linspace(-pi,pi,101),linspace(-pi/2,pi/2,51));
%     f = sphHarm(4,0,LAM,TH) + sqrt(14/11)*sphHarm(4,4,LAM,TH);
%     [XH,YH] =  sph2hammer(LAM,TH);
%     pcolor(XH,YH,f);
%     daspect([1 1 1]);

xh = 2*sqrt(2)*cos(th).*sin(lam/2)./sqrt(1 + cos(th).*cos(lam/2));
yh = sqrt(2)*sin(th)./sqrt(1 + cos(th).*cos(lam/2));

end
