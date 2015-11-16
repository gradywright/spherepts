function x = getCubedSphNodes(p)
%GETCUBEDSPHNODES Computes the equiangluar cubed sphere node set for the 
%   surface of the sphere.
%
%   X = getCubedSphNodes(P) returns an N-by-3 matrix containing the 
%   cubed-sphere nodes generated using an equiangular projection. 
%   The columns corresponds to the (x,y,z) cordinates of the nodes.  
%   Each face of the projection of the inscribed cube has P^2 points
%   corresponding to an equiangular projection.  This gives the total
%   number of unique nodes as:
%                            N = 6*P^2 - 12*P + 8;
%
%   See, for example:
%   Nair, R. D., S. J. Thomas and R. D. Loft, 2005: A discontinuous
%   Galerkin transport scheme on the cubed sphere.  Monthly Weather
%   Review,Vol. 133, pp 814-828
%   For more details on the nodes.
%
%   Example:
%       x = getCubedSphNodes(21);  % Gives N = 2402 nodes

% Author: Grady Wright, 2014

[alp,bet] = meshgrid(linspace(-pi/4,pi/4,p));

xl = sqrt(3)/3*tan(alp); xl = xl(:);
yl = sqrt(3)/3*tan(bet); yl = yl(:);
shft = sqrt(3)/3*ones(size(xl)); shft = shft(:);
r = ones(size(xl));

% Face 5: Cube covering North Pole
xx = xl; yy = yl; zz = shft;
[lam,th] = cart2sph(xx,yy,zz);
[X5,Y5,Z5] = sph2cart(lam,th,r);

% Face 6: Cube covering South Pole
X6 = X5; Y6 = Y5; Z6 = -Z5;

% Face 1: Covering 1/4 of equator:
xx = xl; yy = shft; zz = yl;
[lam,th] = cart2sph(xx,yy,zz);
[X1,Y1,Z1] = sph2cart(lam,th,r);

% Face 3: Covering 1/4 of equator opposite Face 1
X3 = X1; Y3 = -Y1; Z3 = Z1;

% Face 2: Covering 1/4 of equator
xx = shft; yy = xl; zz = yl;
[lam,th] = cart2sph(xx,yy,zz);
[X2,Y2,Z2] = sph2cart(lam,th,r);

% Face 4: Covering 1/4 of equator opposite Face 2
X4 = -X2; Y4 = Y2; Z4 = Z2;

x = [[X1 Y1 Z1];[X2 Y2 Z2];[X3 Y3 Z3];[X4 Y4 Z4];[X5 Y5 Z5];[X6 Y6 Z6]];

% Remove duplicate points on the edges of the projected cube.
% Rounding errors in computing the nodes may not eliminate duplicate
% points.  So force this to happen to a reasonable tolerance.
[~,ia] = unique(round(1e8*x)/1e8,'rows');  

x = x(ia,:);

end

