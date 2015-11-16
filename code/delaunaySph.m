function tri = delaunaySph(x)
% DELAUNAYSPH Computes the Delaunay triangulation for a set of nodes on the
%   surface of the sphere in 3D space.
%
%   TRI = delaunaySph(X) creates a 3D Delaunay triangulation of the points
%   in X, where X is assumed to be an array of size N-by-3, with each row
%   corresponding to the (x,y,z) coordinates of a point on the sphere.
%   This gives a surface triangulation so that no triangle edges go through
%   the sphere.
%
%   Example:
%       p = (1+sqrt(5))/2;
%       x = [[0,p,1];[0,-p,1];[0,p,-1];[0,-p,-1];[1,0,p];[-1,0,p];...
%            [1,0,-p];[-1,0,-p];[p,1,0];[-p,1,0];[p,-1,0];[-p,-1,0]];
%       tri = delaunaySph(x);
%       trisurf(tri,x(:,1),x(:,2),x(:,3),1+0*x(:,1));
%       axis equal;

% Author: Grady Wright, 2014

% Make sure the nodes are on the unit sphere.
x = bsxfun(@rdivide,x,sqrt(sum(x.^2,2)));

% 3D Delaunay "triangulation"
tri = delaunay(x);
% Remove edges through the sphere
tri = freeBoundary(TriRep(tri,x));

end
