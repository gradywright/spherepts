function r = voronoiDiagramSph(x,tri)
%voronoiDiagramSph   Returns the Voronoi diagram of a set of nodes on the
%   surface of the sphere in 3D space.
%
%   R = voronoiDiagramSph(X) Returns a cell array containing the Voronoi
%   region for each node in X.  X is assumed to be a N-by-3 array with the
%   (x,y,z) components of the nodes on the surface of the sphere.  R{j}
%   contains the vertices for the Voronoi region for node X(j,:).
%
%   R = voronoiDiagramSph(X,TRI) uses the triangulation of the nodes TRI
%   instead of computing one internally.
%
%   Example:
%       [x,tri] = getIcosNodes(0,1);
%       r = voronoiDiagramSph(x,tri);
%       for j=1:size(r,1), vv=r{j}; patch(vv(:,1),vv(:,2),vv(:,3),size(vv,1)); end
%       view(3), axis square, colormap(copper);
%
% See also DELAUNAYSPH, VORONOISPH

% Author: Grady Wright, 2014

if ~exist('tri','var')
    tri = delaunaySph(x);
end

% Get the radius of the sphere the nodes lie on.
rad = sqrt(x(1,1).^2 + x(1,2).^2 + x(1,3).^2);

% Do all computations on the unit sphere.
x = bsxfun(@rdivide,x,sqrt(sum(x.^2,2)));

r = cell(size(x,1),1);

[m,n] = size(tri);

% Loop over each vertex of the triangulation and compute its circumcenter.
for j = 1:size(x,1)
    [v1,temp] = find(j == tri(:,1));
    [v2,temp] = find(j == tri(:,2));
    [v3,temp] = find(j == tri(:,3));
    tri_cell = [tri(v1,:);tri(v2,:);tri(v3,:)];
    % Compute the circumcenter of each triangle
    v = computeCircumcenter(tri_cell,x);
    
    %
    % Determine the proper ordering of the circumcenters.  This is done by
    % rotating the nodes to the north pole (where the x(j,:)) is at the
    % north pole and then determining the ordering based on the projection
    % onto the x-y plane.
    %

    % Convert the point to latitude longitude
    [lam,th,~] = cart2sphm(x(j,:));
    
    % Rotate the points to the north pole.
    plam = lam;
    pth = pi-th;
    
    % Rotation matrices
    D = [[cos(plam) sin(plam) 0];[-sin(plam) cos(plam) 0];[0 0 1]];
    C = [[sin(pth) 0 cos(pth)];[0 1 0];[-cos(pth) 0 sin(pth)]];
    vt = (C*(D*v.')).';
    
    % Convert the rotated points to latitude longitude
    lam = cart2sphm(vt);
    
    % Sort the points by their azimuthal angle.
    [temp,id] = sort(lam);
    
    v = v(id,:);
    
    % Scale to the sphere the nodes are on and store the results.
    r{j} = rad*v;
end

% d = zeros(size(tri,1),3);
% for j=1:size(tri,1)
%     cc = computeCircumcenter(tri(j,:),x);
%     vt = x(tri(j,:),:);
%     d(j,:) = acos(cc*vt.');
% end

end

function xc = computeCircumcenter(tri,nodes)

% See http://en.wikipedia.org/wiki/Circumcenter for the details of the
% calculation.

v1 = nodes(tri(:,1),:);
v2 = nodes(tri(:,2),:);
v3 = nodes(tri(:,3),:);
e1 = v2-v1;
e2 = v3-v1;

xc = 0*v1;
xc(:,1) = e1(:,2).*e2(:,3) - e1(:,3).*e2(:,2);
xc(:,2) = e1(:,3).*e2(:,1) - e1(:,1).*e2(:,3);
xc(:,3) = e1(:,1).*e2(:,2) - e1(:,2).*e2(:,1);

%
% Project onto the unit sphere.
%

% Get the radius of the sphere the nodes lie on.
rr = sqrt(xc(:,1).^2 + xc(:,2).^2 + xc(:,3).^2);
xc = bsxfun(@rdivide,xc,rr);

end
