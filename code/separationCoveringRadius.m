function [delta,eta] = separationCoveringRadius(X)
%separationCoveringRadius Computes the separation and covering radius of a node
%   set X on the unit sphere.
%
%   [delta,eta] = separationCovering(X) computes the separation
%   (delta) and covering radius (eta) for the given node set on the unit
%   sphere.  For a set X with N nodes, these quantities are defined as
%   follows:
%  
%        delta = min dist(x_i,x_j), i,j=1,...,N, i~=j,
%          eta = max min dist(x_i,y), i=1,...,N, y \in S^2
%
%   Here dist is the distance on the sphere between two points (i.e. great
%   circle distance).  The points are assumed to be on the unit sphere.  If
%   they are not then you need to multiply the results by the radius of the
%   sphere.
%
%   Note that in the literature that the covering radius is also
%   called the meshnorm in the meshfree methods literature. Also, note that
%   the separation of the point set is sometimes called the separation
%   radius, in which case it is stated as 1/2 of the value delta stated
%   above.

%   Author:  Grady Wright 2017

if isempty(X)
    eta = 0;
    delta = 0;
    return;
end

if size(X,2) ~= 3
    error('SPHEREPTS:SEPARATIONCOVERING',...
        'The input node set must have three columns.');
end

% Normalize to the unit sphere
X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));

% Separation radius is easy, it is the arccos of the dot product of the two
% closest points.  We will use a KD tree to find the nearest points in the
% set.
[nnidx,dist] = findKNearestNeighbors(X,X,2);
[~,id] = min(dist(:,2));
delta = acos(dot(X(id,:),X(nnidx(id,2),:)));

% Covering radius is a bit harder and computationally intensive, so only do
% it if the user request it as output
if nargout <= 1
    return;
end
    
% Triangulate the nodes.
tri = delaunaySph(X);

% Local variables to make the code easier to read.
x = X(:,1);
y = X(:,2);
z = X(:,3);

% Find the circumcenter on the sphere of each triangle
xc = ComputeCircumCenterSpherePoints(tri,X);

mn = zeros(length(tri),3);
% For each point, compute the (spherical) distances between the
% circumcenter and each of the vertices.
mn(:,1) = acos((x(tri(:,1)).*xc(:,1) + y(tri(:,1)).*xc(:,2) + z(tri(:,1)).*xc(:,3)));
mn(:,2) = acos((x(tri(:,2)).*xc(:,1) + y(tri(:,2)).*xc(:,2) + z(tri(:,2)).*xc(:,3)));
mn(:,3) = acos((xc(:,1).*x(tri(:,3)) + xc(:,2).*y(tri(:,3)) + xc(:,3).*z(tri(:,3))));

eta = max(mn(:));

end

% Author: Grady B. Wright
function [xc,cr] = ComputeCircumCenterSpherePoints(tri,X)

% Local variables to make the code easier to read.
x = X(:,1);
y = X(:,2);
z = X(:,3);

%
% See http://en.wikipedia.org/wiki/Circumcenter for the details of the
% calculation.
%

% For each point, compute the distances between each of the vertices.
rd2(:,1) = 2*(1-(x(tri(:,1)).*x(tri(:,2)) + y(tri(:,1)).*y(tri(:,2)) + z(tri(:,1)).*z(tri(:,2))));
rd2(:,2) = 2*(1-(x(tri(:,2)).*x(tri(:,3)) + y(tri(:,2)).*y(tri(:,3)) + z(tri(:,2)).*z(tri(:,3))));
rd2(:,3) = 2*(1-(x(tri(:,1)).*x(tri(:,3)) + y(tri(:,1)).*y(tri(:,3)) + z(tri(:,1)).*z(tri(:,3))));

alpha = rd2(:,2).*dot(X(tri(:,1),:)-X(tri(:,2),:),X(tri(:,1),:)-X(tri(:,3),:),2);
denom = 2*sum(cross(X(tri(:,1),:)-X(tri(:,2),:),X(tri(:,2),:)-X(tri(:,3),:),2).^2,2);
alpha = alpha./denom;

beta = rd2(:,3).*dot(X(tri(:,2),:)-X(tri(:,1),:),X(tri(:,2),:)-X(tri(:,3),:),2);
beta = beta./denom;

gamma = rd2(:,1).*dot(X(tri(:,3),:)-X(tri(:,1),:),X(tri(:,3),:)-X(tri(:,2),:),2);
gamma = gamma./denom;

xc = repmat(alpha,[1 3]).*X(tri(:,1),:) + repmat(beta,[1 3]).*X(tri(:,2),:) + repmat(gamma,[1 3]).*X(tri(:,3),:);

% Calculate the radius for each circumcenter point.
cr = prod(sqrt(rd2),2)./sqrt(denom);

%
% Project onto the sphere.
%

% Get the radius of the sphere the nodes lie on.
rr = sqrt(x(1)^2 + y(1)^2 + z(1)^2);

% Convert the centroids to spherical coordinates
[lam,th] = cart2sph(xc(:,1),xc(:,2),xc(:,3));

% Convert back to cartesian coordinates with with correct reference radius
[xc(:,1),xc(:,2),xc(:,3)] = sph2cart(lam,th,rr);

end
