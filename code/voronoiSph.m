function varargout = voronoiSph(x,varargin)
%voronoiSph Plots the Voronoi diagram for a set of nodes on the surface of
%   the sphere in 3D space.
%
%   voronoiSph(X) produces a plot of the Voronoi diagram.
%
%   voronoiSph(X,TRI) uses the Delaunay triangulation TRI instead of 
%   computing it internally.
%
%   H = voronoiSph(...,'LineSpec') plots the diagram with color and
%   linestyle specified and returns handles to the line objects created in
%   H. 
%
%   Example:
%       x = getIcosNodes(1,1);
%       voronoiSph(x,'b-');
%
% See also DELAUNAYSPH, voronoiDiagramSph

% Author: Grady Wright, 2014

N = size(x,1);

if N <= 2
    error('RBFSPHERE:voronoiSph:NothingToPlot','Not enough nodes for the Voronoi Diagram');
end

computeTri = true;
linespec = 'k-';
lw = 1;

if N <= 400
    ms = 20;
else
    ms = 15;
end

while ~isempty(varargin)
    if ischar(varargin{1})
        linespec = varargin{1};
    elseif isnumeric(varargin{1})
        tri = varargin{1};
        computeTri = false;
    end
    varargin(1) = [];
end

if computeTri
    tri = delaunaySph(x);
end

% Compute the Voronoi diagram
r = voronoiDiagramSph(x,tri);

% Get the radius of the sphere the nodes lie on.
rad = sqrt(x(1,1).^2 + x(1,2).^2 + x(1,3).^2);

%
% Generate a unit sphere.
%
[ll,tt] = meshgrid(linspace(-pi,pi,501),linspace(-pi/2,pi/2,251));
% Color of the sphere will be yellow:
clr = [255 255 102]/255;
% % Color of the sphere will be green:
% clr = [50 205 50]/255;
% Plot the sphere
cax = newplot;

surf(rad*cos(ll).*cos(tt),rad*sin(ll).*cos(tt),rad*sin(tt),1+0*ll,'EdgeColor','None','FaceColor',clr);
hold on;

m = 51;
t = linspace(0,1,m).';

plines = zeros(7*m*(N+1),3);
count = 1;
nodeclr = zeros(N,1);
for j=1:N
    v=r{j};
    dv = diff([v;v(1,:)],1);
    for k=1:size(dv,1);
        pline = bsxfun(@plus,t*dv(k,:),v(k,:));
        % Project to sphere.
        pline = rad*bsxfun(@rdivide,pline,sqrt(sum(pline.^2,2)));
        id = count:(count+size(pline,1)-1);
        plines(id,:) = pline;
        count = count + size(pline,1);
    end
    plines(count,:) = nan(1,3);
    count = count + 1;
    nodeclr(j) = size(v,1);
end
num_lines = find(sum(plines.^2,2)==0,1);
plines = plines(1:num_lines-1,:);
% Plot the Voronoi regions
plot3(plines(:,1),plines(:,2),plines(:,3),linespec,'LineWidth',lw);
axis([-rad rad -rad rad -rad rad]);
view(3), axis square

% Add the nodes
scatter3(x(:,1),x(:,2),x(:,3),ms,nodeclr,'filled');
caxis([min(nodeclr) max(max(nodeclr),6)]);
colormap(flipud(hot));

if nargout == 1
    varargout = cax;
end

hold off;

end