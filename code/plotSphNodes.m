function [varargout] = plotSphNodes(x,mrkrsize,mrkspec)
%PLOTSPHNODES Makes a simple plot of a given set of nodes on the surface of
%   the unit sphere.
%
%   plotSphNodes(X) plots the nodes contained in X on the surface of the
%   unit sphere.  Here X is assumed to be an N-by-3 array with each row
%   consisting of the (x,y,z) Cartesian coordinates for a node.
%
%   cax = plotSphNodes(X) returns a handle to Axes for the plot.
%
%   Example:
%       x = 2*rand(101,3)-1;
%       x = bsxfun(@rdivide,x,sqrt(sum(x.^2,2)));  % Project to the sphere.
%       plotSphNodes(x);
%

% Author: Grady Wright, 2014

if ( nargin < 3 )
    mrkspec = 'b.';
    if ( nargin < 2 )
        mrkrsize = 6;
    end
end

% Make sure the nodes are on the unit sphere.
x = bsxfun(@rdivide,x,sqrt(sum(x.^2,2)));

%
% Generate a unit sphere.
%
[xx,yy,zz] = sphere(101);
% Color of the sphere will be yellow:
clr = [255 255 102]/255;
clr = [1 1 1]*0.90;
% Plot the sphere
cax = newplot;
hold on;
surf(xx,yy,zz,1+0*xx,'EdgeColor','None','FaceColor',clr);

%
% Add the nodes as black solid circles.
%
plot3(x(:,1),x(:,2),x(:,3),mrkspec,'MarkerSize',mrkrsize);
hold off;
daspect([1 1 1]);
view(3);

if nargout == 1
    varargout = cax;
end

end



