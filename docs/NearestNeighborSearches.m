%% Nearest neighbor searching
% Grady Wright

%%
% In this tutorial we look at how to efficiently solve the following
% problem:
%
% Let $X=\{\mathbf{x}_j\}_{j=1}^N$ and $Y=\{\mathbf{y}_i\}_{i=1}^M$ be sets
% of points on the surface of the sphere. For each point $\mathbf{y}_i$,
% determine its $k$ nearest neighbors in $X$.
%
% Here we interpret nearest neighbors to mean the closest points as
% measured by Euclidean distance. 
%
% A Naive solution to this problem is to compute the distance for each
% point in $Y$ to all the points in $X$, then sort these distances to
% determine the $k$ nearest neighbors.  This computation requires 
% $\mathcal{O}(M N)$ operations, which is too expensive for large $N$ and
% $M$. 
% 
% A much more computationally efficient solution is to first build a
% KD-tree of the nodes in $X$ and then use it to search for the $k$ nearest
% neighbors in $Y$.  This operation requires $\mathcal{O}(N\log N)$
% operations to build the KD-tree and $\mathcal{O}(M\log N)$ operations to
% search for the nearest neighbors.
%
% There are functions for building and searching a KD-tree included in the
% MATLAB Statistics Toolbox.  Alternatively, if one does not have this
% package, then several KD-tree codes are available to download for free
% from MatlabCentral from the Mathworks website.  For example, the KD-tree
% package written by Andrea Tagliasacchi is particularly nice.
%
% The |findKNearestNeighbors| function in the *spherepts* provides an easy
% way to do nearest neighbor searching with a KD-tree using either
% functions from the Statistics Toolbox or from the Tagliasacchi package.
% The default is to use the former if available.  
LW = 'linewidth'; lw = 1; FS = 'FontSize'; fs = 12; 
FC = 'FaceColor'; fc = [0.93 0.93 0.93]; EC = 'EdgeColor'; ec = 'none';
MS = 'MarkerSize'; ms = 20; vw = [77 7];

%% Using the |findKNearestNeighbors| function
% Here we give an example for finding the $k$ nearest neighbors in a node 
% set $X$ on the sphere to a given point.  We use the minimum energy nodes
% as $X$ and the point $(1,0,0)$ as $Y$
x = getMinEnergyNodes(3136);
y = [1 0 0];
%%
% The 21 nearest neighbors in $X$ to $Y$ can be found as
idx = findKNearestNeighbors(x,y,21);
%%
% |idx| contains indices for the 21 nearest neighbors of |x| to |y|.  You
% can list these as 
x(idx,:)
%%
% The indices in |idx| are sorted according to the distance each point is
% in $X$ to $Y$. These distances can be obtained by calling
% |findKNearestNeighbors| with a second output argument:
[idx,dist] = findKNearestNeighbors(x,y,21);
dist'

%%
% To get a visual picture, we can plot all the nodes in $X$ and the node in
% $Y$, then highlight the nearest neighbors.
plotSphNodes(x); hold on;
plot3(y(:,1),y(:,2),y(:,3),'c.',MS,ms);       % Plot the node y in cyan
plot3(x(idx,1),x(idx,2),x(idx,3),'r.',MS,ms);
view(vw);
hold off;
%%
% More than one node can be included in $Y$ when using the
% |findKNearestNeighbors| function.  In this case, each row of idx will
% correspond to the indices in $X$ closest to the point in each row of $Y$.
% Here's an example:
y = [[1 0 0];[0 1 0];[0 0 1]];
idx = findKNearestNeighbors(x,y,5);
disp('5 nearest nodes to [1 0 0]');
x(idx(1,:),:)
disp('5 nearest nodes to [0 1 0]');
x(idx(2,:),:)
disp('5 nearest nodes to [0 0 1]');
x(idx(3,:),:)
%%
% To find the 20 nearest neighbors between all points in $X$ we simply use:
idx = findKNearestNeighbors(x,x,20);
%%
% This size of idx in this case is
size(idx)

%% Computing effective resolution of a node set
% The effective resolution of a node set on the sphere is typically defined
% as the maximum distance between all nearest neighbors of the point set.
% In this case the distance is measured as great circle distance. We can
% compute this quantity quite naturally using a KD-tree.
[idx,dist] = findKNearestNeighbors(x,x,2);
maxdist = max(dist(:,2));
sphdist = acos(1-maxdist^2/2);
fprintf('Effective resolution of N=3136 ME nodes on unit sphere is %f\n',sphdist);
%%
% To see what this means for Earth we simply need to scale this quantity by
% the mean radius of Earth: 6371km
R = 6371;
sphdist = R*acos(1-maxdist^2/2);
fprintf('Effective resolution of N=3136 ME nodes on Earth is %f km\n',sphdist);
%%
% This can be interpreted as meaning that there is one node every 235.7 km.
% Here is this quantity for a larger set of MD nodes.
N = 192^2;
x = getMaxDetNodes(N);
[idx,dist] = findKNearestNeighbors(x,x,2);
maxdist = max(dist(:,2));
sphdist = R*acos(1-maxdist^2/2);
fprintf('Effective resolution of N=%d MD nodes on Earth is %f km\n',N,sphdist);
%%
% Here is this quantity for an even larger set of Icosahedral nodes.
x = getIcosNodes(8,0);
N = size(x,1);
[idx,dist] = findKNearestNeighbors(x,x,2);
maxdist = max(dist(:,2));
sphdist = R*acos(1-maxdist^2/2);
fprintf('Effective resolution of N=%d icosahedral nodes on Earth is %f km\n',N,sphdist);
%%
% How many quasi-uniformly nodes are need for a resolution of 1km?  
