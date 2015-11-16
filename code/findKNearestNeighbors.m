function [nnidx,dist] = findKNearestNeighbors(x,y,k)
%FINDNEARESTNEIGHBORS Finds the K nearest neighbors in a given node set
%   using a K-D tree.
%
%   NNIDX = findKNearestNeighbors(X,Y,K) finds the K nodes in X that are
%   nearest to the given node Y using the standard Euclidean distance as
%   the distance measure. X is assumed to be a N-by-d array, where N is the
%   total number of nodes and d is the dimension of the Euclidean space the
%   nodes are in. Y can be a single 1-by-d array or a M-by-d array of
%   nodes.  NNIDX is an M-by-K array with row J containing the K indices of
%   the nodes in X that are closest to the Jth node in Y.  That is,
%   X(NNIDX(j,:),:) gives the K nearest nodes in X to the node Y(J,:).
%
%   [NNIDX,DIST] = findKNearestNeighbors(X,Y,K) returns the distances from
%   the node Y to the K nearest nodes in X.  DIST is an M-by-K array.
%
%   Note that this function is just a wrapper for other existing functions
%   in MATLAB that build K-D tree and search it for nearest neighbors. If
%   the MATLAB Statistics Toolbox is present then the code uses 
%   functions from this toolbox to do the searching.  Otherwise the code
%   uses the KD Tree package from Andrea Tagliasacchi which is made freely
%   available on the mathworks website.  Note that the user will need to
%   install this package, compile the corresponding mex files, and add
%   the location of the package to the MATLAB path.  This package can be
%   downloaded at http://www.mathworks.com/matlabcentral/fileexchange/21512-kd-tree-for-matlab
%
%   Example:
%    x = 2*rand(500,2)-1;  % 500 random nodes in [-1,1]x[-1,1]
%    % Find the 10 nearest neighbors to each of these 500 nodes
%    nnidx = findKNearestNeighbors(x,x,10);
%    plot(x(:,1),x(:,2),'b.');
%    hold on;
%    % Plot the 10 nearest neighbors to the first node in x.
%    plot(x(nnidx(1,:),1),x(nnidx(1,:),2),'ro','MarkerSize',8);
%    plot(x(1,1),x(1,2),'y.','MarkerSize',8);
%    hold off;

% Author: Grady Wright, 2014

%   
% Use the statistics toolbox KD Tree searching feature if the user has a
% license.
%
if license('test','statistics_toolbox')
    treeroot = createns(x);
    [nnidx,dist] = knnsearch(treeroot, y, 'k', k);
else  % Use the KD Tree code from by Andrea Tagliasacchi
    if exist('kdtree_build.m','file') == 0
        error('RBFSPHERE:FINDKNEARESTNEIGHBORS',...
            'Cannot find the KD Tree package from Andrea Tagliasacchi on the path. This package is necessary for this function.');
    end
    nnidx = zeros(size(y,1),k);
    tree = kdtree_build(x);
    for j=1:length(y)
        nnidx(j,:) = kdtree_k_nearest_neighbors(tree,y(j,:)',k);
    end
    %  This version of the KD tree sorts the nodes from farthest distance
    %  to closest distance.  So we reverse this.
    nnidx = fliplr(nnidx);
    kdtree_delete(tree);
    if nargout == 2  % Compute the distances.
        dist = nnidx;
        for j=1:length(y)
            dist(j,:) = sqrt(sum((y(j,:)-x(nnidx(j,:),:)).^2));
        end
        dist = fliplr(dist); % See comment above for flipping nnxidx
    end
end

end


