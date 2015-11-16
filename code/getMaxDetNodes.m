function [x,wghts] = getMaxDetNodes(N)
%GETMAXDETNODES Returns the set of N maximal determinant nodes of Womersley and
%   Sloan.
%
%   X = getMaxDetNodes(N) returns an N-by-3 matrix containing the maximal
%   determinant nodes on the sphere computed by Womersley and Sloan. These
%   are nodes arranged on the sphere so that the determinant of a certain
%   Gram matrix related to spherical harmonics is maximized, which requires
%   solving a non-linear optimization problem. The columns of X corresponds
%   to the (x,y,z) cordinates of the nodes.
% 
%   These nodes have been computed for several values of N and put on the
%   web by Prof. Womersley at:
%        http://web.maths.unsw.edu.au/~rsw/Sphere/Extremal/New/index.html 
%   This code simply provides a way to access these nodes from MATLAB.
%
%   IMPORTANT NOTE: Presently only nodes N that are a perfect square
%   between 4 and 27556 are available, with the additional set of N=36864.
%   An error is returned if the node set is not available.
%
%   [x,wghts] = getMaxDetNodes(N) Also returns a set of quadrature weights
%   for the nodes.  These weights are provided in the node sets that are
%   downloadable from the web.
%
%   For more information on how the nodes are computed see:
%   I. H. Sloan and R. S. Womersley, Extremal systems of points and
%   numerical integration on the sphere, Advances in Computational
%   Mathematics 21 (2004) 107--125.
%
%   Example 1:
%       x = getMaxDetNodes(64^2);  % Retuns 4096 nodes on the sphere.
%       plotSphNodes(x);
%
%   Example 2:
%       [x,wghts] = getMaxDetNodes(64^2);  % Retuns 4096 nodes and weights.
%       area_sphere = wghts'*ones(64^2,1);

% Author: Grady Wright, 2014

% Available values of N, assuming the package was installed correctly.
possibleN = [2:166 192].^2;

if isempty(find(N == possibleN, 1))
    id1 = find(N < possibleN,1);
    id2 = id1 - 1;
    if ~isempty(id1)
        error('RBFSPHERE:GETMAXDETNODES',...
            'The value N=%d is not available, the next smallest available value is %d and the largest is %d',N,possibleN(id2),possibleN(id1));
    else
        error('RBFSPHERE:GETMINENERGYNODES',...
            'The value N=%d is not available, the closes available value is %d',N,possibleN(end));
    end
end

fname = sprintf('md%05d.mat',N);

if exist(fname,'file') == 0
    error('RBFSPHERE:GETMAXDETNODES',...
        'Cannot find the maximum determinant node file, please insure these files are installed on the MATLAB path');
end

x = [];
wghts = [];
load(fname,'x')

end

