function [x,wghts] = getSymTDesignNodes(t)
%GETSYMTDESIGNNODES Returns a set of symmetric t-design nodes on the 
%   sphere computed by Rob Womersley.
%
%   X = getSymTDesignNodes(T) returns an N-by-3 matrix containing the 
%   a symmetric t-design nodes for integrating exactly all sphereical
%   harmonics of degree T.  Here is N will be given by 
%                     N = 0.5*T^2 + T + O(1)
%   The nodes are symmetric about the equator.  All the node sets were
%   computed by Rob Womersley of University of New South Wales and are
%   freely available for download from
%           http://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/
%   This code just provides a simple MATLAB interface.
%
%   The columns of X corresponds to the (x,y,z) cordinates of the nodes.
%   If f_i = f(x_i,y_i,z_i), i = 1,...,N, are samples of some function f at
%   the t-design nodes, then the integral of f can be approximated by
%   4*pi/N*sum(f).  In the case where f is a spherical harmonic of
%   degree T, the formula is exact.
% 
%   [X,W] = getSymTDesignNodes(T) Gives the nodes and the quadrature
%   weights for these nodes, which all take the value W_i = 1/(4*pi*N).
%
%   IMPORTANT NOTE: Presently only nodes for T=1,3,5,...,325 are available.
%   An error is returned if the node set is not available.
%
%   For more information on how the nodes are computed see:
%   Rob Womersley's website http://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/
%
%   Example 1:
%       [x,w] = getSymTDesignNodes(21);  % t-design nodes for t=21.
%       plotSphNodes(x);
%       Y = sphHarm(21,13,x(:,1),x(:,2),x(:,3));
%       intY = w'*Y

% Author: Grady Wright, 2016

% Available values of T, assuming the package was installed correctly.
possibleT = 1:2:325;

if isempty(find(t == possibleT, 1))
    id1 = find(t < possibleT,1);
    id2 = id1 - 1;
    if ~isempty(id1)
        error('RBFSPHERE:GETSYMTDESIGNNODES',...
            'The value t=%d is not available, the next smallest available value is %d and the largest is %d',t,possibleT(id2),possibleT(id1));
    else
        error('RBFSPHERE:GETSYMTDESIGNNODES',...
            'The value t=%d is not available, the closest available value is %d',t,possibleT(end));
    end
end

fname = sprintf('std%03d.mat',t);

if exist(fname,'file') == 0
    error('RBFSPHERE:GETSYMTDESIGNNODES',...
        'Cannot find the symmetric t-design node file, please insure these files are installed on the MATLAB path');
end

x = [];
wghts = [];
load(fname,'x')
load(fname,'wghts')

end

