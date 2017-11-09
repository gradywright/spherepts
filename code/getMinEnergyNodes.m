function x = getMinEnergyNodes(N)
%GETMINENERGYNODES Returns a set of N quasi-minimum energy nodes on the sphere.
%
%   X = getMinEnergyNodes(N) returns an N-by-3 matrix containing
%   quasi-minimum Reisz energy nodes on the sphere.  The nodes are arranged
%   in such a way that Reisz energy of the nodes with respect to the power
%   2 is at a local minimum.  The energy is defined as
%   E(x) = sum_{i=1:N} sum_{j=i+1:N} 1 / ||  x_i - x_j  ||^2
%   where ||x_i-x_j|| is the Euclidean distance between node x_i and x_j.
%   The columns of matrix X corresponds to the (x,y,z) cordinates of the nodes.
%
%   Computing these nodes requires solving a non-linear optimization 
%   that can take significant time. This code does not solve this problem.
%   Instead these nodes were computed offline for several values of N.  
%   This code simply provides a way to access these nodes from MATLAB.
%
%   IMPORTANT NOTE: Presently the only nodes available are for
%         1) Every N from 3 to 3000 
%         2) Every N that is a perfect square between 3000 and 10201
%         3) The additional values of N=(128)^2=16,384 and N=(192)^2=36,864
%         4) The additional values of N=40,000, N=62,500, N=90,000, 
%            N=160,000, and N=250,000
%
%   The nodes from 2) and 3) above were computed by Prof. Rob Womersley and
%   made available at his website:
%       http://web.maths.unsw.edu.au/~rsw/Sphere/Energy/index.html
%   The nodes from 4) above were computed from a code by Prof. Doug Hardin
%   using the ideas from the paper:
%   S. V. Borodachov, D. P. Hardin, and E. Saff, Low complexity methods for
%   discretizing manifolds via Riesz energy minimization. Foundations of
%   Computational Math, Volume 14, Issue 6, pp 1173-1208.
%   These nodes differ slightly in that they use a power 3 in the Reisz
%   energy functional rather than a power of 2.
%
%   For more information on minimum energy nodes for the sphere see
%   D.P. Hardin, E.B. Saff, Discretizing manifolds via minimum energy
%   points, Notices Amer. Math. Soc. 51 (2004) 1186?1194.
%
%   Example 1:
%       x = getMinEnergyNodes(64^2);  % Retuns 4096 nodes on the sphere.
%       plotSphNodes(x);

% Available values of N, assuming the package was installed correctly.
possibleN = unique([3:4770 (2:101).^2 128^2 192^2]);
possibleN = sort(possibleN);

if isempty(find(N == possibleN, 1))
    id1 = find(N < possibleN,1);
    id2 = max(id1 - 1,1);
    if ~isempty(id1)
        error('SPHEREPTS:GETMINENERGYNODES',...
            'The value N=%d is not available, the next smallest available value is %d and the largest is %d',N,possibleN(id2),possibleN(id1));
    else
        error('SPHEREPTS:GETMINENERGYNODES',...
            'The value N=%d is not available, the closes available value is %d',N,possibleN(end));
    end
end

fname = sprintf('me%05d.mat',N);

if exist(fname,'file') == 0
    error('SPHEREPTS:GETMINENERGYNODES',...
        'Cannot find the minimum energy node file, please insure these files are installed on the MATLAB path');
end

x = [];
load(fname,'x')

end

