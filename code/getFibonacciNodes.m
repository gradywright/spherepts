function x = getFibonacciNodes(N)
%% GETFIBONACCINODES 
% Computes Fibonacci (or phyllotaxis) node sets for the surface of the sphere.
%
% X = getFibonacciNodes(N) returns an N-by-3 matrix containing the 
%   Fibonacci nodes, where N must be odd. The columns corresponds
%   to the (x,y,z) cordinates of the nodes.  Note that these nodes are only
%   unique up to a roataion.
%
% There are different ways to generate the Fibonacci nodes. This code uses the 
% method advocated by Swinbank and Pruser:
%
% R. Swinbank and R. James Purser, Fibonacci grids: A novel approach to
% global modelling, Quarterly Journal of the Royal Meteorological Society,
% 132 (2006), pp. 1769-1793.
% 
% This construction avoids placing any nodes at the poles.

% Author: Grady Wright, 2014

N = ceil((N-1)/2);  % Has to be an odd number of points.
x = zeros(2*N+1,3);
lat = zeros(2*N+1,1);
lon = lat;
gr = (1+sqrt(5))/2;
k = 1;
for i = -N:N
   lat(k) = asin(2*i/(2*N+1));
   lon(k) = 2*pi*i/gr;
   k = k+1;
end
[x(:,1),x(:,2),x(:,3)] = sph2cart(lon,lat,1+0*lat);

end