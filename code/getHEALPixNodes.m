function x = getHEALPixNodes(S)
%GETHEALPIXNODES Computes the Hierarchical Equal Area iso-latitude
%   PIXelization nodes.
%
%   X = getHEALPixNodes(S) returns an N-by-3 matrix containing the HEALPix
%   nodes on the sphere, which were developed for applications related to
%   cosmic microwave background (CMB) radiation. For a given S, the total
%   number of nodes returned is
%                            N = 12*S^2;
%   The columns of X corresponds to the (x,y,z) cordinates of the nodes.
%   
%   For more details see
%   K. M. Gorski, E. Hivon, A. J. Banday, B. D. Wandelt, F. K. Hansen, M.
%   Reinecke, and M. Bartelmann. Healpix: A framework for high-resolution
%   discretization and fast analysis of data distributed on the sphere.
%   The Astrophysical Journal, 622(2):759-771, 2005.
%
%   Example:
%       x = getHEALPixNodes(10);  % Gives N = 1200 nodes
%       plotSphNodes(x);

% Author: Grady Wright, 2018

xn = zeros(2*S*(S-1),3); xs = xn;
count = 1;
for j=1:S-1
    k = (0:4*j-1)';
    % Longitude (0 <= lam < 2*pi)
    lam = (k+0.5)/(2*j)*pi;
    % Latitude  (0 <= th < pi)
    th = ones(size(k))*acos(1-1/3*(j/S)^2);
    % Cartesian coordinates for the rings in the northern hemisphere
    xn(count+k,:) = [cos(lam).*sin(th) sin(lam).*sin(th) cos(th)];
    % Cartesian coordinates for the rings in the southern hemisphere
    xs(count+k,:) = flipud(xn(count+k,:)).*repmat([1 1 -1],[4*j 1]);
    count = count + 4*j;
end

% Longitude (0 <= lam < 2*pi)
lam = linspace(0,2*pi,4*S+1); lam = lam(1:end-1); 
% Latitude  (0 <= th < pi)
th = acos(2/3-2*(0:2*S)/(3*S));
[tt,ll] = meshgrid(th,lam);
% Every other ring is offset by pi/(4*S) from zero in longitude
ll(:,1:2:end) = ll(:,1:2:end) + pi/(4*S);
ll = ll(:); tt = tt(:);
% Cartesian coordinates for the rings along the equator.
xe = [cos(ll).*sin(tt) sin(ll).*sin(tt) cos(tt)];

% Assemble the final point set
x = [xn;xe;flipud(xs)];

end



