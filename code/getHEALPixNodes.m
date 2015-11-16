function x = getHEALPixNodes(P)
%GETHEALPIXNODES Computes the Hierarchical Equal Area iso-latitude
%   PIXelization nodes.
%
%   X = getHEALPixNodes(P) returns an N-by-3 matrix containing the HEALPix
%   nodes on the sphere, which were developed for applications related to
%   cosmic microwave background (CMB) radiation. For a given P, the total
%   number of nodes returned is
%                            N = 12*P^2;
%   The columns of X corresponds to the (x,y,z) cordinates of the nodes.
%   
%   For more details see
%   K. M. Gorski, E. Hivon, A. J. Banday, B. D. Wandelt, F. K. Hansen, M.
%   Reinecke, and M. Bartelmann. Healpix: A framework for high-resolution
%   discretization and fast analysis of data distributed on the sphere.
%   The Astrophysical Journal, 622(2):759?771, 2005.
%
%   Example:
%       x = getHEALPixNodes(10);  % Gives N = 1200 nodes
%       plotSphNodes(x);

% Author: Grady Wright, 2014

Nside = P;
Npix = 12*Nside^2;

x = zeros(Npix+1,3);
count = 1;
for p=0:Npix
    ph = (p+1)/2;
    i = floor(sqrt(ph-sqrt(floor(ph))))+1;
    if i >= Nside
        break;
    end
    j = p+1-2*i*(i-1);
    zp = 1 - i^2/(3*Nside^2);
    phi = pi/2/i*(j-1/2);
    xp = cos(phi).*cos(asin(zp));
    yp = sin(phi).*cos(asin(zp));
    x(count,:) = [xp yp zp];
    count = count + 1;
end

for p=0:Npix
    i = floor(p/4/Nside)+Nside;
    if i > 2*Nside
        break;
    end
    j = mod(p,4*Nside)+1;
    zp = 4/3-2*i/3/Nside;
    s = mod(i-Nside+1,2);
    phi = pi/2/Nside*(j-s/2);
    xp = cos(phi).*cos(asin(zp));
    yp = sin(phi).*cos(asin(zp));
    x(count,:) = [xp yp zp];
    count = count + 1;
end

x = x(1:count-1,:);
xf = flipud(x(1:Npix/2-2*Nside,:));
xf(:,3) = -xf(:,3);
x = [x;xf];

end



