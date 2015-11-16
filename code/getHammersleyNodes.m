function x = getHammersleyNodes(N)
%GETHAMMERSLEYNODES Comutes a Hammersley set of nodes on the unit sphere,
%   which are low-discrepancy sequences of nodes.
%
%   X = getHammersleyNodes(N) returns an N-by-3 matrix of Hammersley nodes
%   on the sphere, which form a low-discrepancy sequence for the sphere.
%   The columns of X corresponds to the (x,y,z) cordinates of the nodes.
%
%   For more details on these node sets see
%   J. Cui and W. Freeden. Equidistribution on the sphere. SIAM Journal on
%   Scientific Computing, 18(2):595?609.
%
%   Tien-Tsin Wong and Wai-Shing Luk and Pheng-Ann Heng, 1997, Journal of
%   Graphics Tools , vol. 2, no. 2, 1997, pp 9-24.
%
%   Example:
%       x = getHammersleyNodes(2000);
%       plotSphNodes(x);

% Author: Grady Wright, 2014


% This code uses vdcorput, which was created by Dimitri Shvorob.

t = vdcorput(N,2);

% If the statistical toolbox is available then this code can be used
% alternatively to vdcorput.
% p = haltonset(1);  
% t = net(p,N+1);

t = 2*t(1:end-1) - 1;
phi = 2*pi*((2*(1:N)-1)/2/N).';
x = [sqrt(1-t.^2).*cos(phi) sqrt(1-t.^2).*sin(phi) t];

function s = vdcorput(k,b)
% VDCORPUT   Base-b Van der Corput sequence, elements 0,..,k
% INPUTS   : k - maximum sequence index, non-negative integer
%            b - sequence base, integer exceeding 1
% OUTPUTS  : s - (k+1)*1 array, with s(i) storing element (i+1)
%                of base-b Van der Corput sequence
% AUTHOR   : Dimitri Shvorob
if k ~= floor(k)||(k < 0)
   error('Input argument "k" must be a non-negative integer')
end
if b ~= floor(b)||(b < 2)
   error('Input argument "b" must be a positive integer greater than 1')
end
s = zeros(k+1,1);
for i = 1:k
    a = basexpflip(i,b);
    g = b.^(1:length(a));
    s(i+1) = sum(a./g);
end    

function[a] = basexpflip(k,b) % reversed base-b expansion of positive integer k
j = fix(log(k)/log(b)) + 1;
a = zeros(1,j);
q = b^(j-1);
for ii = 1:j
   a(ii) = floor(k/q);
   k = k - q*a(ii);
   q = q/b;
end
a = fliplr(a);
end

end

end