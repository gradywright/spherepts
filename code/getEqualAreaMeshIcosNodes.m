function [X, tri] = getEqualAreaMeshIcosNodes(m,n)
%% Returns mesh icosahedral nodes with Snyder equal area gnomonic projection. 
%  [X] = getEqualAreaMeshIcosNodes(m,n) Will return total points T in a 
%  T x 3 matrix.
%
%   T = 10(m^2+mn+n^2) + 2
%
% [X,tri] returns a triangulation of the nodes.
%
%  Author: T. Michaels
%
% [1] D.L.D Caspar and A. Klug Physical principles in the construction of
% regular viruses. Cold Springs Harb. Symp. Quant. Biol. 27, 1962
% 
% [2] J.P. Snyder An equal area map projection for polyhedral globes.
% Cartographic 29(1): 10-21, 1992
%
% [3] T. Michaels Equidistributed Icosahedral Configurations on the Sphere,
% submitted


%% Generate the points on the north polar triangles. For each polar triangle,
%rotate the points around one of the non polar vertices twice to generate
%each of the ten equatorial triangles

X = [];

RotNP = [cos(2*pi/5),-sin(2*pi/5),0;
         sin(2*pi/5),cos(2*pi/5),0;
         0          ,0          ,1];

Roty = [1/sqrt(5),0, -2/sqrt(5);
        0        ,1,0          ;
        2/sqrt(5),0,1/sqrt(5)  ];

  Y = TriEqualAreaProj(m,n)';    
for k=0:4
 
 Rotz = [cos(pi/5+k*2*pi/5),sin(pi/5+k*2*pi/5),0;
         -sin(pi/5+k*2*pi/5),cos(pi/5+k*2*pi/5),0;
         0                  ,0                 ,1];
    
 RotVert = Rotz^(-1)*Roty^(-1)*RotNP*Roty*Rotz;

X = [X;(RotNP^k*Y)'];
X = [X;(RotVert*RotNP^k*Y)'];
X = [X;(RotVert^2*RotNP^k*Y)'];

%To generate south polar triangles, rotate base triangle around x
%axis by pi and shift around south pole by pi/5+k*2pi/5.

Rotx = [1,0,0;
        0,-1,0;
        0,0,-1];

X = [X;(Rotz*Rotx*Y)'];    
end

%Remove duplicate points
X = real(X);
X = arrayfun(@round,1e8*X);
X = unique(X,'rows');
X = X/1e8;

%Project back onto the sphere
X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2)));

%Triangulate the nodes
tri = delaunay(X);
tri = freeBoundary(TriRep(tri,X));

end
function [x] = TriEqualAreaProj(m,n)

%The function TriEqualAreaProjection generates the eqaul area projection on
%the face of the icosahedron defined by base triangle
%
% (0,0,1),
% (2/sqrt(5)cos(-pi/5),2/sqrt(5)sin(-pi/5),1/sqrt(5)), 
% (2/sqrt(5)cos(pi/5),2/sqrt(5)sin(pi/5),1/sqrt(5)).

if (n>m)
    k = n;
    n = m;
    m = k;
end

%The icosahedron with surface area 4pi has radius, r and edge length a.
side = sqrt(4*pi/(5*sqrt(3)));
circum = side*sin(2*pi/5); 

%Change of base matrix from (m,n) coordinates to Euclidean
A = [1, 1/2;
     0, sqrt(3)/2];

 
 a = A*[0;0];
 b = A*[m;n];
 c = A*[-n;m+n];
 

%Generate lattice on plane

e1 = [1,0];
e2 = [1/2,sqrt(3)/2];

L = ones(m+n+1,1)*e1;
L=[-n:m;-n:m]'.*L;

for i = 0:m+n
         L1 = L(1:m+n+1,:) + i*ones(m+n+1,1)*e2;
         L = [L;L1];
end

 [sizeL,~] = size(L);

 Mbc = (c(2,1)-b(2,1))/(c(1,1)-b(1,1));
 
 %Change of base matrix sending [m;n] to [1;0] and [-n;m+n] to [0;1]
 MN = 1/(m^2+m*n+n^2)*[m+n,n;
      -n,m];
  

T = circum*[2/sqrt(5)*cos(-pi/5),2/sqrt(5)*cos(pi/5);
     2/sqrt(5)*sin(-pi/5),2/sqrt(5)*sin(pi/5);
      1/sqrt(5)-1,1/sqrt(5)-1];

%Rot represents rotation by 2pi/3 around the center of the triangle
ctao = (5+2*sqrt(5))/sqrt(((5+sqrt(5))^2+(5+2*sqrt(5))^2));

Roty2 = [ctao,0,-sqrt(1-ctao^2);
        0,       1, 0       ;
        sqrt(1-ctao^2),0,ctao ];

Rotz2 = [-1/2,-sqrt(3)/2,0;
        sqrt(3)/2,-1/2,0 ;
        0        ,0   ,1];

Rot = Roty2^(-1)*Rotz2*Roty2;    
RotInv = Roty2^(-1)*Rotz2^(-1)*Roty2;
    
p=1;

%For each point in the lattice, test if it is in the triangle. If so,
%rotate it and map it using the function ThirdTriProj.
for j=1:sizeL
    %Determine which points lie in the triangle third defined by two
    %vertices and the face center. Due to precision errors in Matlab, round
    %the Boolean statements to a reasonable tolerance.
    rL = round(1e10*L(j,2));
    rMbc = round(1e10*(Mbc*(L(j,1)-b(1,1))+b(2,1)));
    rA = round(1e10*norm(L(j,:)-a'));
    rB = round(1e10*norm(L(j,:)-b'));
    rC = round(1e10*norm(L(j,:)-c'));
    if(rL <= rMbc && rA>=rC && rA>=rB)
        %The last column stores which side of the altitude from a the
       %lattice point sits on
        if (rB>=rC)
        XIcos = [(T*MN*A^(-1)*[L(j,1);L(j,2)]+[0;0;circum])',1];
        xy = ThirdTriProj(XIcos(1,:));
        x(p,:) = xy;
        x(p+1,:) = Rot*xy';
        x(p+2,:) = RotInv*xy';
        p=p+3;
        
        else
        XIcos = [(T*MN*A^(-1)*[L(j,1);L(j,2)]+[0;0;circum])',-1];
        xy = ThirdTriProj(XIcos(1,:));
        x(p,:) = xy;
        x(p+1,:) = Rot*xy';
        x(p+2,:) = RotInv*xy';
        p=p+3;
        
        end
    end
end
    function [x] = ThirdTriProj(x)
        %Applies the Icosahedral Equal Area Projection within the 
        %traingle between the face center and two vertices. Input is
        %a 1 x 4 matrix, with the first three columns representing the
        %Cartesian coordinates and the fourth column is + or - 1.
        %Output is a 1 x 3 matrix of Cartesian coordinates.
        
    side2 = sqrt(4*pi/(5*sqrt(3)));
    circum2 = side2*sin(2*pi/5);    
        
    xx = x(1,1);
    y = x(1,2);
    z = x(1,3);
    
    h = sqrt(3)/2*1/sin(2*pi/5)*(circum2-z)/(1-1/sqrt(5))-side2/sqrt(3);

    w = real(sqrt(xx^2+y^2+(circum2-z)^2-(h+side2/sqrt(3))^2));
    psi = h^2*sqrt(3)/2+pi/6;
    q = h*w/2+pi/2;
    ttheta = cos(q)/(cos(psi)/sin(pi/3)-sin(q));
    cg = cos(psi)/sin(pi/3);
    ctao2 = (5+2*sqrt(5))/sqrt(((5+sqrt(5))^2+(5+2*sqrt(5))^2));
    rottao = [ctao2,0,sqrt(1-ctao2^2);
              0,1,0;
              -sqrt(1-ctao2^2),0,ctao2];
 
    tdelta = sqrt(1-cg^2)*sqrt(1+ttheta^2)/cg;
    x = [1/sqrt(1+ttheta^2)*tdelta/sqrt(1+tdelta^2);
        x(1,4)*ttheta/sqrt(1+ttheta^2)*tdelta/sqrt(1+tdelta^2);
        1/sqrt(1+tdelta^2)];
    
    x = rottao*x;
    x = x';
    
%Adjust for cumulative rounding errors at midpoint of Icosahedral vertices.    
    if(w<=10^(-6) && abs(psi-pi/5)<=10^(-6))
   x = 1/(4/5*cos(pi/5)^2+1/5)^(1/2)*[2/sqrt(5)*cos(pi/5),0,1/sqrt(5)];
    end
    end

end