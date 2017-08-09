function [ X, tri ]=getGenSpiralNodes(n)
%% Generates generalized n spiral points on the sphere. Output is an N x 3 
% matrix. 
% 
% [X,tri] returns a triangulation of the nodes.
%
% Author: T. Michaels
% 
% 
% [1] R. Bauer Distribution of Points on a Sphere with Application to Star
% Catalgos. J. Guid. Cont. Dyn., 23(1) 130-137, 2000

%% 

X=zeros(n,3); Theta=zeros(1,n); Phi=zeros(1,n);

H=(-1+1/n):(2/n):(1-1/n);
L=sqrt(n*pi);

Theta(1,:)=acos([H(1,:)]);
Phi(1,:)=L*Theta(1,:);

X(:,1)=sin(Theta(1,:)).*cos(Phi(1,:));
X(:,2)=sin(Theta(1,:)).*sin(Phi(1,:));
X(:,3)=H(1,:);

%Triangulate the nodes
tri = delaunay(X);
tri = freeBoundary(TriRep(tri,X));

end