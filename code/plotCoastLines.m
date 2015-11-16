function plotCoastLines(ptype,lspec,fspec)
%plotCoastLines Plots the coast lines of the land on Earth.
%
%   plotCoastLines(TYPE) Plots the coast lines of the land on Earth
%   according to the plot TYPE selected. The default is to plot a solid
%   black line and to not fill the continents in with a color.
%
%   TYPE=0
%   Plots the coast lines in longitude-latitude coordinates so that a 2D
%   figure results.
%
%   TYPE=1 (default)
%   Plots the coastlines on the surface of the unit sphere so that a 3D
%   figure results.
%
%   TYPE=2
%   Plots the coast lines using the Hammer map projection so that a 2D
%   figure results.
%
%   plotCoastLines(TYPE,LSPEC) Allows the properties of the lines
%   representing the coasts to be specified using LSPEC according to
%   MATLABs plot command.
%
%   plotCoastLines(TYPE,LSPEC,FSPEC) Fills the continents with a color
%   specified in FSPEC using a character string. FSPEC must be one of the
%   following values 'r','g','b', 'c','m','y','w','k', or an RGB row vector
%   triple, [r g b].
%
%   Example:
%   plotCoastLines(1,'c-','y');  % Plots coasts in cyan and continents in yellow
%   axis tight;
%   daspect([1 1 1]);

%   Author: Grady Wright, 2014

fillCont = false;
if nargin == 0
    ptype = 1;
    lspec = 'k-';
elseif nargin == 1
    lspec = 'k-';
elseif nargin == 3
    fillCont = true;
end

% Coast data contains all the information about how to draw the lines connecting
% the coasts of all the coninents.
load CoastData;

lon = ncst(:,1)*pi/180;    % Data is stored in degrees
lat = ncst(:,2)*pi/180;

id = find(isnan(lon));

if ptype == 0
    if fillCont
        for j=1:length(id)-1
            fid = id(j)+1:id(j+1)-1;
            fill(lon(fid),lat(fid),fspec,'EdgeColor','none');
            hold on;
        end
    end
    plot(lon,lat,lspec);
elseif ptype == 2
    [xh,yh] = sph2hammer(lon,lat);
    if fillCont
        for j=1:length(id)-1
            fid = id(j)+1:id(j+1)-1;
            fill(xh(fid),yh(fid),fspec,'EdgeColor','none');
            hold on;
        end
    end
    plot(xh,yh,lspec);
    hold on;
    % Also plot the outline of projection;
    th = linspace(-pi/2,pi/2,101);
    lam = -pi+0*th;
    [xh,yh] = sph2hammer(lam,th);
    plot(xh,yh,'k--');
    lam = pi+0*th;
    [xh,yh] = sph2hammer(lam,th);
    plot(xh,yh,'k--');
else
    [x,y,z] = sph2cart(lon,lat,1.001 + 0*ncst(:,1));
    if fillCont
        for j=1:length(id)-1
            fid = id(j)+1:id(j+1)-1;
            fill3(x(fid),y(fid),z(fid),fspec,'EdgeColor','none');
            hold on;
        end
    end
    plot3(x,y,z,lspec)
end

end
