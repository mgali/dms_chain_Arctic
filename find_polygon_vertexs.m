% Finds the vertexs of a rectangle (e.g. rectangular macropixel in sinusoidal
% equal-area geogrid) given the following inputs:
%
%   xcenter = center of the rectangle, degrees longitude
%   ycenter = center of the rectangle, degrees latitude
%   xlength = base of the rectangle, km
%   ylength = height of the rectangle, km
%
% Vertexs are 4-point lon (xv) and lat (yv) vectors

function [xv,yv] = find_polygon_vertexs(xcenter,ycenter,xlength,ylength)

% % Uncomment for tests
% % xcenter = -179.763470000000;
% % ycenter = 45.1666680000000;
% xcenter = -179.290400000000;
% ycenter = 45.1666680000000;
% xlength = 4.64*8;
% ylength = 4.64*8;

deglatkm = 4.64*24; % length of 1 degree lat in km 
dlat = (ylength/2)/deglatkm; % lat distance from ycenter
yv = [ones(1,2)*(ycenter-dlat) ones(1,2)*(ycenter+dlat)];

% dlon in degrees is different at lower and upper side of rectangle
dlon1 = (xlength/2)/deglatkm/cos(deg2rad(ycenter-dlat));
dlon2 = (xlength/2)/deglatkm/cos(deg2rad(ycenter+dlat));
xv = [xcenter-dlon1 ones(1,2)*(xcenter+dlon2) xcenter-dlon1];