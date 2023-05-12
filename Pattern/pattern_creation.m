% function: create a pattern image
% change value of m for different pattern sizes
% uncomment the last line for first time use, it will save the pattern as myfigure.tif
% 2023/4/28 - Weiyuan Fan

clear all
close all
clc

m = 7; % number of rings
n_x = 1344; % size of the pattern
n_y = 1024;
x = 1:n_x;
y = 1:n_y;
[x,y] = meshgrid(x,y);
pattern = zeros(n_x,n_y);
center_x = (n_x)/2;
center_y = (n_y)/2;
r = 6.2; % radius of the dot in pixel under current zoom ratio
R = r*6; % increment of the radius of the ring
%%
for i = 2:m
    for j = 1:6*(i-1)
        centroid_x(i,j) = center_x + (R*(i-1))*cos(2*pi/(6*(i-1))*j);
        centroid_y(i,j) = center_y + (R*(i-1))*sin(2*pi/(6*(i-1))*j);
    end
end
%%    
centroid_x(centroid_x==0) = [];
centroid_y(centroid_y==0) = [];
centroid_x = [centroid_x,center_x];
centroid_y = [centroid_y,center_y];
% centroids of dots in CAD drawing are not intergers due to trigonometric
% so that dots created are not symmetric 
% choose to use the centroids as they are or round them up 

% centers = [centroid_x(:),centroid_y(:)];
centers = [round(centroid_x(:)),round(centroid_y(:))];
%========================================================================
radii = r*ones(size(centers,1),1);
pattern = circle_mask([n_y,n_x],centers,radii);    
figure(1)
imagesc(pattern)    
title('pattern')
axis equal
%imwrite(pattern, 'myfigure.tif');
