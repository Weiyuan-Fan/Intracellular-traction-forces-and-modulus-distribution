% function: eliminate the image drift for DAPI images
% This has to be run right after  "registration.m"
% MAKE SURE TO NAME IMAGES SO THAT THEY HAVE SUFFIXES AS "01, 02...10, 11..."
% Don't forget to change the directory of the current folder
% 2023/5/4 - Weiyuan Fan

clearvars -except x_shift_total y_shift_total 
close all
clc

nimages = 17; % number of good images

imagefiles = dir('*.tif');      
nfiles = length(imagefiles);    % Number of files found
crop = 50;
for i=1:nimages
   currentfilename = imagefiles(i).name;
   currentimage = imread(currentfilename);
   currentimage1 = double(imread(currentfilename));
   currentimage2 = currentimage1(50:100,50:100);
   average = mean(currentimage2(:));
   rand_y = (average-10)+20*rand(size(currentimage,1),crop);
   rand_x = (average-10)+20*rand(crop,size(currentimage,2)+2*crop);
   images{i} = currentimage;
   images_crop{i} = [rand_x;rand_y,currentimage,rand_y;rand_x];
end
images_n{1} = images{1};
for i = 2:nimages
    currentimage2 = images_crop{i};
    images_n{i} = currentimage2(crop-x_shift_total(i)+1:size(currentimage,1)+crop-x_shift_total(i),...
                               crop-y_shift_total(i)+1:size(currentimage,2)+crop-y_shift_total(i));
end

 for i = 1:nimages
     file_name = sprintf('Nodrift_DAPI%d.tif', i);
     imwrite(images_n{i},file_name,'tif')
 end

