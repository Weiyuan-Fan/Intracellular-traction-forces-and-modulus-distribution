% function: eliminate the image drift for DAPI images
% This has to be run right after  "registration.m"
% MAKE SURE TO NAME IMAGES SO THAT THEY HAVE SUFFIXES AS "01, 02...10, 11..."
% Don't forget to change the directory of the current folder
% 2023/5/2 - Weiyuan Fan

clearvars -except x_shift_total y_shift_total 
close all
clc

nimages = 25; % number of good images

imagefiles = dir('*.tif');      
nfiles = length(imagefiles);    % Number of files found
crop = 50;
for i=1:nimages
   currentfilename = imagefiles(i).name;
   currentimage = imread(currentfilename);
   images{i} = currentimage;
   images_crop{i} = currentimage(crop:size(currentimage,1)-crop,crop:size(currentimage,2)-crop);
end
images_n{1} = images_crop{1};
for i = 2:nimages
    currentimage = images{i};
    images_n{i} = currentimage(crop-x_shift_total(i):size(currentimage,1)-crop-x_shift_total(i),...
                               crop-y_shift_total(i):size(currentimage,2)-crop-y_shift_total(i));
end

 for i = 1:nimages
     file_name = sprintf('Nodrift_DAPI%d.tif', i);
     imwrite(images_n{i},file_name,'tif')
 end

