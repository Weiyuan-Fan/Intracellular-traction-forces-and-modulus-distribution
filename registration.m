% function: eliminate the image drift for fluorescent dot images
% MAKE SURE TO NAME IMAGES SO THAT THEY HAVE SUFFIXES AS "01, 02...10, 11..."
% 2023/5/4 - Weiyuan Fan

clear all
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
   currentimage_top = currentimage1(:,size(currentimage,2)-100:size(currentimage,2));
   currentimage_bottom = currentimage1(:,1:100);
   currentimage_left = currentimage1(1:100,:);
   currentimage_right = currentimage1(size(currentimage,1)-100:size(currentimage,1),:);
   average_top = mean(currentimage_top(:));
   average_bottom = mean(currentimage_bottom(:));
   average_left = mean(currentimage_left(:));
   average_right = mean(currentimage_right(:));
   rand_top = (average_top-2)+4*rand(crop,size(currentimage,2)+2*crop);
   rand_bottom = (average_bottom-2)+4*rand(crop,size(currentimage,2)+2*crop);
   rand_left = (average_left-2)+4*rand(size(currentimage,1),crop);
   rand_right = (average_right-2)+4*rand(size(currentimage,1),crop);
   images{i} = currentimage;
   images_f{i} = double(filter2(ones(7,7)/49,currentimage)); 
   images_crop{i} = [rand_top;rand_left,currentimage,rand_right;rand_bottom];
end
%%
x_shift = zeros(1,nimages);
y_shift = zeros(1,nimages);
x_shift_total = zeros(1,nimages);
y_shift_total = zeros(1,nimages);
images_n{1} = images{1};
for i = 2:nimages
    ncc = normxcorr2(images_f{i},images_crop{i-1});
    ncc_2 = ncc(size(images_f{i},1):size(ncc,1),size(images_f{i},2):size(ncc,2));  
    max_ncc = max(ncc_2(:));
    [xpeak,ypeak] = find(ncc_2==max(ncc_2(:)));
    y_shift(i) = ypeak-crop-1; 
    x_shift(i) = xpeak-crop-1;
    y_shift_total(i) = y_shift_total(i-1)+y_shift(i);
    x_shift_total(i) = x_shift_total(i-1)+x_shift(i);
    currentimage2 = images_crop{i};
    images_n{i} = currentimage2(crop-x_shift_total(i)+1:size(currentimage,1)+crop-x_shift_total(i),...
                               crop-y_shift_total(i)+1:size(currentimage,2)+crop-y_shift_total(i));
end

 for i = 1:nimages
     file_name = sprintf('Nodrift%d.tif', i);
     imwrite(images_n{i},file_name,'tif')
 end

