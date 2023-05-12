clear all
close all
clc

nimages=17; % number of good images


imagefiles = dir('*.tif');      
nfiles = length(imagefiles);    % Number of files found
crop = 50;
for i=1:nimages
   currentfilename = imagefiles(i).name;
   currentimage = imread(currentfilename);
   images{i} = currentimage;
   images_crop{i} = currentimage(crop:size(currentimage,1)-crop,crop:size(currentimage,2)-crop);
end
x_shift = zeros(1,nimages);
y_shift = zeros(1,nimages);
x_shift_total = zeros(1,nimages);
y_shift_total = zeros(1,nimages);
images_n{1} = images_crop{1};
for i = 2:17
    ncc = normxcorr2(images_crop{i},images{i-1});
    ncc_2 = ncc(size(images_crop{i},1)+1:size(images{i-1},1),size(images_crop{i},2)+1:size(images{i-1},2));  
    max_ncc = max(ncc_2(:));
    [xpeak,ypeak] = find(ncc_2==max(ncc_2(:)));
    y_shift(i) = ypeak-crop+1; 
    x_shift(i) = xpeak-crop+1;
    y_shift_total(i) = y_shift_total(i-1)+y_shift(i);
    x_shift_total(i) = x_shift_total(i-1)+x_shift(i);
    currentimage = images{i};
    images_n{i} = currentimage(crop-x_shift_total(i):size(currentimage,1)-crop-x_shift_total(i),...
                               crop-y_shift_total(i):size(currentimage,2)-crop-y_shift_total(i));
end


%  for i = 1:17
%      file_name = sprintf('Noshift%d.tif', i);
%      imwrite(images_n{i},file_name,'tif')
%  end

