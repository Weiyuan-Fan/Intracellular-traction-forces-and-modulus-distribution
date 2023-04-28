% function: find the rotation angle between reference image and leastdeformed image.
% need to change the value for "I" as the name of the least deformed image
% need to change the value for "pattern" as the name of the reference image
% need to modify "threshold_l" so that dot pixels (bright pixels) are extracted,
% one may want to decide the value for "threshold_l" by running to line 31...
% and run imagesc(D1_l) to check 

% How to use:
% Draw a box in figure 3 around the center dot,one may want to imput "ctrl + c"
% in the command window to end this run and run again if the centroid of the 
% center dot is not accurate, but this accutally does not matter as long as we run
% the "static equilibrium" code after.
% the angle will be stored as "angle"

% 2023/3/28 - Weiyuan Fan

clear all
close all
clc

left = 100;
right = 800;
bottom = 350;
top = 1000;
% left and right are actually y direction
% bottom and top are actually x direction
I = imread('Noshift17.tif'); % the least deformed image
pattern = imread('Noshift18.tif'); % the reference image
pattern_2 = pattern(left:right,bottom:top); 
I_1 = double(filter2(ones(7,7)/49,I)); 
I_2 = I_1(left:right,bottom:top); % make I1 the same size as pattern_2 for corr calculation 
%% paremeters for reference image
threshold_r = 0; 
S_D_r = 0; 
D1_r = pattern_2>threshold_r; % extract the dot pixels
D2_r = D1_r; % a copy of D1 for further use
%% paremeters for least defomed image
D = 19; % dot diameter = 19
threshold_l = 315; % threshold for dot intensity is 311 
S_D_l = 1*100; % sum of the intensity for one dot
D1_l = I_2>threshold_l; % extract the dot pixels
%% identify each dot for dot image
i=1;
while (sum(sum(D1_r)))>0  % loop through dots in Image_1
    zero_1_d = zeros(size(D1_r,1),size(D1_r,2));
    [x1_d,y1_d] = find(D1_r>0); % position for all dot pixels

    rand_1_d = randi(size(x1_d,1)); % choose a random dot pixel

    a_1_d = x1_d(rand_1_d);             % position of the random dot pixel
    b_1_d = y1_d(rand_1_d);

    zero_1_d(a_1_d-D:a_1_d+D,b_1_d-D:b_1_d+D) = 1; % make a neighbor array of 1 within dot diameter 

    dot_image_1_d = and(zero_1_d,D1_r);  % Boolean And to extract one dot
    
    if sum(sum(dot_image_1_d))>S_D_r % Don't compute the centroid for a non-complete dot
        [x_d1_d,y_d1_d] = find(dot_image_1_d>0); % position for all dot pixels
        x_c1_d(i) = round(sum(sum(x_d1_d.*pattern_2(x_d1_d,y_d1_d)))/sum(sum(pattern_2(x_d1_d,y_d1_d)))); % centroid
        y_c1_d(i) = round(sum(sum(y_d1_d.*pattern_2(x_d1_d,y_d1_d)))/sum(sum(pattern_2(x_d1_d,y_d1_d))));
        % x_cl is the x position in the image
        % y_cl is the y position in the image
        i = i+1;
    end
    
        D1_r = D1_r-dot_image_1_d;        % subtract this dot from all dots 
 
end

x_d = (1:size(pattern_2,1))'; % y direction in the image
y_d = (1:size(pattern_2,2));  % x direction in the image
%% do the same for I1
i=1;
while (sum(sum(D1_l)))>0  % loop through dots in Image_1
    zero_1_r = zeros(size(D1_l,1),size(D1_l,2));
    [x1_r,y1_r] = find(D1_l>0); % position for all dot pixels

    rand_1_r = randi(size(x1_r,1)); % choose a random dot pixel

    a_1_r = x1_r(rand_1_r);             % position of the random dot pixel
    b_1_r = y1_r(rand_1_r);

    zero_1_r(a_1_r-D:a_1_r+D,b_1_r-D:b_1_r+D) = 1; % make a neighbor array of 1 within dot diameter 

    dot_image_1_r = and(zero_1_r,D1_l);  % Boolean And to extract one dot
    
    if sum(sum(dot_image_1_r))>S_D_l % Don't compute the centroid for a non-complete dot
        [x_d1_r,y_d1_r] = find(dot_image_1_r>0); % position for all dot pixels
        x_c1_r(i) = round(sum(sum(x_d1_r.*I_2(x_d1_r,y_d1_r)))/sum(sum(I_2(x_d1_r,y_d1_r)))); % centroid
        y_c1_r(i) = round(sum(sum(y_d1_r.*I_2(x_d1_r,y_d1_r)))/sum(sum(I_2(x_d1_r,y_d1_r))));
        % x_cl is the x position in the image
        % y_cl is the y position in the image
        i = i+1;
    end
    
        D1_l = D1_l-dot_image_1_r;        % subtract this dot from all dots 
 
end

x_r = (1:size(I_2,1))'; % y direction in the image
y_r = (1:size(I_2,2));  % x direction in the image
%% find the dot whose centroid is the cloest to the centroid of the whole image (dot image)
% centoid of the whole image
center_x_d = sum(sum(x_d.*D2_r(:,:)))/sum(sum(D2_r(:,:)));
center_y_d = sum(sum(y_d.*D2_r(:,:)))/sum(sum(D2_r(:,:)));

dis_d = sqrt((y_c1_d-center_y_d).^2+(x_c1_d-center_x_d).^2);
[val_d, index_d] = min(dis_d);

figure(2)
imagesc(pattern_2)
axis equal
hold on
scatter(y_c1_d(index_d),x_c1_d(index_d),'r')
hold on
scatter(center_y_d,center_x_d,'g')
axis equal
title('centroid of the whole dot image and the centriod of the center dot')
%% pick the center dot for real image
figure(3)
imagesc(I_2)
title('pick the center dot for real image')
axis equal
hold on
scatter(y_c1_r,x_c1_r,'r')
h_rect = imrect();
pos_rect = h_rect.getPosition();
pos_rect = round(pos_rect); 
for j = 1:size(x_c1_r,2)-1
    if (pos_rect(2)<=x_c1_r(j))&&(x_c1_r(j)<=pos_rect(2)+pos_rect(4))...
            &&(pos_rect(1)<=y_c1_r(j))&&(y_c1_r(j)<=pos_rect(1)+pos_rect(3))
        x_c1_r_c = x_c1_r(j);
        y_c1_r_c = y_c1_r(j);
    end
end
hold off
imagesc(I_2)
hold on
axis equal
scatter(y_c1_r_c,x_c1_r_c,'r')
%% location of the center dot in the uncoropped image
% location for the center dot of the uncoropped real image
center_y_r_u = y_c1_r_c+bottom-1;
center_x_r_u = x_c1_r_c+left-1;
center_y_d_u = y_c1_d(index_d)+bottom-1;
center_x_d_u = x_c1_d(index_d)+left-1;

figure(4)
subplot(1,2,1)
imagesc(I)
hold on
scatter(center_y_r_u,center_x_r_u,'r')
axis equal
title('center dot of the uncoropped real image')
subplot(1,2,2)
imagesc(pattern)
hold on
scatter(center_y_d_u,center_x_d_u,'r')
axis equal
title('center dot of the uncropped dot image')
%% pick the patch just big enough to cover the second ring 
% width = 47; % size of the pattern that is just big enough to cover the second ring 
% width = 87; % size of the pattern that is just big enough to cover the third ring
width = 200; 
ring_r = I_2(x_c1_r_c-width:x_c1_r_c+width,y_c1_r_c-width:y_c1_r_c+width);
figure(5)
imagesc(ring_r)
title('ring in real image')
ring_d = pattern_2(x_c1_d(index_d)-width:x_c1_d(index_d)+width,...
                   y_c1_d(index_d)-width:y_c1_d(index_d)+width);
                             
figure(6)
imagesc(ring_d)
title('ring in dot image')
%% rotate the image
subplot(1,2,1)
imagesc(ring_d)
axis equal
hold on
scatter((size(ring_d,1)+1)/2,(size(ring_d,1)+1)/2,'r')
subplot(1,2,2)
imagesc(ring_r)
axis equal
hold on
scatter((size(ring_d,1)+1)/2,(size(ring_d,1)+1)/2,'r')

theta = linspace(0,360,722);
for i = 1:size(theta,2)
    pattern_3 = rotateAroundapoint(ring_d,(size(ring_d,1)+1)/2, (size(ring_d,2)+1)/2, theta(i));
    DIC(i) = corr2(pattern_3,ring_r);
end
[corr, index_corr] = max(DIC);
pattern_final = rotateAroundapoint(ring_d,(size(ring_d,1)+1)/2, (size(ring_d,2)+1)/2,theta(index_corr));
angle = theta(index_corr);

figure(7)
subplot(1,3,1)
imagesc(ring_d)
colormap('gray')
axis equal
axis tight
title('rings from reference image before rotation')
subplot(1,3,2)
imagesc(pattern_final)
colormap('gray')
axis equal
axis tight
title('rings from reference image after rotation')
subplot(1,3,3)
imagesc(ring_r)
colormap('gray')
axis equal
axis tight
title('rings from least deformed image')
figure(8)
z = imfuse(pattern_final,ring_r);
imagesc(z)
axis equal
title('registration of reference image')