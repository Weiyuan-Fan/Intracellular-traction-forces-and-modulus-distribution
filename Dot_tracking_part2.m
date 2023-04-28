% function: dot tracking using combination of dot identification + DIC
% It tells the displacement of the dots from the reference image to the
% most deformed image.
% run this after running Dot_tracking_part1.

% "R" has to be bigger than the radius of the pattern,may need to chage it for
% different patch sizes.

% "good_image" is the total number of images including the refernce image.
% may need to change "D", "threshold", and "S_D" from case to case.
% "m" is the movement limit for a dot.

% How to use:
% Draw a box around the dots you like to add besides the automatic selected ones.
% Draw a box around the dots you don't like. 
% The centoids shown in image 1 do not have to be accurate, we picked the
% dots in least deformed image first to tell which dots should be kept in 
% the reference image, the centroids in the reference image are the ones we
% finally use for computation

% 2023/4/28 - Weiyuan Fan
clearvars -except angle center_y_r_u center_x_r_u center_y_d_u center_x_d_u
close all
clc

R = 340; % Radius of the island 
imagefiles = dir('*Noshift*.tif');      
nfiles = length(imagefiles);    % Number of files found
for i=1:nfiles
   currentfilename = imagefiles(i).name;
   currentimage = imread(currentfilename);
   images{i} = currentimage;
end
good_image = 18; 

D = 20; % dot diameter = 18
threshold = 312; % threshold for dot intensity is 311
S_D = 1*100; % sum of the intensity for one dot

I1 = double(filter2(ones(7,7)/49,images{good_image-1}));
I_pattern = I1;
I1_o = double(filter2(ones(1,1)/1,images{good_image-1}));
I1_o = I1_o(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R); 
I1_p = I1(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);          % choose the part where the pattern is in the middle
I1_reference = I1(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);  
D1 = I1_p>threshold; % extract the dot pixels
%%
i=1;
while (sum(sum(D1)))>0  % loop through dots in Image_1
    zero_1 = zeros(size(D1,1),size(D1,2));
    [x1,y1] = find(D1>0); % position for all dot pixels

    rand_1 = randi(size(x1,1)); % choose a random dot pixel

    a_1 = x1(rand_1);             % position of the random dot pixel
    b_1 = y1(rand_1);

    zero_1(a_1-D:a_1+D,b_1-D:b_1+D) = 1; % make a neighbor array of 1 within dot diameter 

    dot_image_1 = and(zero_1,D1);  % Boolean And to extract one dot
    
    if sum(sum(dot_image_1))>S_D % Don't compute the centroid for a non-complete dot
        [x_d1,y_d1] = find(dot_image_1>0); % position for all dot pixels
        x_c1(i) = sum(sum(x_d1.*I1_p(x_d1,y_d1)))/sum(sum(I1_p(x_d1,y_d1))); % centroid
        y_c1(i) = sum(sum(y_d1.*I1_p(x_d1,y_d1)))/sum(sum(I1_p(x_d1,y_d1)));
        i = i+1;
    end
    
        D1 = D1-dot_image_1;        % subtract this dot from all dots 
  
end

x_fix = round(x_c1);
y_fix = round(y_c1);
m = 20; % movement limit of a dot

x_dis_total = 0;
y_dis_total = 0;

imagesc(I1_o)
axis equal
hold on
scatter(y_fix,x_fix,'r')
s1 = 'Yes';
answer = questdlg('Want to select dots manually?');
while strcmp(answer,s1)==1
    h_rect = imrect();
    pos_rect = h_rect.getPosition();
    pos_rect = round(pos_rect); 
    pattern = I1_p(pos_rect(2):pos_rect(2)+pos_rect(4),pos_rect(1):pos_rect(1)+pos_rect(3));
    y_c1(i) = sum(sum((pos_rect(1):pos_rect(1)+pos_rect(3)).*...
        I1_p(pos_rect(2):pos_rect(2)+pos_rect(4),pos_rect(1):pos_rect(1)+pos_rect(3))))...
        /sum(sum(I1_p(pos_rect(2):pos_rect(2)+pos_rect(4),pos_rect(1):pos_rect(1)+pos_rect(3))));
    x_c1(i) = sum(sum((pos_rect(2):pos_rect(2)+pos_rect(4)).*...
        I1_p(pos_rect(2):pos_rect(2)+pos_rect(4),pos_rect(1):pos_rect(1)+pos_rect(3))'))...
        /sum(sum(I1_p(pos_rect(2):pos_rect(2)+pos_rect(4),pos_rect(1):pos_rect(1)+pos_rect(3))));
    x_fix(i) = round(x_c1(i));
    y_fix(i) = round(y_c1(i));
    i = i+1;
    imagesc(I1_p)
    hold on
    scatter(y_fix,x_fix,'r')
    answer = questdlg('Want to select dots manunally?');
end
answer = questdlg('Want to unselect dots?');
while strcmp(answer,s1)==1
    h_rect = imrect();
    pos_rect = h_rect.getPosition();
    pos_rect = round(pos_rect); 
    %pattern = I1_p(pos_rect(2):pos_rect(2)+pos_rect(4),pos_rect(1):pos_rect(1)+pos_rect(3));
    for j = 1:size(x_fix,2)-1
        if (pos_rect(2)<=x_fix(j))&&(x_fix(j)<=pos_rect(2)+pos_rect(4))...
                &&(pos_rect(1)<=y_fix(j))&&(y_fix(j)<=pos_rect(1)+pos_rect(3))
            x_fix(j) = [];
            y_fix(j) = [];
            x_c1(j) = [];
            y_c1(j) = [];
        end
    end
imagesc(I1_p)
hold on
scatter(y_fix,x_fix,'r')
answer = questdlg('Want to unselect dots?');
end
imagesc(I1_o)
hold on
scatter(y_fix,x_fix,'r')
axis equal
pause(3)
%% rotate dot image
pattern = double(filter2(ones(7,7)/49,images{good_image}));
pattern_1 = pattern(center_x_d_u-R:center_x_d_u+R,center_y_d_u-R:center_y_d_u+R); 
pattern_2 = rotateAroundapoint(pattern_1,R+1, R+1, angle);
figure(3)
axis equal
imagesc(pattern_2)
title('dot image after rotation')
figure(4)
axis equal
h = imfuse(pattern_2,I1_p);
imagesc(h)
title('overlap dot image and real image')
image{good_image} = pattern_2;
%% identify each dot for dot image
threshold_d = 0; % threshold for dot intensity is 311 
S_D_d = 0; % sum of the intensity for one dot
D1_d = pattern_2>threshold_d; % extract the dot pixels
D2_d = D1_d; % a copy of D1 for further use
i=1;
while (sum(sum(D1_d)))>0  % loop through dots in Image_1
    zero_1_d = zeros(size(D1_d,1),size(D1_d,2));
    [x1_d,y1_d] = find(D1_d>0); % position for all dot pixels

    rand_1_d = randi(size(x1_d,1)); % choose a random dot pixel

    a_1_d = x1_d(rand_1_d);             % position of the random dot pixel
    b_1_d = y1_d(rand_1_d);

    zero_1_d(a_1_d-D:a_1_d+D,b_1_d-D:b_1_d+D) = 1; % make a neighbor array of 1 within dot diameter 

    dot_image_1_d = and(zero_1_d,D1_d);  % Boolean And to extract one dot
    
    if sum(sum(dot_image_1_d))>S_D_d % Don't compute the centroid for a non-complete dot
        [x_d1_d,y_d1_d] = find(dot_image_1_d>0); % position for all dot pixels
        x_c1_d(i) = sum(sum(x_d1_d.*pattern_2(x_d1_d,y_d1_d)))/sum(sum(pattern_2(x_d1_d,y_d1_d))); % centroid
        y_c1_d(i) = sum(sum(y_d1_d.*pattern_2(x_d1_d,y_d1_d)))/sum(sum(pattern_2(x_d1_d,y_d1_d)));
        % x_cl is the x position in the image
        % y_cl is the y position in the image
        i = i+1;
    end
    
        D1_d = D1_d-dot_image_1_d;        % subtract this dot from all dots 
 
end
figure(5)
imagesc(pattern_2)
axis equal
hold on
scatter(y_c1_d,x_c1_d,'r')
%% only pick those dots also being picked in the least deformed image
h = zeros(size(x_c1_d));
for i = 1:size(x_c1_d,2)
    r(i) = min(sqrt((x_fix-x_c1_d(i)).^2+(y_fix-y_c1_d(i)).^2));
    if r(i)>1.2*D
        h(i) = 1;
    end
end
k = find(h);
y_c1_d(k) = [];
x_c1_d(k) = [];
y_c1_d_copy = y_c1_d;
x_c1_d_copy = x_c1_d;
figure(6)
imagesc(pattern_2)
axis equal
hold on
scatter(y_c1_d,x_c1_d,'r')   
%% pick the center dot for dot image
figure(11)
imagesc(pattern_2)
title('pick the center dot for dot image')
hold on
scatter(y_c1_d,x_c1_d,'r')
h_rect = imrect();
pos_rect = h_rect.getPosition();
pos_rect = round(pos_rect); 
for j = 1:size(x_c1_d,2)-1
    if (pos_rect(2)<=x_c1_d(j))&&(x_c1_d(j)<=pos_rect(2)+pos_rect(4))...
            &&(pos_rect(1)<=y_c1_d(j))&&(y_c1_d(j)<=pos_rect(1)+pos_rect(3))
        x_c1_d_c = x_c1_d(j);
        y_c1_d_c = y_c1_d(j);
    end
end
hold off
imagesc(pattern_2)
hold on
scatter(y_c1_d_c,x_c1_d_c,'r')
%% displacement from perfect dot image to least deformed image 
x_dis_total = 0;
y_dis_total = 0;
I1 = double(filter2(ones(7,7)/49,images{good_image-1}));
I1 = I1(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);  
I0 = pattern_2;
figure(10)
subplot(1,2,1)
imagesc(I0)
subplot(1,2,2)
imagesc(I1)
n = 18;
for j = 1:size(y_c1_d,2)
    patch_1 = I0(round(x_c1_d(j)-D/2):round(x_c1_d(j)+D/2),round(y_c1_d(j)-D/2):round(y_c1_d(j)+D/2));
    patch_2 = I1(round(x_c1_d(j)-D/2-m):round(x_c1_d(j)+D/2+m),round(y_c1_d(j)-D/2-m):round(y_c1_d(j)+D/2+m));  
    ncc = normxcorr2(patch_1,patch_2);  
    % crop off the part where the patch is not completely inside
    ncc_2 = ncc(size(patch_1,1)+1:size(patch_1,1)+1+2*m,size(patch_1,1)+1:size(patch_1,1)+1+2*m);  
    max_ncc(j) = max(ncc_2(:));
    [xpeak(j),ypeak(j)] = find(ncc_2==max(ncc_2(:)));
    y_dot(n,j) = ypeak(j)-m+y_c1_d(j);  % the matching position for the DIC subset 
    x_dot(n,j) = xpeak(j)-m+x_c1_d(j);

    x_dis(n,j) = x_dot(n,j)-x_c1_d(j);
    y_dis(n,j) = y_dot(n,j)-y_c1_d(j);   
    y_c1_d(j) = y_dot(n,j);
    x_c1_d(j) = x_dot(n,j);
end
x_dis_total = x_dis_total+x_dis(n,:);
y_dis_total = y_dis_total+y_dis(n,:);

y_c1_d_f = y_c1_d;
x_c1_d_f = x_c1_d;
figure(11)
subplot(1,2,1)
imagesc(I0)
axis equal
hold on
scatter(y_c1_d_copy,x_c1_d_copy,'r')
subplot(1,2,2)
imagesc(I1)
axis equal
hold on
scatter(y_c1_d,x_c1_d,'r')
%% displacement from the least deformed image to the last image on the list 
% compute the DIC
axis equal
for n = (good_image-1):-1:2
    I1 = double(filter2(ones(7,7)/49,images{n}));
    I2 = double(filter2(ones(7,7)/49,images{n-1}));
    I1_o = double(filter2(ones(1,1)/1,images{n}));
    I1_o = I1_o(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);  
    I2_o = double(filter2(ones(1,1)/1,images{n-1}));
    I2_o = I2_o(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);  
    I1_p = I1(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);    
    I2_p = I2(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);   
    for j = 1:size(x_fix,2)
        patch_1 = I1_p(round(x_c1_d(j)-D/2):round(x_c1_d(j)+D/2),round(y_c1_d(j)-D/2):round(y_c1_d(j)+D/2));
        patch_2 = I2_p(round(x_c1_d(j)-D/2-m):round(x_c1_d(j)+D/2+m),round(y_c1_d(j)-D/2-m):round(y_c1_d(j)+D/2+m));  
        ncc = normxcorr2(patch_1,patch_2);  
        % crop off the part where the patch is not completely inside
        ncc_2 = ncc(size(patch_1,1)+1:size(patch_1,1)+1+2*m,size(patch_1,1)+1:size(patch_1,1)+1+2*m);  
        max_ncc(j) = max(ncc_2(:));
        [xpeak(j),ypeak(j)] = find(ncc_2==max(ncc_2(:)));
        y_dot(n,j) = ypeak(j)-m+y_c1_d(j);  % the matching position for the DIC subset 
        x_dot(n,j) = xpeak(j)-m+x_c1_d(j);
        
        x_dis(n,j) = x_dot(n,j)-x_c1_d(j);
        y_dis(n,j) = y_dot(n,j)-y_c1_d(j);   
        y_c1_d(j) = y_dot(n,j);
        x_c1_d(j) = x_dot(n,j);
    end
        x_dis_total = x_dis_total+x_dis(n,:);
        y_dis_total = y_dis_total+y_dis(n,:);
    
        figure(7)
        imagesc(I2_o)
        hold on
        scatter(y_c1_d,x_c1_d,'r')
        hold off
        axis equal
        pause(0.5)
end

%%
figure(9)
colormap('gray')
imagesc(I1_reference)
hold on; 
freezeColors;
colormap('jet');
quiver_colorbar(y_c1_d_f,x_c1_d_f,-y_dis_total*0.1613,-x_dis_total*0.1613)
axis ij
axis equal
title('displacement (μm)')

y_fix = y_c1_d_f;
x_fix = x_c1_d_f;
