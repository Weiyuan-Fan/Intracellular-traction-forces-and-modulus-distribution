% Function: draw the boundary of a cell
% I_cell reads the brightfield image of a cell after deformation
crop = 50; % same as in the "registration.m"
I_cell = (imread('BVSMCs_p11_Fn_6p7kPa1_w2Brightfield Widefield_s11_t1.tif')); 
I_cell_1 = I_cell(crop:size(I_cell,1)-crop,crop:size(I_cell,2)-crop);
I_cell_1 = double(I_cell_1);
I_cell_2 = I_cell_1(center_x_r_u-R:center_x_r_u+R,center_y_r_u-R:center_y_r_u+R);  
boundary_points = boundary(I_cell_2,I2_o);
imagesc(I_cell_2)
axis equal
bbox = [100,100;900,900]; % draw a box for randomly distributed points
image = imfuse(I_cell_2,I2_o);
imagesc(image)
axis equal
hold on
scatter(y_fix,x_fix,'r')
