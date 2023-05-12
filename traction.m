% function: static equilibrium
% calculate the traction forces from displacements.
% find a set of forces that are the closest to the measured forces and
% also consistent with static equilibrium

% export locations and values of traction forces to a .txt file for later use in FEniCS 

% run this after Dot_tracking code

% 2023/4/28 - Weiyuan Fan


clearvars -except x_dis_total y_dis_total y_fix x_fix I2_o
close all
clc

G =2000; % pa, shear modulus
nu = 0.445; % Poisson's ratio
% T = (2*pi*G*a*u)/(2-nu)
% T = (pi*E*a)/(1+nu)(2-nu)
% a is the radius of the dot
a = 1e-3;% mm
ratio = 0.1613e-3;% mm

x_dis_total = x_dis_total*ratio; % change dis from pixel to mm
y_dis_total = y_dis_total*ratio;
T1 = (2*pi*G*a*x_dis_total)/(2-nu)*1e3; % change T from mN to nN
T2 = (2*pi*G*a*y_dis_total)/(2-nu)*1e3;

syms alpha1 alpha2 beta 
fex = sym('fex', [1 size(T1,2)]);
fey = sym('fey', [1 size(T1,2)]);

eqn1 = sum(fex) == 0;
eqn2 = sum(fey) == 0;
eqn3 = fex-T1-alpha1-beta*(-y_fix) == 0;
eqn4 = fey-T2-alpha2-beta*(x_fix) == 0;
eqn5 = sum(x_fix.*fey-y_fix.*fex) == 0;

[S] = solve(eqn1,eqn2,eqn3,eqn4,eqn5);
C = struct2cell(S);
A = [C{:}];

alpha1 = double(A(1));
alpha2 = double(A(2));
beta = double(A(3));

for i = 1:size(T1,2)
    feX(i) = double(A(3+i));
    feY(i) = double(A(3+size(T1,2)+i));
end

error_x = (feX-T1)./T1;
error_y = (feY-T1)./T1;

figure(1)
imagesc(I2_o)
title('traction force nN')
hold on
quiver(y_fix,x_fix,-feY,-feX,'r')
axis ij
axis equal

figure(2)
colormap('gray')
imagesc(I2_o)
hold on; 
freezeColors;
colormap('autumn');
quiver_colorbar(y_fix,x_fix,-feY,-feX)
hold on; 
freezeColors;
colormap('winter');
axis ij
axis equal
title('traction force nN')

y_ref = [y_fix,400];
x_ref = [x_fix,100];
yyy_ref = [feY,-20];
xxx_ref = [feX,0];
figure(111)
quiver(y_ref,x_ref,-yyy_ref,-xxx_ref,'r')
hold on
quiver(y_fix,x_fix,-T2, -T1,'b')
text(500,150,'1nN')
axis ij
axis equal
legend('with reference correction','without reference correction')
title('traction force nN')

%  export data to .txt
data = [(y_fix*1.613.*10.^(-1))',(x_fix*1.613.*10.^(-1))',(feY./1e3)',(feX./1e3)'];
writematrix(data,'celldata.txt','Delimiter','tab')

