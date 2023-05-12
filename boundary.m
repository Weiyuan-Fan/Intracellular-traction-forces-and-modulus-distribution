% BC Code
% BC stands for Boundary Click...
function [boundary_points] = boundary(image1,image2)
image = imfuse(image1,image2);
figure
imagesc(image);     % display image
xold = 0;
yold = 0;
k = 0;
hold on;           % and keep it there while we plot
while 1
    [xi, yi, button] = ginput(1);      % get a point
    if ~isequal(button, 1)             % stop if not left click
        break
    end
    k = k + 1;
    boundary_points(1,k) = round(xi);
    boundary_points(2,k) = round(yi);
    if xold                                 % the point link process excludes the foirst point
        plot([xold xi], [yold yi], 'ro-');  % link the current point with the previous one 
    else
        plot(xi, yi, 'ro');         
    end
    xold = xi;
    yold = yi;
  end
hold off;
end

