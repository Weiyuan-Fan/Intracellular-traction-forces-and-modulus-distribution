% function to compute the signed distantce
function dis = fd(points,boundary_points)
points = points';
x = points(1,:);
y = points(2,:);
x_b = boundary_points(1,:);
y_b = boundary_points(2,:);
x_b_c = [boundary_points(1,:),boundary_points(1,1)]; % closed boundary points
y_b_c = [boundary_points(2,:),boundary_points(2,1)];
in = inpolygon(x,y,boundary_points(1,:),boundary_points(2,:));
d = zeros(1,size(x_b,2));
for i = 1:size(points,2)
    for j = 1:size(x_b,2)
        a = [x_b_c(j+1)-x_b_c(j),y_b_c(j+1)-y_b_c(j)];
        b = [x_b_c(j)-x(i),y_b_c(j)-y(i)];
        c = [x_b_c(j+1)-x(i),y_b_c(j+1)-y(i)];
        if dot(a,b)*dot(a,c)<0
            A = y_b_c(j)-y_b_c(j+1);
            B = x_b_c(j+1)-x_b_c(j);
            C = x_b_c(j)*y_b_c(j+1)-x_b_c(j+1)*y_b_c(j);
            d(j) = abs(A*x(i)+B*y(i)+C)/sqrt(A^2+B^2);
        else
            d(j) = min(sqrt(b(1)^2+b(2)^2),sqrt(c(1)^2+c(2)^2));
        end 
    end
    dis(i) = min(d);
    d = zeros(1,size(x_b,2)); 

end

dis(in) = -dis(in);
dis = dis';
end

