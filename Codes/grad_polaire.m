function grad = grad_polaire(x,y)
r = sqrt(x^2+y^2);
theta = mod(atan2(y,x),2*pi);
b = sin(theta);
a = cos(theta);
grad = [a -b;b a]*[2*sin(2*theta/3)/(3*r^(1/3)); 2*cos(2*theta/3)/(3*r^(1/3))];
end