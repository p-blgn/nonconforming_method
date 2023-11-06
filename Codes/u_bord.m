function val = u_bord(x,y)

r_carre = x^2+y^2;
theta = atan2(y,x);
theta = mod(theta,2*pi);
val = (r_carre^(1/3))*sin(2*theta/3);
end