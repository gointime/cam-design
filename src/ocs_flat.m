function [R_b,min_l,beta,c,rho_p,x_s,y_s,x_c,y_c] = ocs_flat(x,s,v,a,omega,c,f)
% RADIAL CAM WITH OSCILLATING FLAT-FACED FOLLOWER
% Offset f of the follower face from its centerline is positive if directed
% away from the cam center and negative if toward the cam center.
% xac dinh khoang cach tam cam va tam tay quay c

      % -1 neu khop quay nam ben phai truc cam va nguoc lai      
rho_p=0;
%if R_b ~= 0, R_b = R_b -0.1; end
beta = 0;
while min(rho_p) < 0.2 && beta < 70/180*pi
    
    beta = beta + 0.001;
    rho_p = a.*cos(beta+s)./(1+v).^3 + (1+2*v).*sin(beta+s)./(1+v).^2 + f/c;
        
end

c_c = c + f./sin(beta + s);
R_b = c_c(1)*sin(beta);

v_c = v.*(c*cos(beta+s) + f./tan(beta+s) - f./tan(beta+s).^2);
l = c_c.*cos(beta+s)-v_c;
m = sqrt((c_c.*sin(beta + s)).^2 + v_c.^2);
psi = pi/2 - (beta+s) - sign(omega).*asin(v_c./m);
min_l = max(l);
% Cam Surface Coordinates
x_s = m.*cos(psi - sign(omega)*x);
y_s = m.*sin(psi - sign(omega)*x);
% Cam Cutter Coordinates for Manufaturer
R_c = 2; %radius of cutting tool
sigma = psi + asin(v.*l./m);
x_c = x_s + R_c.*cos(sigma);
y_c = y_s + R_c.*sin(sigma);
