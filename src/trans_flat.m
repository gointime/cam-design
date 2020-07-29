function [R_b,rho_p,x_s,y_s,x_c,y_c] = trans_flat(x,s,v,a,omega,R_b)
% RADIAL CAM WITH TRANSLATING FLAT-FACED FOLLOWER

rho_p = 0;
if R_b ~= 0, R_b = R_b -0.1; end
while min(rho_p)<1
    % Radius of Curvature
    R_b = R_b + 0.1;
    rho_p = R_b + s + a;
end

% Cam Surface Coordinates
sigma = x + atan(v./(R_b + s));
x_s = sqrt((R_b + s).^2 + v.^2).*cos(sigma);
y_s = -sign(omega)*sqrt((R_b + s).^2 + v.^2).*sin(sigma);
% Cam Cutter Coordinates for Manufaturer
R_c = 1; %radius of cutting tool
sigma_p = atan(v./(R_b + s + R_c)) + x;
x_c = sqrt((R_b + s + R_c).^2 + v.^2).*cos(sigma_p);
y_c = -sign(omega)*sqrt((R_b + s + R_c).^2 + v.^2).*sin(sigma_p);