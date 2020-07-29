function [R_b,beta,l,phi,rho_p,x_s,y_s,x_f,y_f] = ...
    ocs_roller(x,s,v,a,omega,R_f,c,envelope)
% RADIAL CAM WITH OSCILLATING ROLLER FOLLOWER
 
if nargin < 8,  envelope = -1; end      
l=c*0.5;
rot = -1;
l_beta = [];

while l <=c*1.5
    step_beta = 0.01; beta=30/180*pi;
    while step_beta > 0.005
        beta=beta+step_beta;
        if beta >= 65/180*pi, break; end
        delta = beta + s;  % s, v, a theo rad
        b=sqrt(l^2 + c^2 - 2*l*c.*cos(delta));                                  
        psi = acos((c^2 + b.^2 - l^2)./(2*c.*b));
        
        phi = (pi/2 - asin(c./b.*sin(delta)))*sign(omega) ...
            + atan(1./(b.^2./(l*c.*sin(delta).*v)...
            -(c^2 - b.^2 - l^2)./(2*c.*b.*sin(psi))));
        if max(abs(phi))> 30/180*pi, continue; end
        
        rho_p = (c^2.*sin(delta).^2)...
            ./(cos(phi).*(c.*sin(delta)-l.*cos(phi)...
            .*(a.*cos(phi) + v.*(v-1).*sin(phi))));
        if min(abs(rho_p))< 2*R_f, continue; end
        beta=beta-step_beta;
        step_beta = step_beta/2;
    end
    l = l +0.1;
    if beta >= 65/180*pi, continue; end
    beta = beta + step_beta*2;
    R_p = sqrt(l^2 + c^2 - 2*l*c*cos(beta)); 
    l_beta=[l_beta;l beta R_p]; %#ok<AGROW>
end
l_beta = sortrows(l_beta,3);
l=l_beta(1,1);
R_p = round(l_beta(1,3)+0.005,2);
R_b = R_p -R_f;
beta = acos((l^2 + c^2 - R_p^2)/(2*l*c));
delta = beta + s;  % s, v, a theo rad
b=sqrt(l^2 + c^2 - 2*l*c.*cos(delta));                                    
psi = asin(l./b.*sin(delta));
disp(b(1))
% Pressure Angle
lamda = atan(1./(b.^2./(l*c.*sin(delta).*v)...
    -(c^2 - b.^2 - l^2)./(2*c.*b.*sin(psi))));
phi = (pi/2 - asin(c./b.*sin(delta)))*sign(omega) + lamda;

% Radius of Curvature
rho_p = (c^2.*sin(delta).^2)...
    ./(cos(phi).*(c.*sin(delta)-l.*cos(phi)...
    .*(a.*cos(phi) + v.*(v-1).*sin(phi))));

% Roller Centerline Coordinates
x_f = b.*cos(-sign(omega).*x - rot.*psi);
y_f = b.*sin(-sign(omega).*x - rot.*psi);
% Cam Surface Coordinates
sigma = -sign(omega)*x - pi + psi + sign(omega)*lamda;
x_s = x_f - envelope*R_f.*cos(sigma); 
y_s = y_f - envelope*R_f.*sin(sigma);

% Cam Cutter Coordinates for Manufaturer
R_c = 0.2; % radius of cutting tool
x_c = x_f - envelope*(R_f - R_c).*cos(sigma);
y_c = y_f - envelope*(R_f - R_c).*sin(sigma);
% -1 for inner and +1 for outer                                      
% xac dinh chieu dai tay quay l, input
% gamma la vi tri ban dau hop c va truc x
% -1 neu khop quay nam ben phai truc cam va nguoc lai
% Prime cicle radius, toi uu, co gioi han
%beta = 0.570; 