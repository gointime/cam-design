function [R_b,eps,phi,rho_p,x_s,y_s,x_c,y_c] = ...
    trans_roller(x,s,v,a,omega,R_f,eps_min,eps_max,R_b,envelope)
% RADIAL CAM WITH TRANSLATING ROLLER FOLLOWER
% omega duong neu quay nguoc chieu kdh

if nargin < 10,  envelope = -1; end % -1 for inner and +1 for outer
eps = eps_min;
eps_rp = []; %eps and R_p
if R_b ~= 0, R_b = R_b -0.05; end
while eps <= eps_max
    s0 = R_b+R_f;
    step_s0 = 0.05;
    while step_s0 > 0.01
        s0 = s0 + step_s0;
        R_p = sqrt(s0^2 + eps^2);
        phi = atan((v-eps)./(s + s0));
        if max(abs(phi))> 30/180*pi, continue; end
        rho_p = ((R_p + s).^2 + v.^2 ).^(3/2)...
            ./((R_p+s).^2 + 2*v.^2 - a.*(R_p + s));
        if min(abs(rho_p)) < 1.7*R_f, continue; end
        s0 = s0 - step_s0;
        step_s0 = step_s0/2;
    end
    eps_rp = [eps_rp;eps R_p]; %#ok<AGROW>
    eps = eps +0.1;
end
eps_rp = sortrows(eps_rp,2);
eps = round(eps_rp(1,1),1);
R_p = round(eps_rp(1,2),1);
R_b = R_p - R_f;

% Presure Angle, toi uu, khong dc > pi/6
phi = atan((v-eps)./ (s + sqrt(R_p^2 - eps^2)));
% Radius of Curvature
rho_p = ((R_p + s).^2 + v.^2 ).^(3/2)...
    ./((R_p+s).^2 + 2*v.^2 - a.*(R_p + s));
% Radius of Follower

% Roller Centerline Coordinates
lamda = atan(eps./(sqrt(R_p^2 - eps^2) + s)) + x;
x_f = sqrt((sqrt(R_p^2-eps^2)+s).^2 + eps^2).*cos(lamda);
y_f = -sign(omega)*sqrt((sqrt(R_p^2-eps^2)+s).^2 + eps^2).*sin(lamda);
% Cam Surface Coordinates
x_s = x_f + envelope*R_f.*cos(-phi+x); 
y_s = y_f - sign(omega)*envelope*R_f.*sin(-phi+x);
% Cam Cutter Coordinates for Manufaturer
R_c = 10; % Input
x_c = x_f + envelope*(R_f-R_c).*cos(-phi+x);
y_c = y_f - sign(omega)*envelope*(R_f-R_c).*sin(-phi+x);
