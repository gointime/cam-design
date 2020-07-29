function [S,V,A,J] = cycloidal(X,s1,s2)
% CYCLOIDAL RISE ACCELERATION

if nargin < 3,error('Ham co it nhat 3 doi so.');end

h = abs(s2-s1);
beta = X(length(X))-X(1);

if s2-s1>0
    x = (X-X(1))/beta;
    S = s1+h*(x-1/(2*pi)*sin(2*pi*x));
    V = h/beta*(1-cos(2*pi*x));
    A = h/beta^2*2*pi*sin(2*pi*x);
    J = h/beta^3*4*pi^2*cos(2*pi*x);
else
    x = 1-(X-X(1))/beta;
    S = s2+h*(x-1/(2*pi)*sin(2*pi*x)); 
    V = -h/beta*(1-cos(2*pi*x));
    A = h/beta^2*2*pi*sin(2*pi*x);
    J = -h/beta^3*4*pi^2*cos(2*pi*x);
end