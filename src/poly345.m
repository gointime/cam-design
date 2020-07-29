function [S,V,A,J] = poly345(X,s1,s2)
% THE 3-4-5 POLYNOMIAL

if nargin < 3,error('Ham co it nhat 3 doi so.');end

h=abs(s2-s1);
beta = X(length(X))-X(1);

if s2-s1>0
    x = X/beta;
    S = s1+h*(10*x.^3-15*x.^4+6*x.^5);
    V = h/beta*(30*x.^2-60*x.^3+30*x.^4);
    A = h/beta^2*(60*x-180*x.^2+120*x.^3);
    J = h/beta^3*(60-360*x+360*x.^2);
else
    x = 1-(X-X(1))./beta;
    S = s2+h*(10*x.^3-15*x.^4+6*x.^5);
    V = -h/beta*(30*x.^2-60*x.^3+30*x.^4);
    A = h/beta^2*(60*x-180*x.^2+120*x.^3);
    J = -h/beta^3*(60-360*x+360*x.^2);
end