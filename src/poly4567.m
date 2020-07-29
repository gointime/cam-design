function [S,V,A,J] = poly4567(X,s1,s2)
% THE 4-5-6-7 POLYNOMIAL
if nargin < 3,error('Ham co it nhat 3 doi so.');end

h=abs(s2-s1);
beta = X(length(X))-X(1);

if s2-s1>0
    x = (X-X(1))/beta;
    S = s1+h*(35*x.^4-84*x.^5+70*x.^6-20*x.^7); 
    V = h/beta*(140*x.^3-420*x.^4+420*x.^5-140*x.^6);
    A = h/beta^2*(420*x.^2-1680*x.^3+2100*x.^4-840*x.^5);
    J = h/beta^3*(840*x-5040*x.^2+8400*x.^3-4200*x.^4);
else
    x = 1-(X-X(1))/beta;
    S = s2+h*(35*x.^4-84*x.^5+70*x.^6-20*x.^7); 
    V = -h/beta*(140*x.^3-420*x.^4+420*x.^5-140*x.^6);
    A = h/beta^2*(420*x.^2-1680*x.^3 +2100*x.^4-840*x.^5);
    J = -h/beta^3*(840*x-5040*x.^2+8400*x.^3-4200*x.^4);
end