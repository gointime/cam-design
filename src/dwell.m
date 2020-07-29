function [S,V,A,J] = dwell(X,h)
    S = ones(size(X))*h;
    V = S*0;
    A = S*0;
    J = S*0;