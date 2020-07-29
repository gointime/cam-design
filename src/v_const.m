function [S, V, A, J] = v_const(X,s1,s2)
 % VELOCITY CONSTANT
    if nargin < 3,error('Ham co it nhat 3 doi so.');end
    alpha = (s2-s1)/(X(end)-X(1));
    beta = s1-alpha*X(1);
    S = alpha*X+beta;
    V = alpha*ones(size(X));
    A = X*0;
    J = X*0;
end