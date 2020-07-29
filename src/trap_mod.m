function [S,V,A,J] = trap_mod(X,s1,s2)
% MODIFIED TRAPEZOIDAL ACCELERATION
if nargin < 3, error('Ham co it nhat 3 doi so.'); end

S = zeros(1,length(X)); V = S; A = S; J = S;
h = abs(s2-s1);
beta = X(length(X))-X(1);

if s2-s1>0
    for i = 1 : length(X)
        x=(X(i)-X(1))/beta;
        if i>=1 && i<= (length(X)-1)/8+1
            S(i) = s1+h*(.38898448*x-.0309544*sin(4*pi*x)); 
            V(i) = h/beta*.38898448*(1-cos(4*pi*x));
            A(i) = h/beta^2*4.888124*sin(4*pi*x);
            J(i) = h/beta^3*61.425769*cos(4*pi*x);
        elseif i>(length(X)-1)/8+1 && i<=3/8*(length(X)-1)+1
            S(i) = s1+h*(2.44406184*x^2-.22203097*x+.00723407);
            V(i) = h/beta*(4.888124*x-.22203097);
            A(i) = h/beta^2*4.888124;
            J(i) = 0;
        elseif i>3/8*(length(X)-1)+1 && i<=5/8*(length(X)-1)+1
            S(i) = s1+h*(1.6110154*x-0.0309544*sin(4*pi*x-pi)-.3055077);
            V(i) = h/beta*(1.6110154-.38898448*cos(4*pi*x-pi));
            A(i) = h/beta^2*4.888124*sin(4*pi*x-pi);
            J(i) = h/beta^3*61.425769*cos(4*pi*x-pi);
        elseif i>5/8*(length(X)-1)+1 && i<=7/8*(length(X)-1)+1
            S(i) = s1+h*(-2.44406184*x^2+4.6660917*x-1.2292648);
            V(i) = h/beta*(-4.888124*x+4.6660917);
            A(i) = -h/beta^2*4.888124;
            J(i) = 0;
        elseif i>7/8*(length(X)-1)+1 && i<=length(X)
            S(i) = s1+h*(.6110154+.38898448*x+.0309544*sin(4*pi*x-3*pi));
            V(i) = h/beta*.38898448*(1+cos(4*pi*x-3*pi));
            A(i) = -h/beta^2*4.888124*sin(4*pi*x-3*pi);
            J(i) = -h/beta^3*61.425769*cos(4*pi*x-3*pi);
        end
    end
else
    for i = 1 : length(X)
        x=1-(X(i)-X(1))/beta;
        if i>=1 && i<= (length(X)-1)/8+1
            S(i) = s2+h*(.6110154+.38898448*x+.0309544*sin(4*pi*x-3*pi));
            V(i) = -h/beta*.38898448*(1+cos(4*pi*x-3*pi));
            A(i) = -h/beta^2*4.888124*sin(4*pi*x-3*pi);
            J(i) = h/beta^3*61.425769*cos(4*pi*x-3*pi);
        elseif i>(length(X)-1)/8+1 && i<=3/8*(length(X)-1)+1
            S(i) = s2+h*(-2.44406184*x^2+4.6660917*x-1.2292648);
            V(i) = -h/beta*(-4.888124*x+4.6660917);
            A(i) = h/beta^2*-4.888124;
            J(i) = 0;
        elseif i>3/8*(length(X)-1)+1 && i<=5/8*(length(X)-1)+1
            S(i) = s2+h*(1.6110154*x-.0309544*sin(4*pi*x-pi)-.3055077);
            V(i) = -h/beta*(1.6110154-.38898448*cos(4*pi*x-pi));
            A(i) = h/beta^2*4.888124*sin(4*pi*x-pi);
            J(i) = -h/beta^3*61.425769*cos(4*pi*x-pi);
        elseif i>5/8*(length(X)-1)+1 && i<=7/8*(length(X)-1)+1
            S(i) = s2+h*(2.44406184*x^2-.22203097*x+.00723407);
            V(i) = -h/beta*(4.888124*x-.22203097);
            A(i) = h/beta^2*4.888124;
            J(i) = 0;
        elseif i>7/8*(length(X)-1)+1 && i<=length(X)
            S(i) = s2+h*(.38898448*x-.0309544*sin(4*pi*x)); 
            V(i) = -h/beta*.38898448*(1-cos(4*pi*x));
            A(i) = h/beta^2*4.888124*sin(4*pi*x);
            J(i) = -h/beta^3*61.425769*cos(4*pi*x);
        end
    end
end
