function [S,V,A,J] = sinu_mod(X,s1,s2)
% MODIFIED SINUSOIDAL ACCELERATION
 
if nargin < 3,error('Ham co it nhat 3 doi so.');end

S = zeros(1,length(X)); V = S; A = S; J = S;
h=abs(s2-s1);
beta = X(length(X))-X(1);

if s2-s1>0
    for i = 1 : length(X)
        x=(X(i)-X(1))/beta;
        if i>=1 && i<= (length(X)-1)/8+1
            S(i) = s1+h*(.43990085*x - .0350062*sin(4*pi*x)); 
            V(i) = h/beta*.43990085*(1 - cos(4*pi*x)); 
            A(i) = h/beta^2*5.5279571*sin(4*pi*x);
            J(i) = h/beta^3*69.4663577*cos(4*pi*x);
        elseif i>(length(X)-1)/8+1 && i<=7/8*(length(X)-1)+1
            S(i) = s1+h*(.28004957+.43990085*x-.31505577*cos(4/3*pi*x-pi/6));
            V(i) = h/beta*.43990085*(1+3*sin(4/3*pi*x-pi/6));
            A(i) = h/beta^2*5.5279571*cos(4/3*pi*x-pi/6);
            J(i) = h/beta^3*-23.1553*sin(4/3*pi*x-pi/6);
        elseif i>7/8*(length(X)-1)+1 && i<=length(X)
            S(i) = s1+h*(.56009915+.43990085*x-.0350062*sin(4*pi*x-2*pi));
            V(i) = h/beta*.43990085*(1-cos(4*pi*x-2*pi));
            A(i) = h/beta^2*5.5279571*sin(4*pi*x-2*pi);
            J(i) = h/beta^3*69.4663577*cos(4*pi*x-2*pi);
        end
    end
else 
    for i = 1 : length(X)
        x=1-(X(i)-X(1))/beta;
        if i>=1 && i<= (length(X)-1)/8+1
            S(i) = s2+h*(.56009915 + .43990085*x-.0350062*sin(4*pi*x-2*pi));
            V(i) = -h/beta*.43990085*(1 - cos(4*pi*x-2*pi));
            A(i) = h/beta^2*5.5279571*sin(4*pi*x-2*pi);
            J(i) = -h/beta^3*69.4663577*cos(4*pi*x-2*pi);
        elseif i>(length(X)-1)/8+1 && i<=7/8*(length(X)-1)+1
            S(i) = s2+h*(.28004957+.43990085*x-.31505577*cos(4/3*pi*x-pi/6));
            V(i) = -h/beta*.43990085*(1 + 3*sin(4/3*pi*x - pi/6));
            A(i) = h/beta^2*5.5279571*cos(4/3*pi*x - pi/6);
            J(i) = -h/beta^3*-23.1553*sin(4/3*pi*x - pi/6);
        elseif i>7/8*(length(X)-1)+1 && i<=length(X)
            S(i) = s2+h*(.43990085*x-.0350062*sin(4*pi*x)); 
            V(i) = -h/beta*.43990085*(1-cos(4*pi*x)); 
            A(i) = h/beta^2*5.5279571*sin(4*pi*x);
            J(i) = -h/beta^3*69.4663577*cos(4*pi*x);
        end
    end
end