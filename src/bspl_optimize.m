function [S,V,A,J]= bspl_optimize(X,s,v,a,j,n)
% The B-spline

svaj_c = [s;v;a;j];
svaj_l = [size(s,1) size(v,1) size(a,1) size(j,1)];
BCs = sum(svaj_l);
u = BCs-n;
t_dot = ones(1,n+BCs);
t_dot(1:n) = t_dot(1:n)*X(1);
t_dot(BCs+1:end) = t_dot(BCs+1:end)*X(end);
Mc=repmat(1:BCs,BCs,1);
Mc(1,2:end) = 0;
Mc(size(s,1),1:end-1) = 0;
if ~isempty(v)
    if v(1,1) == t_dot(1)
        Mc(size(s,1)+1,3:end) = 0;
    end
    if v(end,1) == t_dot(end)
        Mc(size(s,1)+size(v,1),1:end-2) = 0; 
    end
end

if ~isempty(a)
    if a(1,1) == t_dot(1)
        Mc(size(s,1)+size(v,1)+1,4:end) = 0;
    end
    if a(end,1) == t_dot(end)
        Mc(size(s,1)+size(v,1)+size(a,1),1:end-3) = 0; 
    end
end
if ~isempty(j)
    if j(1,1) == t_dot(1)
        Mc(BCs-size(j,1)+1,5:end) = 0;
    end
    if j(end,1) == t_dot(end)
        Mc(BCs,1:end-4) = 0; 
    end
end

if u~=0
    t_bien = ones(u,2);
    for i = 1:u
        t_bien(i,1) = t_bien(i,1)*(X(1)+(X(end) - X(1))/180);
        t_bien(i,2) = t_bien(i,2)*(X(end)-(X(end) - X(1))/180);
    end
    tau = double((sqrt(5)-1)/2);
    t = t_bien(:,1)+(1-tau).*(t_bien(:,2)-t_bien(:,1));
    t = [t t_bien(:,1)+tau.*(t_bien(:,2)-t_bien(:,1))];
    while (abs(t_bien(1,1)-t_bien(1,2))) > (X(end) - X(1))/360
        %disp(t_bien)
        k = [];
        t_sort = dec2bin(0:(2^u-1)) - 47;
        for i = 1:size(t_sort,1)
           for j = 1:size(t_sort,2)
               if j~=size(t_sort,2) &&...
                       t(j,t_sort(i,j))>t(j+1,t_sort(i,j+1))
                   k=[k i]; break %#ok<AGROW>
               else, t_sort(i,j) = t(j,t_sort(i,j));
               end
            end
        end
        t_sort(k,:) = [];
        t_opti = zeros(size(t_sort,1),4);
        j = 0;
        for i = 1:size(t_sort,1)
            t_dot(n+1:n+u) = t_sort(i-j,:);
            [S,~,A] = K_bspl_curve(X,Mc,svaj_c,svaj_l,n,t_dot);
            if isempty(S)
                t_opti(i-j,:) = [];
                t_sort(i-j,:) = [];
                j=j+1;
            else
                t_opti(i-j,:) = [i-j abs(min(S)) max(S) max(A)-min(A)];
            end
        end
        t_opti = sortrows(t_opti,[2 4 3]);
        t_sort = t_sort(t_opti(1,1),:);
        % Loại bỏ và thay mới các giá trị chưa tối ưu
        for i = 1:u
            if(t(i,1) == t_sort(i))
                t_bien(i,2)=t(i,2);
                t(i,2) = t(i,1);
                t(i,1) = t_bien(i,1) + (1-tau)*(t_bien(i,2)-t_bien(i,1));
            else
                t_bien(i,1) = t(i,1);
                t(i,1) = t(i,2);
                t(i,2) = t_bien(i,1) + tau*(t_bien(i,2)-t_bien(i,1));
            end
        end
    end
    t_dot(n+1:n+u) = t_sort;
end
[S,V,A,J] = K_bspl_curve(X,Mc,svaj_c,svaj_l,n,t_dot);
    
function [fx,der1,der2,der3] = K_bspl_curve(X,Mc,con,con_index,n,t)
% The B-spline
BCs = sum(con_index);
for i=1:con_index(1)
    for k = Mc(i,:)
        if k == 0, continue; end
        Mc(i,k) = K_bspline_basis(k-1,n,t,con(i,1));
    end
end

% BCs is Velocity
for i=i+1:i+ con_index(2)
    for k = Mc(i,:)
        if k == 0 
            continue
        elseif k==1
            Mc(i,k) = -(n-1)...
                /(t(1+n)-t(2))*K_bspline_basis(1,n-1,t,con(i,1));
        elseif k==BCs
            Mc(i,k) = (n-1)/(t(BCs+n-1)-t(BCs))...
                *(K_bspline_basis(BCs-1,n-1,t,con(i,1)));
        else
            Mc(i,k) = (n-1)...
                *(K_bspline_basis(k-1,n-1,t,con(i,1))/(t(k+n-1)-t(k))...
                -K_bspline_basis(k,n-1,t,con(i,1))/(t(k+n)-t(k+1)));
        end
    end
end

% BCs is Accelarate
for i = i+1 : i+con_index(3)
    for k = Mc(i,:)
        if k == 0 
            continue
        elseif k == 1
            Mc(i,1) = (n-1)*(n-2)...
               /((t(1+n)-t(2))*(t(k+n)-t(3)))...
               *K_bspline_basis(2,n-2,t,con(i,1));
        elseif k == 2
             Mc(i,2) = (n-1)*(n-2)/((t(2+n)-t(4))*(t(2+n)-t(2+1)))...
                *K_bspline_basis(3,n-2,t,con(i,1))...
                -(n-1)*(n-2)*( t(n+1)+t(2+n)-t(2)-t(3) )...
                /((t(n+1)-t(3))*(t(n+1)-t(2))*(t(2+n)-t(3)))...
                *K_bspline_basis(2,n-2,t,con(i,1));
        elseif k == BCs-1
             Mc(i,k) = - (n-1)*(n-2)*(t(k+n-1)+t(k+n)-t(k)-t(k+1))...
                /((t(k+n-1)-t(k+1))*(t(k+n-1)-t(k))*(t(k+n)-t(k+1)))...
                *K_bspline_basis(k,n-2,t,con(i,1))...
                +(n-1)*(n-2)/((t(k+n-1)-t(k))*(t(k+n-2)-t(k)))...
                *K_bspline_basis(k-1,n-2,t,con(i,1));
        elseif k == BCs
            Mc(i,k) = (n-1)*(n-2)/((t(k+n-1)-t(k))*(t(k+n-2)-t(k)))...
               *K_bspline_basis(k-1,n-2,t,con(i,1));
        else
            Mc(i,k) = (n-1)*(n-2)/( (t(k+n)-t(k+2))*(t(k+n)-t(k+1)))...
               *K_bspline_basis(k+1,n-2,t,con(i,1))...
               -(n-1)*(n-2)*( t(k+n-1)+t(k+n)-t(k)-t(k+1))...
               /((t(k+n-1)-t(k+1))*(t(k+n-1)-t(k))*(t(k+n)-t(k+1)))...
               *K_bspline_basis(k,n-2,t,con(i,1))...
               +(n-1)*(n-2)/((t(k+n-1)-t(k))*(t(k+n-2)-t(k)))...
               *K_bspline_basis(k-1,n-2,t,con(i,1));
        end
    end
end
% BCs is Jerk

for i = i+1 : i+con_index(4)
    for k = Mc(i,:)
        if k == 0 
            continue
        elseif k == 1
            Mc(i,1) = (n-1)*(n-2)*(n-3)*...
               (1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))*(t(k+3) - t(n+k))))*...
               K_bspline_basis(3,n-3,t,con(i,1));
        elseif k == 2
             Mc(i,2) = (n-1)*(n-2)*(n-3)*...
                 (-(1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))) + (1/(t(k) - t(n+k-1)) +...
                 1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+2) - t(n+k-1))) *...
                 K_bspline_basis(3,n-3,t,con(i,1)) +...
                 (n-1)*(n-2)*(n-3)*...
                 (1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))*(t(k+3) - t(n+k))))*...
                 K_bspline_basis(4,n-3,t,con(i,1));
        elseif k == 3
             Mc(i,k) = (n-1)*(n-2)*(n-3)*...
                 ((1/((t(k) - t(n+k-2))*(t(k) - t(n+k-1))) + (1/(t(k) - t(n+k-1)) +...
                 1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+1) - t(n+k-2))) *...
                 K_bspline_basis(3,n-3,t,con(i,1)) +...
                 (n-1)*(n-2)*(n-3)*...
                 (-(1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))) + (1/(t(k) - t(n+k-1)) +...
                 1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+2) - t(n+k-1))) *...
                 K_bspline_basis(4,n-3,t,con(i,1)) +...
                 (n-1)*(n-2)*(n-3)*...
                 (1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))*(t(k+3) - t(n+k))))*...
                 K_bspline_basis(5,n-3,t,con(i,1));
        elseif k == BCs
            Mc(i,k) = (n-1)*(n-2)*(n-3) *...
                (-1/((t(k) - t(n+k-3))*(t(k) - t(n+k-2))*(t(k) - t(n+k-1)))) *...
                K_bspline_basis(k-1,n-3,t,con(i,1));
        elseif k == BCs - 1
            Mc(i,k) = (n-1)*(n-2)*(n-3) *...
                ((1/((t(k) - t(n+k-2))*(t(k) - t(n+k-1))) + (1/(t(k) - t(n+k-1)) +...
                1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+1) - t(n+k-2)))*...
                K_bspline_basis(k,n-3,t,con(i,1)) +...
                (n-1)*(n-2)*(n-3) *...
                (-1/((t(k) - t(n+k-3))*(t(k) - t(n+k-2))*(t(k) - t(n+k-1)))) *...
                K_bspline_basis(k-1,n-3,t,con(i,1));
        elseif k == BCs - 2
            Mc(i,k) = (n-1)*(n-2)*(n-3) *...
                (-(1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))) + (1/(t(k) - t(n+k-1)) +...
                1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+2) - t(n+k-1)))*...
                K_bspline_basis(k+1,n-3,t,con(i,1)) +...
                (n-1)*(n-2)*(n-3) *...
                ((1/((t(k) - t(n+k-2))*(t(k) - t(n+k-1))) + (1/(t(k) - t(n+k-1)) +...
                1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+1) - t(n+k-2)))*...
                K_bspline_basis(k,n-3,t,con(i,1)) +...
                (n-1)*(n-2)*(n-3) *...
                (-1/((t(k) - t(n+k-3))*(t(k) - t(n+k-2))*(t(k) - t(n+k-1)))) *...
                K_bspline_basis(k-1,n-3,t,con(i,1));
        else
            Mc(i,k) = (n-1)*(n-2)*(n-3) *...
                (-1/((t(k) - t(n+k-3))*(t(k) - t(n+k-2))*(t(k) - t(n+k-1)))) *...
                K_bspline_basis(k-1,n-3,t,con(i,1)) + ...
                (n-1)*(n-2)*(n-3) *...
                ((1/((t(k) - t(n+k-2))*(t(k) - t(n+k-1))) + (1/(t(k) - t(n+k-1)) +...
                1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+1) - t(n+k-2))) *...
                K_bspline_basis(k,n-3,t,con(i,1)) + ...
                (n-1)*(n-2)*(n-3) *...
                (-(1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))) + (1/(t(k) - t(n+k-1)) +...
                1/(t(k+1) - t(n+k)))/(t(k+1) - t(n+k-1)))/(t(k+2) - t(n+k-1))) *...
                K_bspline_basis(k+1,n-3,t,con(i,1)) +...
                (n-1)*(n-2)*(n-3) *...
                (1/((t(k+1) - t(n+k))*(t(k+2) - t(n+k))*(t(k+3) - t(n+k))))*...
                K_bspline_basis(k+2,n-3,t,con(i,1));
        end
    end
end
%disp(Mc)
if any(isnan(Mc),'all')
    fx = [];
    der1 = [];
    der2 = []; 
    der3 = [];
else
    Q = Mc^-1*con(:,2);
    %Q = Q';
    fx = K_bspl(Q,n,t,X,0);
    der1 = K_bspl(Q,n,t,X,1);
    der2 = K_bspl(Q,n,t,X,2);
    der3 = K_bspl(Q,n,t,X,3);
end


function y = K_bspl(c,n,t,x,devir)
% Derivation of B-spline
y=zeros(size(x));
for l = 1:devir
    for j = 1:length(c)-1
            c(j) = (n-l)*(c(j+1)-c(j))/(t(j+n)-t(j+l));
    end
    c(end) = [];
end
for i = 1:length(c)
    y = y + c(i)*K_bspline_basis(i-1+devir,n-devir,t,x);
end

function y = K_bspline_basis(j,n,t,x)

y = zeros(size(x));
if n > 1
    b = K_bspline_basis(j,n-1,t,x);
    dn = x - t(j+1);
    dd = t(j+n) - t(j+1);
    if dd ~= 0  % indeterminate forms 0/0 are deemed to be zero
        y = y + b.*(dn./dd);
    end
    b = K_bspline_basis(j+1,n-1,t,x);
    dn = t(j+n+1) - x;
    dd = t(j+n+1) - t(j+1+1);
    if dd ~= 0
        y = y + b.*(dn./dd);
    end
elseif t(j+2) < t(end)  % treat last element of knot vector as a special case
	y(t(j+1) <= x & x < t(j+2)) = 1;
else
	y(t(j+1) <= x) = 1;
end