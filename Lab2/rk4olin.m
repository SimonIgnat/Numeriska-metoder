function [t,y] = rk4olin(u0,T,n,alpha)

%u0 = initial
% T = total time
% n = num of points??
% h = sample period (h = T/N),
%m*phi'' + a*phi' + mg/L * sin(phi) = 0;
%first order. x1 = phi, x2 = phi'.
% x1' = x2, x2' = (a*x2 + mg/L*x1)/M

x1 = zeros(1,n+1); % phi
x2 = zeros(1,n+1); % phi'

x1(1) = u0(1);
x2(1) = u0(2);

a = alpha; L = 1; m = 0.1; g = 9.81; h = T/n; 

for i= 1:n;
    k0 = h*x2(i);
    l0 = h*(-(a*x2(i)/m)-(g/L)*sin(x1(i)));
    
    k1 = h*(x2(i)+(1/2)*l0);
    l1 = h*(-a/m*(x2(i)+(1/2)*l0)-(g/L)*sin(x1(i)+(1/2)*k0));
    
    k2 = h*(x2(i)+(1/2)*l1);
    l2 = h*(-a/m*(x2(i)+(1/2)*l1)-(g/L)*sin(x1(i)+(1/2)*k1));
    
    k3 = h*(x2(i)+l2);
    l3 = h*(-a/m*(x2(i)+l2)-(g/L)*sin(x1(i)+k2));
    
    x1(i+1) = x1(i) + (1/6)*(k0+2*k1+2*k2+k3);
    x2(i+1) = x2(i) + (1/6)*(l0+2*l1+2*l2+l3);
    
    
end

t = 0:h:T;
y = [x1', x2'];


end

