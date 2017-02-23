function [t,y] = feuler(u0,T,n,alpha)
%alluh akbar

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

a = alpha; L = 1; m = 0.1; g = 9.81; h = T/n; t = [];
%x' = A*x + Bu = 0.
t =[]; 
for i= 1:n;
    x1(i+1) = x1(i) + h*x2(i);
    x2(i+1) = x2(i) + h*(-(a*x2(i)/m) - (g/L)*x1(i));
    %t = [t (i-1)*h];
end

%t = [t n*h];
t = 0:h:T;
y = [x1', x2'];
end

