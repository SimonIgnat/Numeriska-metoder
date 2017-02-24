function [t,y] = rk4olin(u0,T,n,alpha)
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


% x1' = x2;
% x2' = -(a*x2(i)/m) - (g/L)*x1(i);
%dydx = f(x).


for i= 1:n;
    k11 = h*x2(i);
    k21 = h*(x2(i)+(1/2)*k11);              % due to linearity, f(x+h) = f(x)+f(h).
    k31 = h*(x2(i)+(1/2)*k21);
    k41 = h*(x2(i)+k31);
    x1(i+1) = x1(i) + (k11+2*k21+2*k31+k41)/6;
    
    k12 = h*(-(a*x2(i)/m)-(g/L)*sin(x1(i)*pi/180));
    k22 = h*(x2(i) + (1/2)*k12);
    k32 = h*(x2(i) + (1/2)*k22);
    k42 = h*(x2(i) + k32);
    x2(i+1) = x2(i) +(k12+2*k22+2*k32+k42)/6;
    %t = [t (i-1)*h];
end

%t = [t n*h];
t = 0:h:T;
y = [x1', x2'];


end

