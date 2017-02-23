%% Test for feuler task a.
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 1000;
[t,y] = feuler(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');

%% 
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 1000;

i = 4;
x = zeros(1,i);

for k = 1:i
    [t,y] = feuler(u0,T,(n*k),alpha);
    x(k) = y(end,1);
end




%












%%