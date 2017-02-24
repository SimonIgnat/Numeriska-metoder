%% Test for feuler task a.
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 200000;
[t,y] = feuler(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');

%% 
alpha = 2;
u0 = [0.1 0];
T = 10;

N = [100, 500, 1000, 2000];
i = 1;
for n = N
    factorsOfn = [n, 2*n, 4*n];
    x = zeros(1,length(N));

    j = 1;
    for k = factorsOfn
        [t,y] = feuler(u0,T,(n*k),alpha);
        x(i) = y(end,1);
        j = j+1;
    end
    kvot(i) =  (x(3) - x(2))/(x(2) - x(1));    
end
%%
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 1000;
[t,y] = rk4(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');

%%
alpha = 0.05;
u0 = [1.5 0];
T = 10;
n = 1000;
[t,y] = rk4olin(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');




%%
