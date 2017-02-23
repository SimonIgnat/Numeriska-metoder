%% Test for feuler task a.
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 10*1000;
[t,y] = feuler(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');

%% 
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 1000;
x = [];
s
for k = 1:4
    
    [t,y] = feuler(u0,T,n*k,alpha);
    x = [x (y((T*(n*k)+1),1))]
end




%












%%