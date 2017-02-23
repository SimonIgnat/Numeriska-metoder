%% Test for feuler task a.
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 100000;
[t,y] = feuler(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');

%% 

format long
alpha = 2;
u0 = [0.1 0];
T = 10;

N = [1000, 2000];

kvot = zeros(1, length(N));
i = 1;
for n = N
    factorsOfn = [n, 2*n, 4*n];
    x = zeros(1,length(factorsOfn));

    j = 1;
    for k = factorsOfn
        [t,y] = feuler(u0,T,(n*k),alpha);
        x(j) = y(end,1);
        j = j+1;
    end
    kvot(i) =  (x(3) - x(2))/(x(2) - x(1));   
    i = i + 1;
end

%%



for i = 1:length(N)
    disp(sprintf('N = %f, kvot = %f', N(i), kvot(i)));
end




%












%%