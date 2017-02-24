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
u0 = [0.5 0];
T = 10;
n = 1000;
[t,y] = rk4olin(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');




%%
alpha = 2;
alpha_nonlinear = 0.05;
T = 10;
n = 1000;

phis = [0.1, 1, 1.5, 3.1];

i = 1;
figure(1)
for phi = phis
    u0 = [phi,0];
    [t_lin, y_lin] = rk4(u0,T,n,alpha);
    [t_nlin, y_nlin] = rk4olin(u0,T,n,alpha);
    
    subplot(2,2,i)
    hold on
    plot(t_lin,y_lin(:,1), 'DisplayName', 'linear');
    plot(t_nlin, y_nlin(:,1), 'DisplayName', 'non-linear');
    
    legend('show')
    title(sprintf('Phi_0 = %.f',phi));
    hold off
    
    i = i +1;
end










