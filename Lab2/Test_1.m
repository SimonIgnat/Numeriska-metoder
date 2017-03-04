%% Test for feuler task a.
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 1000;
[t,y] = feuler(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');
y(end,1)


%% 
alpha = 2;
u0 = [0.1 0];
T = 10;

N = [10, 20,40,80, 160, 320, 640, 1280, 2560];
fprintf('Framm?t euler nogranhetstabell\n')
fprintf('__________________________________________________________________\n')
fprintf('h            u_h           u_h-u_(h/2)    kvot          ordning\n')
for n = N
    
    factorsOfn = [n, 2*n, 4*n];
    x = zeros(1,length(N));

    j = 1;
    for k = factorsOfn
        [t,y] = feuler(u0,T,(n*k),alpha);
        x(j) = y(end,1);
        j = j + 1;
    end
   
    kvot = (x(1) - x(2))/(x(2) - x(3));
    fprintf('%f     %f      %f       %f      %f\n', T/n,x(1), x(1)-x(2), kvot, log2(kvot));
end 
%%
alpha = 2;
u0 = [0.1 0];
T = 10;
n = 1000;
[t,y] = rk4(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');
y(end,1)

%%
alpha = 2;
u0 = [0.1 0];
T = 10;

N = [10, 20,40,80, 160, 320, 640, 1280, 2560];
fprintf('Runge-Kutta4 nogranhetstabell\n')
fprintf('__________________________________________________________________\n')
fprintf('h            u_h           u_h-u_(h/2)    kvot          ordning\n')
for n = N
    
    factorsOfn = [n, 2*n, 4*n];
    x = zeros(1,length(factorsOfn));

    j = 1;
    for k = factorsOfn
        [t,y] = rk4(u0,T,(n*k),alpha);
        x(j) = y(end,1);
        j = j + 1;
    end
   
    kvot = (x(1) - x(2))/(x(2) - x(3));
    fprintf('%f     %f      %f       %f      %f\n', T/n,x(1), x(1)-x(2), kvot, log2(kvot));
end 


%%
alpha = 0.05;
u0 = [0.5 0];
T = 10;
n = 1000;
[t,y] = rk4olin(u0,T,n,alpha);
plot(t,y(:,1),'b',t,y(:,2),'r');


%%

clear all, close all
alpha = 2;
alpha_nonlinear = 0.05;
T = 10;
n = 1000;

phis = [0.01, 1, 1.5, 3.1];

i = 1;

figure(1)
for phi = phis
    u0 = [phi,0];
    [t_lin, y_lin] = rk4(u0,T,n,alpha_nonlinear);
    [t_nlin, y_nlin] = rk4olin(u0,T,n,alpha_nonlinear);
    
  
    subplot(2,2,i)
    hold on
    plot(t_lin,y_lin(:,1), 'DisplayName', 'linear');
    plot(t_nlin, y_nlin(:,1), 'DisplayName', 'non-linear');
    
    legend('show')
    title(sprintf('Phi_0 = %.f',phi));
    hold off
    
    i = i +1;
end

%% d)

clear all, close all

alpha_nonlinear = 0.05;
T = 10;
n = 100;
L = 1;

phis = [0.1]

i = 1;

figure(1)
for phi = phis
    u0 = [phi,0];
    [t_nlin, y_nlin] = rk4olin(u0,T,n,alpha_nonlinear);
    P = y_nlin(:,1);
    
    subplot(2,2,i)
    for k=1:length(P) 
        plot([-L L],[0,0],[0 L*sin(P(k))],[0 -L*cos(P(k))],'-o'); 
        axis equal 
        axis(1.2*[-L L -L L]); 
        drawnow 
    end 
    title(sprintf('Phi_0 = %.f',phi));
    i = i + 1;
end


%%



















