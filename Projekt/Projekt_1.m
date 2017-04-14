% L = L0/(1 + I^2);
%  U = LdI/dt;
% I = -CdU/dt
% Anv�nd RK-4
% Visa att f�ljande uttryck kan h�rledas.
% d^2I/dt = 2I/(1+I^2)(dI/dt)^2 - I(1+I^2)/L0C
% t=0 ,I = 0 dI/dt = U0/L0.
% E(t) = CU(t)^2 + L0ln(1+I(t)^2) = konstant



%% RK4
clear all, close all
L0 = 1;
C = 1*10^-6;

T = 0.3*2;
n = 5000*2;

uVec = [240 1200 2400];
plotCounter = 0;
for u0 = uVec
    plotCounter = plotCounter + 1;
    subplot(3,2,plotCounter)
    I0 = [0 u0/L0]
    [t,y] = RK4_projekt(I0,T,n, L0,C);
    plot(t,y(:,1),'b');
    title(sprintf('I(t) with U0 = %.f',u0));
    plotCounter = plotCounter + 1;
    xlim([0 0.01])
   
    subplot(3,2,plotCounter)
    plot(t,C*(L0./(1+y(:,1).^2).*y(:,2)).^2+L0*log(1+y(:,1).^2),'g')
    title(sprintf('E(t) with U0 = %.f',u0)); 
end


% ta fram frekvensen p� sinusen --> T f�r att ber�kna fouriertransform.


% a_k = 2/T int I(t) sin(kwt)dt'

%%
k = 14;
a = zeros(1,k);

for k = 1:14
    a(k) = 2*(y(,1)*sin(k*2*pi/T))
end





