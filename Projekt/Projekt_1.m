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
    xlim([0 0.015])
   
    subplot(3,2,plotCounter)
    plot(t,C*(L0./(1+y(:,1).^2).*y(:,2)).^2+L0*log(1+y(:,1).^2),'g')
    title(sprintf('E(t) with U0 = %.f',u0)); 
end


% ta fram frekvensen p� sinusen --> T f�r att ber�kna fouriertransform.

 
% a_k = 2/T int I(t) sin(kwt)dt'


%% Estimating max current and period at transient

T = 0.3*0.05;
n = 5000*0.05;

uVec = [240 1200 2400];
plotCounter = 0;

for u0 = uVec
    fprintf('U_0 = %.0f', u0)
    I0 = [0 u0/L0];
    [t,y] = RK4_projekt(I0,T,n, L0,C);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Imax
    maxIdx = find(y(:,1)==max(y(:,1)));
    % Fit second degree polynomial around maximum point
    poly = polyfit(t(maxIdx-1:maxIdx+1)', y(maxIdx-1:maxIdx+1,1),2);
    Imax1 = max(polyval(poly,t))
    
    % Increase step length to evaluate error
    poly = polyfit(t(maxIdx-2:maxIdx+2)', y(maxIdx-2:maxIdx+2,1),2);
    Imax2 = max(polyval(poly,t));
    
    Etrunk_imax = abs(Imax1-Imax2)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate T
    notFoundT = true
    
    
    i = 2; % To avoid first crossing of zero (sign(0) = 0)
    while notFoundT
        if(sign(y(i,1)) ~= sign(y(i+1,1)))
            period = t(i) + t(i+1);
            notFoundT = false;
        end
        i = i + 1;
    end
    period
    
    
    [t,y] = RK4_projekt(I0,T,n/2, L0,C);
    
    notFoundT = true
    
    i = 2; % To avoid first crossing of zero (sign(0) = 0)
    while notFoundT
        if(sign(y(i,1)) ~= sign(y(i+1,1)))
            period2h = t(i) + t(i+1);
            notFoundT = false;
        end
        i = i + 1;
    end
    
    EtrunkT = abs(period - period2h)

end


 




