% L = L0/(1 + I^2);
%  U = LdI/dt;
% I = -CdU/dt
% Använd RK-4
% Visa att följande uttryck kan härledas.
% d^2I/dt = 2I/(1+I^2)(dI/dt)^2 - I(1+I^2)/L0C
% t=0 ,I = 0 dI/dt = U0/L0.
% E(t) = CU(t)^2 + L0ln(1+I(t)^2) = konstant
L0 = 1;
C = 1*10^-6;
