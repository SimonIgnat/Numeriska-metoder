%% Task 1, finitia differensmetoden
L = 2; k = 3; t0 = 290; tL = 320; n = 4; h = L /(n+1); x = 0:h:L;

%Q(x) = -d/dx(k*dT/dx), T(0) = t0, T(L) = t1.
%Q(x) = 5000*exp(-500(x-0.25L)^2) + 1200exp(-10(x-0.7L)^2), 0\geq \geq L.
% centrala diffkvot: -T(i-1) + 2Ti - T(i+1) / h^2 = 1/k Q(x). i = 1,2,..,n.
%    2 -1  0  0 0 
%   -1 2 -1  0  0 
%    0 -1 2 -1  0     
%    0  0 -1 2 -1   *  h^2, i = 1,2,3,4.
%    0  0  0 -1 2 
%


Q = 5000*exp(-500*((x-0.25*L).^2))+ 1200*exp(-10*((x-0.7*L).^2));
A = diag(ones(1,n-1)*(-1),-1) + diag(ones(1,n)*2,0) + diag(ones(1,(n-1))*(-1),1);
A = A*(1/(h^2));
b = ([(Q(2)+(t0*k)) Q(3:(end-2)) (Q(end-1)+ (tL*k))]*(1/k))';

T = A\b+273.15 % just to get K.
plot(x(2:end-1),T)
%%

clear all, close all, clc

L = 2; k = 3; t0 = 290; tL = 320; 
N = [10, 20, 40, 80, 160, 320];

h = L/(N(1)+1);
x = 0:h:L;

pointsToEvaluate = 


figure

hold on
for n = N
   fprintf('n = %.1f\n', n)
   h = L /(n+1); 
   x = 0:h:L;
   
   
  
   Q = 5000*exp(-500*((x-0.25*L).^2))+ 1200*exp(-10*((x-0.7*L).^2));
   
   A = diag(ones(1,n-1)*(-1),-1) + diag(ones(1,n)*2,0) + diag(ones(1,(n-1))*(-1),1);
   A = A*(1/(h^2));
   b = ([(Q(2)+(t0*k)) Q(3:(end-2)) (Q(end-1)+ (tL*k))]*(1/k))';
    
   T = A\b + 273.15;
   T = [t0 T' tL];
    
   maxVal = max(T');
   minVal = min(T');
   meanVal = mean(T');
   
   fprintf('max = %.1f\n', maxVal)
   fprintf('min = %.1f\n', minVal)
   fprintf('mean = %.1f\n', meanVal)
     
   disp('_________________________')
   
   plot(x,T, 'DisplayName', sprintf('n = %.1f', n))
       
end

xlabel('m')
ylabel('K')

legend('show')
hold off

%%


