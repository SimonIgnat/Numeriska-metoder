u0 = [2 0 0.8 0];

T= 5;
n = 1000;
h = T/n
phi1 = zeros(1,n);
phi2 = zeros(1,n);

[T_out, Y_out] = ode45(@fpendel,[0 50],u0);


%plot(T_out,Y_out(:,1),'r',T_out,Y_out(:,3));

L1 = 1.5;
L2 = 1;
M1 = 1;
M2 = 1.5;

for k=1:length(Y_out)
plot([-L1 L1],[0,0],[0 L1*sin(Y_out(k,1))],[0 -L1*cos(Y_out(k,1))],'-o');hold on;
plot([L1*sin(Y_out(k,1)) (L1*sin(Y_out(k,1))+L2*sin(Y_out(k,3)))],[-L1*cos(Y_out(k,1)) -L1*cos(Y_out(k,1))-L2*cos(Y_out(k,3))],'-o');
hold off;
axis(1.2*[-L1-L2 L1+L2 -L1-L2 L1+L2]);
drawnow
end