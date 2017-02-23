function ut=fpendel(t,u)

g = 9.81;

L1=1.5;
L2=1;
m1=1;
m2=1.5;

s = sin(u(1)-u(3));
c = cos(u(1)-u(3));

f1 =(-m2*L2*u(4)^2*s-g*(m1+m2)*sin(u(1))-c*m2*(L1*u(2)^2*s-g*sin(u(3))))/...
    (L1*(m1+m2*s^2));
f2 =((m1+m2)*L1*u(2)^2*s-g*(m1+m2)*sin(u(3))+c*(m2*L2*u(4)^2*s+g*(m1+m2)*sin(u(1))))/...
    (L2*(m1+m2*s^2));

ut = [u(2); f1; u(4); f2];

