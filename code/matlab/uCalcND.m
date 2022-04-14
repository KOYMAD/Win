function a = uCalcND(rho,u,v,T,p,i,j,tau,hx,hy)

nu = 0.1;
Gamma = 1.4;
T0 = 300;
dT0 = 0.3;
L0 = 6;
u0 = 250;
rho0 = 1.2;
R = 8.3;
mu = 0.000017;
Cp = 1000;
k = 0.025;
Re = u0*L0*rho0/mu;
Pr = mu*Cp/k;
M = u0/(Gamma*R*T0)^0.5;
a = tau*(-1*u(i,j)*((u(i+1,j) - u(i-1,j) )/(2*hx)) - v(i,j)*(u(i,j+1) - u(i,j-1))/(2*hy)) -1/rho(i,j)*((p(i+1,j) - p(i-1,j))/(2*hx)) + mu/(rho(i,j))*(4/3*(16*(u(i+1,j) + u(i-1,j)) - (u(i+2,j) + u(i-2,j)) - 30*u(i,j))/(12*hx^2) + (16*(u(i,j+1) + u(i,j-1)) - (u(i,j+2) + u(i,j-2)) - 30*u(i,j))/(12*hy^2) + 1/3*((v(i+1,j+1) - v(i+1,j-1) - 2*(v(i,j+1) - v(i,j-1)) + v(i-1,j+1) -v(i-1,j-1))/(4*hx*hy))) + 350000*nu*((16*(u(i+1,j) + u(i-1,j)) - (u(i+2,j) + u(i-2,j)) - 30*u(i,j))/(12) + (16*(u(i,j+1) + u(i,j-1)) - (u(i,j+2) + u(i,j-2)) - 30*u(i,j))/(12)));
end