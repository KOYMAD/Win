function a = rhoCalcND(rho,u,v,T,p,i,j,tau,hx,hy)


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
a =  tau*(-rho(i,j)*((-u(i+2,j) + 8*u(i+1,j) - 8*u(i-1,j) + u(i-2,j))/(12*hx) + (-v(i,j+2) + 8*v(i,j+1) - 8*v(i,j-1) + v(i,j-2))/(12*hy)) - u(i,j)*((-rho(i+2,j) + 8*rho(i+1,j) - 8*rho(i-1,j) + rho(i-2,j))/(12*hx)) - v(i,j)*((-rho(i,j+2) + 8*rho(i,j+1) - 8*rho(i,j-1) + rho(i,j-2))/(12*hy)) + 150000*nu*((16*(rho(i+1,j) + rho(i-1,j)) - (rho(i+2,j) + rho(i-2,j)) - 30*rho(i,j))/12 + ((16*(rho(i,j+1) + rho(i,j-1)) - (rho(i,j+2) + rho(i,j-2)) - 30*rho(i,j))/12)));
end
