function a = TCalcND(rho,u,v,T,p,i,j,tau,hx,hy)

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
Cv = Cp/Gamma;
a = tau*(-1*u(i,j)*((-T(i+2,j) + 8*T(i+1,j) - 8*T(i-1,j) + T(i-2,j))/(12*hx)) - v(i,j)*((-T(i,j+2) + 8*T(i,j+1) - 8*T(i,j-1) + T(i,j-2))/(12*hy)) - p(i,j)/(rho(i,j)*Cv)*((-u(i+2,j) + 8*u(i+1,j) - 8*u(i-1,j) + u(i-2,j))/(12*hx) + (-v(i,j+2) + 8*v(i,j+1) - 8*v(i,j-1) + v(i,j-2))/(12*hy)) + k/(Cv*rho(i,j))*((16*(T(i+1,j) + T(i-1,j)) - (T(i+2,j) + T(i-2,j)) - 30*T(i,j))/(12*hx^2) + (16*(T(i,j+1) + T(i,j-1)) - (T(i,j+2) + T(i,j-2)) - 30*T(i,j))/(12*hy^2)) + 10*nu*((T(i+1,j) - 2*T(i,j) + T(i-1,j)) + (T(i,j+1) - 2*T(i,j) + T(i,j-1))));
end
