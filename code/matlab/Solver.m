 clear all 
 close all
%  задание шага по координатам
 hx = 0.005;
 hy = 0.01;
%  задание шага по времени
 tau = 0.0000007;
%  задание количества расчЄтных €чеек
 Nx = 200;
 Ny = 30;
%  задание констант и вычисление безразмерных к-тов
 nu = 0.1;
 Gamma = 1.4;
 T0 = 300;
 dT0 = 0.3;
 L0 = 6;
 u0 = 50;
 rho0 = 1.2;
 R = 8.3;
 mu = 0.000017;
 Cp = 2000;
 k = 0.025;

% объ€вление необходиных массивов


 u= zeros(Nx, Ny);
 v= zeros(Nx, Ny);
 rho = zeros(Nx, Ny);
 T = zeros(Nx, Ny);
 p = zeros(Nx, Ny);
 a = zeros(Nx, Ny);
 b= zeros(Nx, Ny);
 c= zeros(Nx, Ny);
 d= zeros(Nx, Ny);
M = 0.028;	
% задание начальных условий
		for j = 1:Ny           
			for i=1:Nx/2
                T(i,j) = 300;
				u(i,j) = 0;
				v(i,j) = 0;
				rho(i,j) = 6;
            end
 			for i=Nx/2:Nx
                T(i,j) = 300;
				u(i,j) = 0;
				v(i,j) = 0;
				rho(i,j) = 1.2;
            end           
        end
        p = rho.*R.*T./M;
x = 1:1:200;
figure(1)
h = line(x(1:200), rho(:,15));
set(h,'color', 'r');
set(h, 'erasemod', 'xor');
title('плотность');	

% цикл по времени
	for k=1:600
%       цикл по у - координате  
        for j=3:Ny-2
% 		уикл по х - координате
			for i=3:Nx-2
  
				 a(i,j) = rho(i,j) +  rhoCalcND(rho,u,v,T,p,i,j,tau,hx,hy);
				 b(i,j) = u(i,j) + uCalcND(rho,u,v,T,p,i,j,tau,hx,hy);
				 c(i,j) = v(i,j) + vCalcND(rho,u,v,T,p,i,j,tau,hx,hy);
				 d(i,j) = T(i,j) + TCalcND(rho,u,v,T,p,i,j,tau,hx,hy);
                 
           

            end

        end
    for j = 3 : Ny-2
        for i = 3 : Nx-2
 
         		rho(i,j) = a(i,j);
 				u(i,j) = b(i,j);
 				v(i,j) = c(i,j);
 				T(i,j) = d(i,j);
            
  
      
        end
            rho(1,j) = rho(3,j);
            rho(2,j) = rho(3,j);
            rho(Nx-1,j) = rho(Nx-2,j);
            rho(Nx,j) = rho(Nx-2,j);
            T(1,j) = T(3,j);
            T(2,j) = T(3,j);
            T(Nx-1,j) = T(Nx-2,j);
            T(Nx,j) = T(Nx-2,j);
    end
        

        k
        set(h,'XData',x(1:200),'YData',rho(:,15));
        drawnow;
        p = rho.*R.*T./M;
     end


%“ест —ода
% t = 0.00042
% p0 = R*5*1.2*300/0.028;
% pEnd = R*1.2*300/0.028;
% rhoEnd = 1.2;
% rho0 = 5*1.2;
% c1 = sqrt(Gamma*p0/rho0);
% u3 = 203;
% p3 = p0*(1 - (Gamma-1)*u3/(2*c1))^(2*Gamma/(Gamma-1));
% rho3 = rhoEnd*((Gamma - 1)*pEnd + (Gamma + 1)*p3)/((Gamma + 1)*pEnd + (Gamma - 1)*p3);
% rho2 = rho0*(p3/p0)^(1/Gamma);
% X0 = 0.5;
% u1 = 2/(Gamma+1)*(c1 - (X0 - x/200)/t);
% D = (p3 - pEnd)/(rhoEnd*u3);
% Xsw = X0 + D*t;
% Xrw1 = X0 - c1*t;
% Xrw2 = X0 - (c1 - (Gamma + 1)*u3/2)*t;
% X0t = X0 + u3*t;
% for i = 1:fix(Xrw1*200)+1
%     RHO(i) = 6;
%     U(i) = 0;
%     P(i) = p0;
% end
% for i = 1+fix(Xrw1*200) : fix(Xrw2*200)
%     RHO(i) = 6*(1 - (Gamma - 1)/2*u1(i)/c1)^(2/(Gamma-1));
%     U(i) = u1(i);
% 
% end
% for i = fix(Xrw2*200) : fix(X0t*200)
%     RHO(i) = rho2;
%     U(i) = u3;
%  
% end
% for i = fix(X0t*200) : fix(Xsw*200)
%     RHO(i) = rho3;
%     U(i) = u3;
%     P(i) = p3;
% end
% for i = fix(Xsw*200) : 200
%     RHO(i) = 1.2;
%     U(i) = 0;
%     P(i) = pEnd;
% end
% for i = 1+fix(Xrw1*200) : fix(Xrw2*200)+4
% 
%     P(i) = p0*(1 - (Gamma - 1)/2*u1(i)/c1)^(2/(Gamma-1));
% end
% for i = fix(Xrw2*200)+4 : fix(X0t*200)
% 
%     P(i) = p3;
% end
% figure(1)
% hold on
% 
% plot(x, RHO);
% 
% plot(x, rho(:,15));
% hold off
% 
% figure(2)
% hold on
% plot(x, U);
% 
% plot(x, u(:,15));
% hold off
% figure(3)
% hold on
% plot(x, P);
% 
% plot(x, p(:,15));
% hold off


	