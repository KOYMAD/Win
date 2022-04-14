#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>

double eos(double rho, double T, double R, double M)
{
	double p = rho * R * T / M;
	return (p);
}

int main(void)
{
    // TODO: у всех величин должны быть указаны размерности  --- Done
    //Задание шагов по координатам и времени, соответственно
	const double hx = 0.005 /* м*/;
    const double hy = 0.01 /* м*/;
    const double tau = 0.0000007 /* с*/;
    // Задание количества ячеек по координатам
    const int Nx = 200;
    const int Ny = 30;
    // Объявление массивов переменных, задание постоянных
    double nuRho = 15000.0 /* 1/с */;
    double nuV = 35000.0 /* 1/с */;
    double nuU = 20000.0 /* 1/с */;
    double nuT = 1.0 /* 1/с */;
    double Gamma = 1.4;
    double T0 = 300.0 /* К */;
    double dT0 = 0.3 /* К */;
    double L0 = 6.0 /* м */;
    double u0 = 50.0 /* м/с */;
    double rho0 = 1.2 /* кг/м^3 */;
    double R = 8.3 /* Дж * К^-1 * Моль^-1 */;
    double mu = 0.000017 /* Па * с */;
    double Cp = 1000 /* Дж * кг^-1 * К^-1 */;
    double k = 0.025 /* Вт * м^-1 * К^-1 */;
    double u[Nx][Ny];
    double v[Nx][Ny];
    double rho[Nx][Ny];
    double T[Nx][Ny];
    double p[Nx][Ny];
    double rho_new[Nx][Ny];
    double u_new[Nx][Ny];
    double v_new[Nx][Ny];
    double T_new[Nx][Ny];
    double M = 0.028 /* кг / Моль */;
    double Cv = Cp / Gamma;
    FILE *result;
    result = fopen("result.txt", "w");
    // Задание начальных условий - х-разрыв плотности
    for (int j = 0; j < Ny; j++)
    {
        for (int i = 0;i < Nx;i++)
        {
            if ( i < Nx / 2 )
            {
	        T[i][j] = 300.0;
	        u[i][j] = 0.0;
		v[i][j] = 0.0;
		rho[i][j] = 6.0;
                // TODO: сделать отдельную функцию eos, которая возвращает давление по температуре и плостности --- Done
		p[i][j] = eos(rho[i][j], T[i][j], R, M);
            }
            else
            {
	        T[i][j] = 300.0;
		u[i][j] = 0.0;
		v[i][j] = 0.0;
	    rho[i][j] = 1.2;
		p[i][j] = eos(rho[i][j], T[i][j], R, M);
	    }
        }
    }
    // Цикл по времени
    // TODO: очень неудачное обозначение для номера шага по времени o, легко спутать с нулем --- Done
    for ( int TimeStep = 0; TimeStep < 600; TimeStep++ )
    {
        // Цикл по y - координате
        // TODO: совсем отсутсвтует понятие граничного условия и его обработка --- processing
	for ( int j = 2; j < Ny - 2; j++ )
	{
	    // Цикл по x - координате
	    for (int i = 2;i < Nx - 2;i++)
	    {
                // Формула 1) из (16) из \trunk\chuprov\docs\ns_algorithm.docx
                // TODO: это не формула (1), это формула без номера на стр. 5, аналогично остальное ---Done
                // TODO: обозначения a, b, c ... - неудачные, без смысла. Логично называть их rho_new или rho_next и так далее ---Done
                rho_new[i][j] = rho[i][j] + tau * ( -rho[i][j] * ((u[i + 1][j] - u[i - 1][j]) / (2.0 * hx)
                                                          + (v[i][j + 1] - v[i][j - 1]) / (2.0 * hy) )
                                             - u[i][j] * ( (rho[i + 1][j] - rho[i - 1][j]) / (2.0 * hx) )
                                             - v[i][j] * ( (rho[i][j + 1] - rho[i][j - 1]) / (2.0 * hy) )
                                             + nuRho * ( (rho[i + 1][j] - 2.0 * rho[i][j] + rho[i - 1][j])
                                                      +  (rho[i][j + 1] - 2.0 * rho[i][j] + rho[i][j - 1]) ) );
				// Формула 2) из (16) из \trunk\chuprov\docs\ns_algorithm.docx
                // TODO: в документе нет подобного разностного уравнения ---Done
                u_new[i][j] = u[i][j] + tau * ( - u[i][j] * ((u[i + 1][j] - u[i - 1][j]) / (2.0 * hx))
                                            - v[i][j] * ((u[i][j + 1] - u[i][j - 1]) / (2.0 * hy))
                                            - 1 / rho[i][j] * ((p[i + 1][j] - p[i - 1][j]) / (2 * hx)) 
                + mu / (rho[i][j])*(4 / 3 * ((u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) / pow(hx, 2)) + ((u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / pow(hy, 2)) + 1 / 3 * ((2 * v[i + 1][j + 1] - 2 * v[i - 1][j + 1] 
                - v[i + 1][j] + v[i - 1][j]) / (8 * hx*hy))) + nuV * ((u[i + 1][j] - 2 * u[i][j] + u[i - 1][j]) + (u[i][j + 1] - 2 * u[i][j] + u[i][j - 1])));
				// Формула 3) из (16) из \trunk\chuprov\docs\ns_algorithm.docx
                    v_new[i][j] = v[i][j] + tau * (-1 * u[i][j] * ((v[i + 1][j] - v[i - 1][j]) / (2 * hx)) - v[i][j] * ((v[i][j + 1] - v[i][j - 1]) / (2 * hy)) - 1 / rho[i][j] * ((p[i][j + 1] - p[i][j - 1]) / (2 * hy)) 
                    + mu / (rho[i][j])*(4 / 3 * ((v[i][j + 1] - 2 * v[i][j] + v[i][j + 1]) / pow(hy, 2)) + ((v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) / pow(hx, 2)) + 1 / 3 * ((2 * u[i + 1][j + 1] - 2 * u[i - 1][j + 1] 
                    - u[i + 1][j] + u[i - 1][j]) / (8 * hx*hy))) + nuU * ((v[i + 1][j] - 2 * v[i][j] + v[i - 1][j]) + (v[i][j + 1] - 2 * v[i][j] + v[i][j - 1])));
					// Формула 4) из (16) из \trunk\chuprov\docs\ns_algorithm.docx
                    T_new[i][j] = T[i][j] + tau * (-1 * u[i][j] * ((T[i + 1][j] - T[i - 1][j]) / (2 * hx)) - v[i][j] * ((T[i][j + 1] - T[i][j - 1]) / (2 * hy)) - p[i][j]/(Cv*rho[i][j])*((u[i + 1][j] - u[i - 1][j]) / (2 * hx) 
                    + (v[i][j + 1] - v[i][j - 1]) / (2 * hy)) + k/(Cv*rho[i][j])*((T[i + 1][j] - 2 * T[i][j] + T[i - 1][j])/(pow(hx,2)) + (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1])/(pow(hy, 2))) 
                    + nuT * ((T[i + 1][j] - 2 * T[i][j] + T[i - 1][j]) + (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1])));
		}
	    }
	    for (int j = 2;j < Ny - 2;j++)
	    {			
	        for (int i = 2;i < Nx - 2;i++)
			{
				rho[i][j] = rho_new[i][j];
				u[i][j] = u_new[i][j];
				v[i][j] = v_new[i][j];
				T[i][j] = T_new[i][j];
				p[i][j] = rho[i][j] * R*T[i][j] / M;
				// граничные условия по оси Х
				u[i][1] = u[i][2];
				rho[i][1] = rho[i][2];
				v[i][1] = v[i][2];
				T[i][1] = T[i][2];

				u[i][Ny - 2] = u[i][Ny - 3];
				rho[i][Ny - 2] = rho[i][Ny - 3];
				v[i][Ny - 2] = v[i][Ny - 3];
				T[i][Ny - 2] = T[i][Ny - 3];

			}
			//  граничные условия по оси У
			u[1][j] = u[2][j];
			rho[1][j] = rho[2][j];
			v[1][j] = v[2][j];
			T[1][j] = T[2][j];

			u[Nx - 2][j] = u[Nx -3][j];
			rho[Nx - 2][j] = rho[Nx - 3][j];
			v[Nx - 2][j] = v[Nx - 3][j];
			T[Nx - 2][j] = T[Nx - 3][j];

	    }
		
	    printf("%d \n" , TimeStep);
	}
	// Запись результатов
	for (int i = 2;i < Nx - 1;i++)
	{
	    //fprintf(result, " %d", i);
	    fprintf(result, " %5.5e \n", rho[i][15]);
	}
	return 0;
}

