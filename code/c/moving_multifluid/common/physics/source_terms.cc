// source_terms.cc
// Расчет силы межфазного трения
// (c) Уткин Павел, 2014, Порошина Ярослава, 2018
// Создан: 7 января 2014 г.

#include "source_terms.h"

// Расчет силы межфазного трения
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает силу трения
double calc_friction_force( struct ParametersCommon *paramsc, double v_ncons[M]) {

    // для реализации формулы HOUIM
    double Ksg; // коэффициент пропорциональности между силой трения и v_rel
    double alpha_g = 1.0 - v_ncons[B_DISP];
    double alpha_s = v_ncons[B_DISP];
    double rho_g = v_ncons[R_GAS];

    if ( !paramsc->is_friction_force )
        // силу трения учитывать не просят
        return 0.0;
    
    double v_rel = v_ncons[V_DISP] - v_ncons[V_GAS]; // относительная скорости движения фаз

    if ( fabs( v_rel ) < paramsc->eps_general )
        // равна нулю относительная скорость фаз, сила трения ноль
        return 0.0;

    // содержательные формулы для силы трения
    switch ( paramsc->fric_force_formula_num ) {
        
        case ROGUE:
            // формула для засыпки сферических частиц из
            // Rogue et al. Experimental and numerical investigation of the shock-induced fluidization of a particles bed // Shock Waves. - 1998. - V. 8. - P. 29 - 45.
            return ( 0.75 * 0.6 /* постоянное, допускающее варьирование значение Cd */ * rho_g * alpha_s / paramsc->particle_diameter * fabs( v_rel ) * v_rel );
        
        case HOUIM:
            // формула для засыпки сферических частиц из
            // Houim R.W., Oran E.S. A multiphase model for compressible granular-gaseous flows: formulation and initial tests // J. Fluid Mech. - 2016. - V. 789. - P. 166 - 220.
            if ( alpha_g >= 0.8 ) {
                Ksg = 0.75 * calc_Cd( paramsc, v_ncons ) * rho_g * alpha_g * alpha_s * fabs( v_rel ) /
                    paramsc->particle_diameter / pow( alpha_g, 2.65 );
            }
            else {
                Ksg = 150.0 * pow( alpha_s, 2.0 ) * paramsc->mu2 / alpha_g / pow( paramsc->particle_diameter, 2.0 ) +
                    1.75 * rho_g * alpha_s * fabs( v_rel ) / paramsc->particle_diameter;
            }
            return Ksg * v_rel;
        
        case TANINO:
            // формула для системы цилиндров из
            // Tanino Y., Nepf H.M. Laboratory investigation of mean drag in a random array of rigid, emergent cylinders // J. Hydraul. Eng. - 2008. - V.134, No. 1. - P. 34 - 41.
            return ( ( 106.6975 * alpha_s * paramsc->mu2 / alpha_g / pow( paramsc->particle_diameter, 2.0 ) +
                ( 0.58569 + 4.83831 * alpha_s ) * alpha_s * rho_g * fabs( v_rel ) / alpha_g / paramsc->particle_diameter ) * v_rel );
        case SCHWENDEMAN:
            // формула из
            // D. W. Schwendeman , C. W. Wahle & A. K. Kapila. A study of detonation evolution and structure for a model of compressible two-phase reactive flow // Combustion Theory and Modelling. - 2008. - 12:1. - P. 159 - 204.
            return (paramsc->interface_drag_coef * v_rel) ;

             case HOMENKO:
            // формула из
            // 1999 Ю.П. Хоменко. Математическое моделирование внутрибаллистических процессов в ствольных системах.
            {
            double X, lambda, mu, psi; // Геометрически характеристики пороха
            if (v_ncons[Z0] <= 1.0){
                X = paramsc->X_coef1;
                lambda = paramsc->lambda_coef1;
                mu = paramsc->mu_coef1;
            }
            else{
                X = paramsc->X_coef2;
                lambda = paramsc->lambda_coef2;
                mu = paramsc->mu_coef2;
            }
            psi = X * v_ncons[Z0] * (1 + lambda * v_ncons[Z0] + mu * v_ncons[Z0] * v_ncons[Z0]);           
            double Re = v_ncons[R_GAS] * fabs( v_ncons[V_GAS] - v_ncons[V_DISP] ) * paramsc->particle_diameter / paramsc->mu2;// Число Рейнольдса посчитанное по формуле из Хоменко
            double Sp = 200 * 1 / (1 - psi) * X * (1 + 2 * lambda * v_ncons[Z0] + 3 * v_ncons[Z0] * v_ncons[Z0] ); // Площадь межфазной поверхности
            double A, B, C1, C3, Cd_zern, n;// коэффициенты в формуле
            C3 = 2.33 + 200 * alpha_s / (alpha_g * Re); 
            C1 = 24/Re + 4.4/sqrt(Re) + 0.42;
            if (alpha_g > 0.75){
                if (alpha_g < 0.92)
                    Cd_zern = ((0.92 - alpha_g) * C3 + (alpha_g - 0.75) * C1) / 0.17;
                else
                    Cd_zern = C1;
            

                return 0.125 * Sp* Cd_zern * rho_g * fabs(v_rel) * v_rel;
            }
            else 
            {
                if (alpha_g > 0.4){
                    A = 25;
                    B = 0.29;
                    n = 1;
                }
                else
                    if (alpha_g > 0.12){
                        A = 8.3 / (alpha_g  - 0.068);
                        B = 0.087 / (alpha_g - 0.1);
                        n = 1;
                    }
                    else
                    {
                        A = 10.64 / (alpha_g - 0.054);
                        B = 0.28 / (alpha_g - 0.054);
                        n = 1 + 7.07 * (alpha_g - 0.12);
                    }
            return rho_g * alpha_s * Sp * ( A * pow(alpha_s / ( alpha_g * Re) , n) + B ) * fabs(v_rel) * v_rel;
            }
            }

            break;
        default:
            printf( "\ncalc_friction_force -> fric_force_formula_num should be 1, 2, 3 or 4 only.\n\n" );
            exit( EXIT_FAILURE );
    }
    
}

// Расчет коэффициента сопротивления
// Houim R.W., Oran E.S. A multiphase model for compressible granular-gaseous flows: formulation and initial tests // J. Fluid Mech. - 2016. - V. 789. - P. 166 - 220.
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает коэффициент сопротивления
double calc_Cd( struct ParametersCommon *paramsc, double v_ncons[M] ) {

    double alpha_g = 1 - v_ncons[B_DISP];
    
    double Re = v_ncons[R_GAS] * fabs( v_ncons[V_GAS] - v_ncons[V_DISP] ) * paramsc->particle_diameter / paramsc->mu2;

    if ( alpha_g * Re < 1000.0 ) {
        return 24.0 * ( 1.0 + 0.15 * pow( alpha_g * Re, 0.687 ) ) / ( alpha_g * Re );
    }
    else {
        return 0.44;
    }

}   

// Возвращает в виде вектора начальные параметры в конкретном блоке
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами для 1d2phc случая (in)
// params2d - структура с дополнительными параметрами для 2d2phc случая (in)
// number_of_block[2] - массив с номером блока. В 1d случае имеет значение только number_of_block[0] (in)
// v_ncons[M] - вектор наачльных примитивных переменных в блоке (out)
void calc_initial_parameters_values(const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, int *number_of_block, double v_ncons[M]){
    
    if (paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC){
        for (int i = 0; i < M; i++)
            v_ncons[i] = params1d->block_values[number_of_block[0]][i];

    }
    else if (paramsc->program_name == TWOD2PHC){
        v_ncons[B_DISP] = params2d->block_values[number_of_block[0]][number_of_block[1]][B_DISP_2D];
        v_ncons[R_DISP] = params2d->block_values[number_of_block[0]][number_of_block[1]][R_DISP_2D];
        v_ncons[P_DISP] = params2d->block_values[number_of_block[0]][number_of_block[1]][P_DISP_2D];
        v_ncons[R_GAS] = params2d->block_values[number_of_block[0]][number_of_block[1]][R_GAS_2D];
        v_ncons[P_GAS] = params2d->block_values[number_of_block[0]][number_of_block[1]][P_GAS_2D];
        if (dir == X_DIRECTION){
            v_ncons[V_DISP] = params2d->block_values[number_of_block[0]][number_of_block[1]][V_DISP_2D];
            v_ncons[V_GAS] = params2d->block_values[number_of_block[0]][number_of_block[1]][V_GAS_2D];
        }
        else{
            v_ncons[V_DISP] = params2d->block_values[number_of_block[0]][number_of_block[1]][U_DISP_2D];
            v_ncons[V_GAS] = params2d->block_values[number_of_block[0]][number_of_block[1]][U_GAS_2D];
        }
    }
}


// Расчет скорости компактирования
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами для 1d2phc случая (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает скорость компактирования
double calc_compaction_rate( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block ){
	
    if (paramsc->is_physics && paramsc->is_compaction){
        double beta = calc_configuration_pressure(paramsc, params1d, params2d, dir, v_ncons, number_of_block);
	return v_ncons[B_DISP] * (1.0 - v_ncons[B_DISP]) / paramsc->compaction_viscosity * (v_ncons[P_DISP] - v_ncons[P_GAS] - beta);
    }
    else
	return 0.0;
}

// Расчет конфигурационного давления
// params1d - структура с дополнительными параметрами для 1d2phc случая (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает configuration pressure Beta
double calc_configuration_pressure ( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block ){
    if (paramsc->is_physics && paramsc->is_compaction) { 
        double v_ncons0[M];

        calc_initial_parameters_values(paramsc, params1d, params2d, dir, number_of_block, v_ncons0);

       if (paramsc->pressure_relaxation_compaction && (v_ncons[B_DISP] < paramsc->volume_fraction_compaction))
            return 0.0;
        else
            return - (v_ncons0[P_GAS] - v_ncons0[P_DISP]) * v_ncons[B_DISP] * v_ncons[R_DISP] / v_ncons0[B_DISP] / v_ncons0[R_DISP] * 
	            pow( (2.0 - v_ncons0[B_DISP]) / (2.0 - v_ncons[B_DISP]) , 2.0) * log(1.0 - v_ncons[B_DISP]) / log(1.0 - v_ncons0[B_DISP]);
    }
    else 
        return 0.0;
}

double compaction_energy( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block ){

    if (paramsc->is_physics && (paramsc->is_compaction || paramsc->pressure_relaxation_compaction)){

        double v_ncons0[M];

        calc_initial_parameters_values(paramsc, params1d, params2d, dir, number_of_block, v_ncons0);

        if (paramsc->pressure_relaxation_compaction && (v_ncons[B_DISP] < paramsc->volume_fraction_compaction))
            return 0.0;

        double B =  (v_ncons0[P_GAS] - v_ncons0[P_DISP]) * (2.0 - v_ncons0[B_DISP]) * (2.0 - v_ncons0[B_DISP]) / v_ncons0[B_DISP] / v_ncons0[R_DISP] / log(1.0 - v_ncons0[B_DISP]) * 
	        log( (2.0 - v_ncons0[B_DISP]) / (2.0 - v_ncons[B_DISP]) * pow( (1.0 - v_ncons[B_DISP]), (1.0 - v_ncons[B_DISP]) / (2.0 - v_ncons[B_DISP]) ) / pow( (1.0 - v_ncons0[B_DISP]), (1.0 - v_ncons0[B_DISP]) / (2.0 - v_ncons0[B_DISP])) );
        return B;
    }
    else{
        
        return 0.0;
    }
}

// Расчет химической реакции
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных в центре ячейки (in)
double calc_chemical_reaction(const struct ParametersCommon *paramsc, double v_ncons[M], int ignition_flag){
	
    if (paramsc->is_physics && paramsc->is_burning ){
        if (paramsc->burning_formula_num == SCHWENDEMAN_BURNING && (v_ncons[P_GAS] > paramsc->p_ignition))
            return - paramsc->reaction_rate_prefactor * v_ncons[B_DISP] * v_ncons[R_DISP] * (v_ncons[P_GAS] - paramsc->p_ignition);
        else if (paramsc->burning_formula_num == SEREBRYAKOV_BURNING ){ // нужно ввести температуру воспламенения T_ign
            double psi; // относительная масса сгоревшего пороха
            double psi_s; // относительная масса сгоревшего пороха при z = 1
            double d_psi;// производная psi
            psi_s = paramsc->X_coef1* ( 1 + paramsc->lambda_coef1 + paramsc->mu_coef1);
            double X, lambda, mu;
            if (v_ncons[Z0] <= 1.0){
                X = paramsc->X_coef1;
                lambda = paramsc->lambda_coef1;
                mu = paramsc->mu_coef1;
                psi = X * v_ncons[Z0] * (1 + lambda * v_ncons[Z0] + mu * v_ncons[Z0] * v_ncons[Z0]);
                d_psi = X * (1.0 + 2.0 * v_ncons[Z0] * lambda + 3.0 * mu * v_ncons[Z0] * v_ncons[Z0]);

            }
            else{
                X = paramsc->X_coef2;
                lambda = paramsc->lambda_coef2;
                mu = paramsc->mu_coef2;
                psi = psi_s + X * (v_ncons[Z0] - 1) * ( 1 + lambda * (v_ncons[Z0] - 1));
                d_psi = X * (1 + 2* lambda * (v_ncons[Z0] -1));
            }
            
            double zk = 1.472533333;//zk зашито?
            if (v_ncons[Z0] <= zk){				
		if (paramsc->ignition_condition == 2){ // воспламенение по температуре
		    if (T_gas(paramsc, v_ncons[P_GAS], v_ncons[R_GAS]) > paramsc->T_ignition || ignition_flag == 1){
			return - v_ncons[B_DISP] * v_ncons[R_DISP] * paramsc->U_coef / paramsc->e_coef * pow ( v_ncons[P_GAS], paramsc->nu_coef) * 
			       d_psi / (1.0 - psi);
		    }
		    else
			return 0.0;									
		    }
		if (paramsc->ignition_condition == 1){ // воспламенение по давлению
		    if (v_ncons[P_GAS] > paramsc->p_ignition || ignition_flag == 1){
			return - v_ncons[B_DISP] * v_ncons[R_DISP] * paramsc->U_coef / paramsc->e_coef * pow ( v_ncons[P_GAS], paramsc->nu_coef) * 
			d_psi / (1.0 - psi);

		    }
		    else
			return 0.0;							 

		}				
					
	    }
            else
                return 0.0;
    }
        else{
            printf("\n calc_chemical_reaction: wrong paramsc->burning_formula_num is used. \n");
            exit(EXIT_FAILURE);
        }
    }
    else
	return 0.0;

}

// Расчет теплопереноса
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных в центре ячейки (in)
double calc_heat_transfer(const struct ParametersCommon *paramsc, double v_ncons[M]){
	
    if (paramsc->is_physics && paramsc->is_heat_transfer){
	return paramsc->heat_transfer_coefficient * (T_gas(paramsc, v_ncons[P_GAS], v_ncons[R_GAS]) - T_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]));
	printf("\n %e",paramsc->heat_transfer_coefficient );
	}
    else
	return 0.0;

}
//Расчёт вязкого нагрева
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных в центре ячейки (in)