// physics_solver.cc
// Учет всей "физики" в ячейке
// (c) Уткин Павел, 2014
// Создан: 7 января 2014 г.

#include "physics_solver.h"
#include "utils.h"
#include "source_terms.h"
#include "math_utils.h"


// Решение системы обыкновенных дифференциальных уравнений, описывающих межфазное взаимодействие и иную "физику"
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// dt - текущий шаг по времени (in)
// v_ncons[M] - вектор примитивных переменных в ячейке, на выходе - измененный в результате учета "физики" (in/out)
// v_ncons_prev[M] - вектор примитивных переменных в предыдущей ячейке, на выходе - измененный в результате учета "физики" (in/out)
// v_ncons_next[M] - вектор примитивных переменных в следующей ячейке, на выходе - измененный в результате учета "физики" (in/out)
void physics_solver( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double dt, double *v_ncons, double *v_ncons_next, double *v_ncons_prev, int *number_of_block, double *beta, double *B, double *C, int n, int ignition_flag, int initial_ignition, int cell_number, int left_ghost ) {
    
    int i_component; // индекс компонент векторов решения и правых частей
    int i_step; // индекс подшагов интегрирования по времени
    double sub_dt = dt / paramsc->substeps_num; // подшаг интегрирования

    // если дисперсная фаза отсутствует, то "физику" рассчитывать не нужно
    if ( v_ncons[B_DISP] < paramsc->eps_disp_abs )
        return;

    double v_cons[M]; // вектор "консервативных" переменных в ячейке 
    double right_hand_side_terms[M]; // вектор правых частей СОДУ

    // инициализация вектора правых частей нулями
    for ( i_component = 0; i_component < n; i_component++ )
        right_hand_side_terms[i_component] = 0.0;


    if (paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC){


        // переходим к "консервативным" переменным, поскольку члены, описывающие межфазное взаимодействие,
        // записываются в системе, имеющей "консервативную" форму
        convert_noncons_to_cons( paramsc, v_ncons, v_cons, params1d->number_of_scalars );

        for ( i_step = 0; i_step < paramsc->substeps_num; i_step++ ) { // цикл интегрирования по времени
	    calc_right_hand_side_terms( paramsc, params1d, params2d, dir, v_ncons, right_hand_side_terms, number_of_block, M, ignition_flag, initial_ignition ); // расчет текущего вектора правых частей
            if (paramsc->numerical_method_physics == EXPLICIT_EULER){
                for ( i_component = 0; i_component < n; i_component++ ) { // цикл по компонентам вектора
                    // явный метод Эйлера
	            v_cons[i_component] = v_cons[i_component] + sub_dt * right_hand_side_terms[i_component];	
	            //runge_kutta_predictor_2_nd_order(paramsc, params1d, params2d, v_cons, right_hand_side_terms, sub_dt, number_of_block);
                }
                if (cell_number == left_ghost -2){
                            double friction_force = calc_friction_force( paramsc, v_ncons ); // силовое межфазное взаимодействие

    

                            double compaction_rate = calc_compaction_rate(paramsc, params1d, params2d, dir, v_ncons, number_of_block); // скорость компактирования
	
                            double beta = calc_configuration_pressure(paramsc, params1d, params2d, dir, v_ncons, number_of_block);
	
                            double B = compaction_energy(paramsc, params1d, params2d, dir, v_ncons, number_of_block);
        if (!paramsc->is_compaction){
            beta = 0.0;
            B = 0.0;
        }

        double C = calc_chemical_reaction(paramsc, v_ncons, ignition_flag);
        printf("\n phys P %lf, fr %lf, C %lf", v_cons[P_GAS], friction_force, C); 
                }
            }
            else
                impicit_euler_full_1d(paramsc, params1d, v_cons, sub_dt, number_of_block[0], M1D);

            // обновляем вектор примитивных переменных
            convert_cons_to_noncons( paramsc, v_cons, v_ncons, params1d->number_of_scalars );
        }
    }
    else{
   
        // переходим к "консервативным" переменным, поскольку члены, описывающие межфазное взаимодействие,
        // записываются в системе, имеющей "консервативную" форму
        double b2 = 1 - v_ncons[B_DISP_2D]; // объемная доля газовой фазы 
        v_cons[B_DISP_2D] = v_ncons[B_DISP_2D]; // b1, объемная доля дисперсной фазы 
        v_cons[R_DISP_2D] = v_ncons[B_DISP_2D] * v_ncons[R_DISP_2D]; // b1 * r1
        v_cons[V_DISP_2D] = v_cons[R_DISP_2D] * v_ncons[V_DISP_2D]; // b1 * r1 * v1
        v_cons[U_DISP_2D] = v_cons[R_DISP_2D] * v_ncons[U_DISP_2D]; // b1 * r1 * u1
        v_cons[P_DISP_2D] = v_ncons[B_DISP_2D] * v_ncons[R_DISP_2D] * ( 0.5 * pow( v_ncons[V_DISP_2D], 2.0 ) + 0.5 * pow( v_ncons[U_DISP_2D], 2.0 ) + v_ncons[P_DISP_2D] / v_ncons[R_DISP_2D] /
            ( paramsc->g1 - 1.0 ) + paramsc->g1 * paramsc->p01 / ( paramsc->g1 - 1.0 ) / v_ncons[R_DISP_2D] ); // b1 * r1 * E1 
        v_cons[R_GAS_2D] = b2 * v_ncons[R_GAS_2D]; // b2 * r2
        v_cons[V_GAS_2D] = v_cons[R_GAS_2D] * v_ncons[V_GAS_2D]; // b2 * r2* v2
        v_cons[U_GAS_2D] = v_cons[R_GAS_2D] * v_ncons[U_GAS_2D]; // b2 * r2* u2
        v_cons[P_GAS_2D] = v_cons[R_GAS_2D] * ( 0.5 * pow( v_ncons[V_GAS_2D], 2.0 ) + 0.5 * pow( v_ncons[U_GAS_2D], 2.0 ) + ( v_ncons[P_GAS_2D] + paramsc->g2 * paramsc->p02 )/ v_ncons[R_GAS_2D] /
		( paramsc->g2 - 1.0 ) / (1.0 + paramsc->b_virial * v_ncons[R_GAS_2D])) ; // b2 * r2 * E2
        
        for ( i_step = 0; i_step < paramsc->substeps_num; i_step++ ) { // цикл интегрирования по времени
	    calc_right_hand_side_terms( paramsc, params1d, params2d, dir, v_ncons, right_hand_side_terms, number_of_block, M2D, ignition_flag, initial_ignition ); // расчет текущего вектора правых частей
            if (paramsc->numerical_method_physics == EXPLICIT_EULER){
                for ( i_component = 0; i_component < n; i_component++ ) { // цикл по компонентам вектора
                    // явный метод Эйлера
	            v_cons[i_component] = v_cons[i_component] + sub_dt * right_hand_side_terms[i_component];	
	            //runge_kutta_predictor_2_nd_order(paramsc, params1d, params2d, v_cons, right_hand_side_terms, sub_dt, number_of_block);
                }
            }
            else
                impicit_euler_full_2d(paramsc, params2d, v_cons, sub_dt, number_of_block, n);
            // обновляем вектор примитивных переменных
            b2 = 1 - v_cons[B_DISP_2D]; // объемная доля газовой фазы 
            v_ncons[B_DISP_2D] = v_cons[B_DISP_2D]; // b1, объемная доля дисперсной фазы
            v_ncons[R_DISP_2D] = v_cons[R_DISP_2D] / v_cons[B_DISP_2D]; // r1, плотность дисперсной фазы
            v_ncons[V_DISP_2D] = v_cons[V_DISP_2D] / v_cons[R_DISP_2D]; // v1, скорость дисперсной фазы
            v_ncons[U_DISP_2D] = v_cons[U_DISP_2D] / v_cons[R_DISP_2D]; // u1, скорость дисперсной фазы
            v_ncons[P_DISP_2D] = ( paramsc->g1 - 1.0 ) * ( v_cons[P_DISP_2D] / v_cons[B_DISP_2D] - 0.5 * ( pow( v_cons[V_DISP_2D], 2.0 ) + pow( v_cons[U_DISP_2D], 2.0 ) ) / v_cons[B_DISP_2D] / v_cons[R_DISP_2D] ) -
                paramsc->g1 * paramsc->p01; // p1, давление дисперсной фазы
            v_ncons[R_GAS_2D] = v_cons[R_GAS_2D] / b2; // r2, плотность газовой фазы
            v_ncons[V_GAS_2D] = v_cons[V_GAS_2D] / v_cons[R_GAS_2D]; // v2, скорость газовой фазы
            v_ncons[U_GAS_2D] = v_cons[U_GAS_2D] / v_cons[R_GAS_2D]; // u2, скорость газовой фазы
            v_ncons[P_GAS_2D] = ( paramsc->g2 - 1.0 ) * (1.0 + paramsc->b_virial * v_cons[R_GAS_2D] / b2) * ( v_cons[P_GAS_2D] / b2 - 0.5 * ( pow( v_cons[V_GAS_2D], 2.0 ) + pow( v_cons[U_GAS_2D], 2.0 ) ) / b2 / v_cons[R_GAS_2D] ) -
                paramsc->g2 * paramsc->p02; // p2, давление газовой фазы 
        }
    
    }
    if (!paramsc->pressure_relaxation_compaction ){
        *beta = calc_configuration_pressure(paramsc, params1d, params2d, dir, v_ncons, number_of_block);
        *B = compaction_energy(paramsc, params1d, params2d, dir, v_ncons, number_of_block);
    }
    *C = calc_chemical_reaction(paramsc, v_ncons, ignition_flag);
}

// Расчет вектора правых частей
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// right_hand_side_terms[M] - текущий вектор правых частей (out)
void calc_right_hand_side_terms( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double *v_ncons, double *right_hand_side_terms, int *number_of_block, int n, int ignition_flag, int initial_ignition ) {

    if (paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC ) {

        double friction_force = calc_friction_force( paramsc, v_ncons ); // силовое межфазное взаимодействие

      //  friction_force = 0.5 * friction_force;

        double compaction_rate = calc_compaction_rate(paramsc, params1d, params2d, dir, v_ncons, number_of_block); // скорость компактирования
	
        double beta = calc_configuration_pressure(paramsc, params1d, params2d, dir, v_ncons, number_of_block);
	
        double B = compaction_energy(paramsc, params1d, params2d, dir, v_ncons, number_of_block);
        if (!paramsc->is_compaction){
            beta = 0.0;
            B = 0.0;
        }

        double C = calc_chemical_reaction(paramsc, v_ncons, ignition_flag);

        double heat_transfer = calc_heat_transfer(paramsc, v_ncons);
	double additional_heat;//энергия поджига
	double mass_input;// внесённая при поджиге масса
	if (initial_ignition == 1 ){
	    additional_heat = params1d->additional_heat;
            mass_input = params1d->mass_input;
	}
	else{
	    additional_heat = 0;
	    mass_input = 0;
	}
        // правая часть уравнения компактирования
        right_hand_side_terms[B_DISP] = compaction_rate + C / v_ncons[R_DISP];

        // правая часть закона сохранения массы дисперсной фазы
        right_hand_side_terms[R_DISP] = C;
        // правая часть закона сохранения импульса дисперсной фазы
        if ( paramsc->porous_body != true ) // для уменьшения подвижности дисперсной фазы в случае пористого каркаса не учитываем силу трения 
	    right_hand_side_terms[V_DISP] = - friction_force + C * v_ncons[V_DISP] + 0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP]);

        // правая часть закона сохранения энергии дисперсной фазы
        if ( paramsc->porous_body != true ){ // для уменьшения подвижности дисперсной фазы в случае пористого каркаса не учитываем силу трения
		
	    right_hand_side_terms[P_DISP] = - (v_ncons[P_GAS] + beta) * compaction_rate + (- friction_force + 0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP])) * calc_u_i( paramsc, v_ncons ) 
                + ( e_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]) + 0.5 * v_ncons[V_DISP] * v_ncons[V_DISP]+ B + beta / v_ncons[R_DISP] ) * C  + heat_transfer ;
        //    right_hand_side_terms[P_DISP] = ( 0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP])) * calc_u_i( paramsc, v_ncons ) + ( e_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]) + 0.5 * v_ncons[V_DISP] * v_ncons[V_DISP] ) * C ;

        }
        // правая часть закона сохранения массы газовой фазы
        right_hand_side_terms[R_GAS] = - right_hand_side_terms[R_DISP] + mass_input;

        // правая часть закона сохранения импульса газовой фазы
        right_hand_side_terms[V_GAS] = - right_hand_side_terms[V_DISP];

        // правая часть закона сохранения энергии газовой фазы
		// условие воспламенения
		if (paramsc->ignition_condition == 2){//воспламенение по температуре
		    if (T_gas(paramsc, v_ncons[P_GAS], v_ncons[R_GAS]) > paramsc->T_ignition || ignition_flag == 1)
		        right_hand_side_terms[P_GAS] = v_ncons[P_GAS] * compaction_rate - (- friction_force + 0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP])) * calc_u_i( paramsc, v_ncons ) 
			- ( e_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]) + 0.5 * v_ncons[V_DISP] * v_ncons[V_DISP] + B + beta / v_ncons[R_DISP] + paramsc->heat_release) * C + additional_heat*mass_input  - heat_transfer;	
		    else
			right_hand_side_terms[P_GAS] = v_ncons[P_GAS] * compaction_rate - (- friction_force + 0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP])) * calc_u_i( paramsc, v_ncons ) 
			- ( e_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]) + 0.5 * v_ncons[V_DISP] * v_ncons[V_DISP]) * C + additional_heat*mass_input   - heat_transfer;							  
		}
		if (paramsc->ignition_condition == 1){ // воспламенение по давлению
		    if (v_ncons[P_GAS] > paramsc->p_ignition || ignition_flag == 1)
			right_hand_side_terms[P_GAS] = v_ncons[P_GAS] * compaction_rate - (- friction_force + 0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP])) * calc_u_i( paramsc, v_ncons ) 
			- ( e_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]) + 0.5 * v_ncons[V_DISP] * v_ncons[V_DISP] + B + beta / v_ncons[R_DISP] + paramsc->heat_release ) * C + additional_heat*mass_input  - heat_transfer;	
		    else
			right_hand_side_terms[P_GAS] = v_ncons[P_GAS] * compaction_rate - (- friction_force + 0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP])) * calc_u_i( paramsc, v_ncons ) 
			- ( e_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]) + 0.5 * v_ncons[V_DISP] * v_ncons[V_DISP] + B + beta / v_ncons[R_DISP]) * C + additional_heat*mass_input  - heat_transfer;
						  }
		
  //      right_hand_side_terms[P_GAS] =  -  0.5 * C * (v_ncons[V_GAS] - v_ncons[V_DISP]) * calc_u_i( paramsc, v_ncons ) - ( e_disp(paramsc, v_ncons[P_DISP], v_ncons[R_DISP]) + 0.5 * v_ncons[V_DISP] * v_ncons[V_DISP] + paramsc->heat_release) * C;
		// перенос скаляров
        for (int i = 0; i < params1d->number_of_scalars; i++){ 
            if (i == 0 && v_ncons[Z0] <= paramsc->zk)
		// условие воспламенения
	        if (paramsc->ignition_condition == 2){ //воспламенение по температуре
		    if (T_gas(paramsc, v_ncons[P_GAS], v_ncons[R_GAS]) > paramsc->T_ignition || ignition_flag == 1)
		        right_hand_side_terms[Z0 + i] = C * v_ncons[Z0 + i] + v_ncons[B_DISP] * v_ncons[R_DISP] * ( paramsc->U_coef / paramsc->e_coef * pow ( v_ncons[P_GAS], paramsc->nu_coef) );	
		        else
			    right_hand_side_terms[Z0 + i] = C * v_ncons[Z0 + i] ;
			}
			    if (paramsc->ignition_condition == 1) // воспламенение по давлению
			        if (v_ncons[P_GAS] > paramsc->p_ignition || ignition_flag == 1)
				    right_hand_side_terms[Z0 + i] = C * v_ncons[Z0 + i] + v_ncons[B_DISP] * v_ncons[R_DISP] * ( paramsc->U_coef / paramsc->e_coef * pow ( v_ncons[P_GAS], paramsc->nu_coef) );		
				else
				    right_hand_side_terms[Z0 + i] = C * v_ncons[Z0 + i] ;				

        }

    }
    else if (paramsc->program_name == TWOD2PHC){
    
        double v_ncons_x[M];
        double v_ncons_y[M];
    
        convert_2D_to_1D_vector(v_ncons, X_DIRECTION, v_ncons_x);
        convert_2D_to_1D_vector(v_ncons, Y_DIRECTION, v_ncons_y);

        double friction_force_1 = calc_friction_force( paramsc, v_ncons_x ); // силовое межфазное взаимодействие
        double friction_force_2 = calc_friction_force( paramsc, v_ncons_y ); // силовое межфазное взаимодействие

        if (friction_force_1 > 0.0 || friction_force_2 > 0.0){
     //       printf("friciton force\n");
        }

        double compaction_rate = calc_compaction_rate(paramsc, params1d, params2d, dir, v_ncons_x, number_of_block); // скорость компактирования
	
        double beta = calc_configuration_pressure(paramsc, params1d, params2d, dir, v_ncons_x, number_of_block);
	
        double B = compaction_energy(paramsc, params1d, params2d, dir, v_ncons_x, number_of_block);
        if (!paramsc->is_compaction){
            beta = 0.0;
            B = 0.0;
            compaction_rate = 0.0;
        }

        double C = calc_chemical_reaction(paramsc, v_ncons_x, ignition_flag);

        double heat_transfer = calc_heat_transfer(paramsc, v_ncons_x);
		printf("\n %e", heat_transfer);
        // правая часть уравнения компактирования
        right_hand_side_terms[B_DISP_2D] = compaction_rate + C / v_ncons[R_DISP];

        // правая часть закона сохранения массы дисперсной фазы
        right_hand_side_terms[R_DISP_2D] = C;

        // правая часть закона сохранения импульса дисперсной фазы
        if ( paramsc->porous_body != true ) { // для уменьшения подвижности дисперсной фазы в случае пористого каркаса не учитываем силу трения 
	    right_hand_side_terms[V_DISP_2D] = - friction_force_1 + C * v_ncons[V_DISP_2D] + 0.5 * C * (v_ncons[V_GAS_2D] - v_ncons[V_DISP_2D]);
            right_hand_side_terms[U_DISP_2D] = - friction_force_2 + C * v_ncons[U_DISP_2D] + 0.5 * C * (v_ncons[U_GAS_2D] - v_ncons[U_DISP_2D]);
        }
        // правая часть закона сохранения энергии дисперсной фазы
        if ( paramsc->porous_body != true ){ // для уменьшения подвижности дисперсной фазы в случае пористого каркаса не учитываем силу трения
		
	    right_hand_side_terms[P_DISP_2D] = - (v_ncons[P_GAS_2D] + beta) * compaction_rate + (- friction_force_1 + 0.5 * C * (v_ncons[V_GAS_2D] - v_ncons[V_DISP_2D])) * calc_u_i( paramsc, v_ncons_x ) +
                (- friction_force_2 + 0.5 * C * (v_ncons[U_GAS_2D] - v_ncons[U_DISP_2D])) * calc_u_i( paramsc, v_ncons_y ) + 
                ( e_disp(paramsc, v_ncons[P_DISP_2D], v_ncons[R_DISP_2D]) + 0.5 * v_ncons[V_DISP_2D] * v_ncons[V_DISP_2D] + 0.5 * v_ncons[U_DISP_2D] * v_ncons[U_DISP_2D]) * C  + heat_transfer;

        }

        // правая часть закона сохранения массы газовой фазы
        right_hand_side_terms[R_GAS_2D] = - right_hand_side_terms[R_DISP_2D];

        // правая часть закона сохранения импульса газовой фазы
        right_hand_side_terms[V_GAS_2D] = - right_hand_side_terms[V_DISP_2D];
        right_hand_side_terms[U_GAS_2D] = - right_hand_side_terms[U_DISP_2D];

        // правая часть закона сохранения энергии газовой фазы
        right_hand_side_terms[P_GAS_2D] = v_ncons[P_GAS_2D] * compaction_rate - (- friction_force_1 + 0.5 * C * (v_ncons[V_GAS_2D] - v_ncons[V_DISP_2D])) * calc_u_i( paramsc, v_ncons_x ) 
            - (- friction_force_2 + 0.5 * C * (v_ncons[U_GAS_2D] - v_ncons[U_DISP_2D])) * calc_u_i( paramsc, v_ncons_y ) 
            - ( e_disp(paramsc, v_ncons[P_DISP_2D], v_ncons[R_DISP_2D]) + 0.5 * v_ncons[V_DISP_2D] * v_ncons[V_DISP_2D] + 0.5 * v_ncons[U_DISP_2D] * v_ncons[U_DISP_2D] + B + beta / v_ncons[R_DISP_2D] + paramsc->heat_release) * C  - heat_transfer;
    }
}

// Решение СОДУ методом Рунге-Кутты (предиктор-корректор) 2-го порядка
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// v_cons[M] - текущий вектор консервативных переменных в ячейке (in / out)
// right_hand_side_terms[M] - текущий вектор правых частей (in)
// sub_dt - шаг по времени (in)
// i_step_current - номер текущего шага по пространству (in)
void runge_kutta_predictor_2_nd_order (struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double v_cons[M], double right_hand_side_terms[M], double sub_dt, int *number_of_block, int ignition_flag, int initial_ignition ){

    double v_cons_predicted[M];

    for (int i_step = 0; i_step < M; i_step++)
        v_cons_predicted[i_step] = v_cons[i_step] + sub_dt * right_hand_side_terms[i_step];

    double right_hand_side_from_prediction[M];
    double v_ncons_predicted[M];
    convert_cons_to_noncons(paramsc, v_cons_predicted, v_ncons_predicted, params1d->number_of_scalars );
    calc_right_hand_side_terms(paramsc, params1d, params2d, dir, v_ncons_predicted, right_hand_side_from_prediction, number_of_block, M1D, ignition_flag, initial_ignition);

    for (int i_step = 0; i_step < M; i_step++)
	v_cons[i_step] = v_cons[i_step] + sub_dt / 2.0 * (right_hand_side_terms[i_step] + right_hand_side_from_prediction[i_step]);
}

// Решение полной системы уравнений (7) по неявной схеме Эйлера методом Ньютона - задача 1d2phc
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// v_cons[M1D] - текущий вектор консервативных переменных в ячейке (in / out)
// sub_dt - шаг по времени (in)
// number_of_block - номер блока (in)
// n - реальный размер векторов (in)
void impicit_euler_full_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double v_cons[M], double sub_dt, int number_of_block, int n ){

    double a, b, c, d, e, f, g; // Величины, отвечающие консервативным переменным (по порядку от B_DISP до P_GAS)
    double fi[M]; // вектор fi = 0
    double matrix[M][M]; // матрица частных производных

    double F, dF[M]; // Скорость коспактирования F, ее производная по всем переменным a ... g
    double momentum, dM[M]; // Residual momentum exchange M, ее производная по всем переменным a ... g
    double E_disp, dE_disp[M]; // Residual energy exchange E дисперсной фазы , ее производная по всем переменным a ... g 
    double E, dE[M]; // Residual energy exchange E газовой фазы , ее производная по всем переменным a ... g
    double p, dp[M]; // давление газа, производная давления по всем переменным a ... g
    double beta, dbeta[M];
    double C, dC[M];
    double H, dH[M];
    double pressure_disp, dp_disp[M];
    double temp, temp_disp;
    double B, dB[M];

    double du[M], b_vector[M], u[M], u_prev[M]; // вектор решение du, ветор правой части b_vector, u - вектор переменных a ... g на текущем шаге, u_prev - вектор переменных a ... g на предыдущем шаге

    // Присваиваем начальные значения векторам u и u_prev
    for (int i = 0; i < n; i++){
	u[i] = v_cons[i];
	u_prev[i] = v_cons[i];
    }

    // Переприсваиваем значения параметров более удобным для использования переменным
    double muc = paramsc->compaction_viscosity; // взякость компактирования
    double gamma_disp = paramsc->g1; // показатель адиабаты дисперсной фазы (УРС)
    double gamma = paramsc->g2; // показатель адиабаты газовой фазы (УРС)
    double P0_disp = paramsc->p01; // параметр P0 дисперсной фазы (УРС)
    double p0 = params1d->block_values[number_of_block][P_GAS]; // начальное значение давления газовой фазы
    double p0_disp = params1d->block_values[number_of_block][P_DISP]; // начальное значение давления дисперсной фазы
    double r0_disp = params1d->block_values[number_of_block][R_DISP]; // начальное значение плотности дисперсной фазы
    double b0_disp = params1d->block_values[number_of_block][B_DISP]; // начальное значение объемной доли дисперсной фазы
    double A_constant = (p0 - p0_disp) / b0_disp / r0_disp * (2.0 - b0_disp) * (2.0 - b0_disp) / log(1.0 - b0_disp); // константа, присутствующая в определении скорости компактирования F, в частности в параметре бета
    double delta = paramsc->interface_drag_coef;
    double b_vir = paramsc->b_virial;
    double q = paramsc->heat_release;
    double sigma = paramsc->reaction_rate_prefactor;

    double local_epsilon;

    do{

	local_epsilon = 0.0;

	a = u[0];
	b = u[1];
	c = u[2];
	d = u[3];
	e = u[4];
	f = u[5];
	g = u[6];
			
        beta = -A_constant * b * log(1.0 - a) / (2.0 - a) / (2.0 - a);
        p = (g - f * f / 2.0 / e) * (gamma - 1.0) / (1.0 - a) * (1.0 + b_vir * e / (1.0 - a));
        pressure_disp = (d / b - 0.5 * c * c / b / b) * (gamma_disp - 1.0) * b / a - gamma_disp * P0_disp;
        F = a * (1.0 - a) / muc * (pressure_disp - p - beta);
        B = A_constant * log ((2.0 - b0_disp) / (2.0 - a) * pow( (1.0 - a), (1.0 - a) / (2.0 - a)) / pow((1.0 - b0_disp), (1.0 - b0_disp) / (2.0 - b0_disp ) ));

        if (!paramsc->is_compaction){
            B = 0.0;
            beta = 0.0;
            F = 0.0;
        }

        if (paramsc->is_burning && p > paramsc->p_ignition)
            C = -sigma * b * (p - paramsc->p_ignition);
        else
            C = 0.0;

        temp = T_gas(paramsc, p, e / (1.0 - a));
        temp_disp = T_disp(paramsc, pressure_disp, b / a);
        H = paramsc->heat_transfer_coefficient * (temp - temp_disp);
        if (!paramsc->is_heat_transfer)
            H = 0.0;

        momentum = (delta + 0.5 * C) * (f/e - c/b) + C * c / b;
        if (!paramsc->is_friction_force)
            momentum = 0.0;
	E_disp = (momentum - C * c / b) * c / b + d / b * C + H;
        E = (momentum - C * c / b) * c / b + (d / b + B + beta * a / b + q) * C + H;

        dp[B_DISP] =   (g - f * f / 2.0 / e) * (gamma - 1.0) / (1.0 - a) / (1.0 - a) * (1.0 + 2.0 * b_vir * e / (1.0 - a));
	dp[R_DISP] =   0.0;
	dp[V_DISP] =   0.0;
	dp[P_DISP] =   0.0;
	dp[R_GAS]  =   f * f / 2.0 / e / e * (gamma - 1.0) / (1.0 - a) * (1.0 + b_vir * e / (1.0 - a)) + 
            (g - 0.5 * f * f / e ) * (gamma - 1.0) * b_vir / (1.0 - a) / (1.0 - a);
	dp[V_GAS]  = - f / e * (gamma - 1.0) / (1.0 - a) * (1.0 + 1.0 * b_vir * e / (1.0 - a));
	dp[P_GAS]  =   (gamma - 1.0) / (1.0 - a) * (1.0 + 1.0 * b_vir * e / (1.0 - a));

        dbeta[0] = A_constant * b /  (2.0 - a) / (2.0 - a) * ( 1.0 / (a - 1.0) + 2.0 * log(1.0 - a) / (2.0 - a));
	dbeta[1] = A_constant * log(1.0 - a) / (2.0 - a) / (2.0 - a);
	dbeta[2] = 0.0;
	dbeta[3] = 0.0;
	dbeta[4] = 0.0;
	dbeta[5] = 0.0;
        dbeta[6] = 0.0;

        dp_disp[B_DISP] = - (d / b - 0.5 * c * c / b / b) * (gamma_disp - 1.0) * b / a / a;
        dp_disp[R_DISP] = 0.5 * c * c / b / b * (gamma_disp - 1.0) / a;
        dp_disp[V_DISP] = - c / b * (gamma_disp - 1.0) / a;
        dp_disp[P_DISP] = (gamma_disp - 1.0) / a;
        dp_disp[R_GAS]  = 0.0;
        dp_disp[V_GAS]  = 0.0;
        dp_disp[P_GAS]  = 0.0;

        dF[B_DISP] = (1.0 - 2.0 * a) / muc * (pressure_disp - p - beta) + a * (1.0 - a) / muc * (dp_disp[B_DISP] - dp[B_DISP] - dbeta[B_DISP]);
	dF[R_DISP] = a * (1.0 - a) / muc * (dp_disp[R_DISP] - dp[R_DISP] - dbeta[R_DISP]);
	dF[V_DISP] = a * (1.0 - a) / muc * (dp_disp[V_DISP] - dp[V_DISP] - dbeta[V_DISP]);
	dF[P_DISP] = a * (1.0 - a) / muc * (dp_disp[P_DISP] - dp[P_DISP] - dbeta[P_DISP]);
	dF[R_GAS]  = a * (1.0 - a) / muc * (dp_disp[R_GAS] - dp[R_GAS] - dbeta[R_GAS]);
	dF[V_GAS]  = a * (1.0 - a) / muc * (dp_disp[V_GAS] - dp[V_GAS] - dbeta[V_GAS]);
	dF[P_GAS]  = a * (1.0 - a) / muc * (dp_disp[P_GAS] - dp[P_GAS] - dbeta[P_GAS]);

        dB[B_DISP] = beta / b;
        dB[R_DISP] = 0.0;
	dB[V_DISP] = 0.0;
	dB[P_DISP] = 0.0;
	dB[R_GAS]  = 0.0;
	dB[V_GAS]  = 0.0;
	dB[P_GAS]  = 0.0;

        dC[B_DISP] = - sigma * b * dp[B_DISP];
        dC[R_DISP] = - sigma * b * dp[R_DISP] - sigma * (p - paramsc->p_ignition);
	dC[V_DISP] = - sigma * b * dp[V_DISP];
	dC[P_DISP] = - sigma * b * dp[P_DISP];
	dC[R_GAS]  = - sigma * b * dp[R_GAS];
	dC[V_GAS]  = - sigma * b * dp[V_GAS];
	dC[P_GAS]  = - sigma * b * dp[P_GAS];

        /*
        dH[B_DISP] =   paramsc->heat_transfer_coefficient * P0_disp / paramsc->specific_heat_cv1 / b;
        dH[R_DISP] =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv1 / b / b * (d - c * c / b - P0_disp * a);
	dH[V_DISP] =   paramsc->heat_transfer_coefficient * c / paramsc->specific_heat_cv1 / b / b;
	dH[P_DISP] = - paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv1 / b;
	dH[R_GAS]  = - paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e / e * (g - f * f / e);
	dH[V_GAS]  = - paramsc->heat_transfer_coefficient * f / paramsc->specific_heat_cv2 / e / e;
	dH[P_GAS]  =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e;
        */
        dH[B_DISP] = - paramsc->heat_transfer_coefficient * (pressure_disp + P0_disp) / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b 
            - paramsc->heat_transfer_coefficient * dp_disp[B_DISP] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;
        dH[R_DISP] =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv1 / b / b * (pressure_disp + P0_disp) * a
            - paramsc->heat_transfer_coefficient * dp_disp[R_DISP] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;
	dH[V_DISP] = - paramsc->heat_transfer_coefficient * dp_disp[V_DISP] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;
        dH[P_DISP] = - paramsc->heat_transfer_coefficient * dp_disp[P_DISP] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;
	dH[R_GAS]  = - paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e / e * (g - (f * f) / e)
            - paramsc->heat_transfer_coefficient * dp_disp[R_GAS] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;
	dH[V_GAS]  = - paramsc->heat_transfer_coefficient * f / paramsc->specific_heat_cv2 / e / e
            - paramsc->heat_transfer_coefficient * dp_disp[V_GAS] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;
        dH[P_GAS]  =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e 
            - paramsc->heat_transfer_coefficient * dp_disp[P_GAS] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

        for (int i = 0; i < n; i++){
            if (!paramsc->is_compaction){
                dbeta[i] = 0.0;
                dB[i] = 0.0;
                dF[i] = 0.0;
            }
            if (!paramsc->is_burning || C == 0.0 ){
                dC[i] = 0.0;
            }
            if (!paramsc->is_heat_transfer)
                dH[i] = 0.0;
        }
  
        dM[B_DISP] =   dC[B_DISP] * c / b + 0.5 * dC[B_DISP] * (f / e - c / b);
        dM[R_DISP] =   dC[R_DISP] * c / b - C * c / b / b + delta * c / b / b + 0.5 * dC[R_DISP] * (f / e - c / b) + 0.5 * C * c / b / b;
        dM[V_DISP] =   C / b + dC[V_DISP] * c / b - delta / b + 0.5 * dC[V_DISP] * (f / e - c / b) - 0.5 * C / b;
        dM[P_DISP] =   dC[P_DISP] * c / b + 0.5 * dC[P_DISP] * (f / e - c / b);
        dM[R_GAS]  =   dC[R_GAS] * c / b - delta * f / e / e + 0.5 * dC[R_GAS] * f / e - 0.5 * C * f / e / e - 0.5 * dC[R_GAS] * c / b;
        dM[V_GAS]  =   dC[V_GAS] * c / b + delta / e + 0.5 * dC[V_GAS] * f / e + 0.5 * C / e - 0.5 * dC[V_GAS] * c / b;
	dM[P_GAS]  =   dC[P_GAS] * c / b + 0.5 * dC[P_GAS] * (f / e - c / b);

        for (int i = 0; i < n; i++){
            if (!paramsc->is_friction_force)
                dM[i] = 0.0;
        }

        dE_disp[B_DISP] =   d / b * dC[B_DISP] + dM[B_DISP] * c / b - dC[B_DISP] * c * c / b / b + dH[B_DISP];
        dE_disp[R_DISP] = - d * C / b / b + d / b * dC[R_DISP] + dM[R_DISP] * c / b - momentum * c / b / b - dC[R_DISP] * c * c / b / b + C * 2.0 * c * c / b / b / b + dH[R_DISP];
        dE_disp[V_DISP] =   d / b * dC[V_DISP] + dM[V_DISP] * c / b + momentum / b - dC[V_DISP] * c * c / b / b - C * 2.0 * c / b / b + dH[V_DISP];
        dE_disp[P_DISP] =   C / b + d / b * dC[P_DISP] + dM[P_DISP] * c / b - dC[P_DISP] * c * c / b / b + dH[P_DISP];
        dE_disp[R_GAS] =    d / b * dC[R_GAS] + (dM[R_GAS] - dC[R_GAS] * c / b) * c / b + dH[R_GAS];
	dE_disp[V_GAS] =    d / b * dC[V_GAS] + (dM[V_GAS] - dC[V_GAS] * c / b) * c / b + dH[V_GAS];
	dE_disp[P_GAS] =    d / b * dC[P_GAS] + (dM[P_GAS] - dC[P_GAS] * c / b) * c / b + dH[P_GAS];

        dE[B_DISP] =   d / b * dC[B_DISP] + B * dC[B_DISP] + C * dB[B_DISP] + beta * C / b + dbeta[B_DISP] * a / b * C + a * beta / b * dC[B_DISP] + q * dC[B_DISP] + dM[B_DISP] * c / b - dC[B_DISP] * c * c / b / b + dH[B_DISP];
        dE[R_DISP] = - d * C / b / b + d / b * dC[R_DISP] + B * dC[R_DISP] + C * dB[R_DISP] + dbeta[R_DISP] * a / b * C - beta * a / b / b * C + beta  * a / b * dC[R_DISP] + dM[R_DISP] * c / b - momentum * c / b / b - dC[R_DISP] * c * c / b / b + C * 2.0 * c * c / b / b / b + dH[R_DISP] + q * dC[R_DISP];
        dE[V_DISP] =   d / b * dC[V_DISP] + B * dC[V_DISP] + C * dB[V_DISP] + dbeta[V_DISP] * a / b * C + beta * a / b * dC[V_DISP] + dM[V_DISP] * c / b + momentum / b - dC[V_DISP] * c * c / b / b - C * 2.0 * c / b / b + dH[V_DISP] + q * dC[V_DISP];
        dE[P_DISP] =   C / b + d / b * dC[P_DISP] + B * dC[P_DISP] + C * dB[P_DISP] + dbeta[P_DISP] * a / b * C + beta * a / b * dC[P_DISP]  + dM[P_DISP] * c / b - dC[P_DISP] * c * c / b / b + dH[P_DISP] + q * dC[P_DISP];
        dE[R_GAS] =    d / b * dC[R_GAS] + B * dC[R_GAS] + C * dB[R_GAS] + dbeta[R_GAS] * a / b * C + beta * a / b * dC[R_GAS]  + dM[R_GAS] * c / b - dC[R_GAS] * c * c / b / b + dH[R_GAS] + q * dC[R_GAS];
	dE[V_GAS] =    d / b * dC[V_GAS] + B * dC[V_GAS] + C * dB[V_GAS] + dbeta[V_GAS] * a / b * C + beta * a / b * dC[V_GAS]  + dM[V_GAS] * c / b - dC[V_GAS] * c * c / b / b + dH[V_GAS] + q * dC[V_GAS];
	dE[P_GAS] =    d / b * dC[P_GAS] + B * dC[P_GAS] + C * dB[P_GAS] + dbeta[P_GAS] * a / b * C + beta * a / b * dC[P_GAS]  + dM[P_GAS] * c / b - dC[P_GAS] * c * c / b / b + dH[P_GAS] + q * dC[P_GAS];

        matrix[0][B_DISP] =   1.0 / sub_dt - dF[B_DISP] - C / b - a / b * dC[B_DISP];
        matrix[0][R_DISP] = - dF[R_DISP] + C * a / b / b - a / b * dC[R_DISP];
	matrix[0][V_DISP] = - dF[V_DISP] - a / b * dC[V_DISP];
	matrix[0][P_DISP] = - dF[P_DISP] - a / b * dC[P_DISP];
	matrix[0][R_GAS]  = - dF[R_GAS] - a / b * dC[R_GAS];
	matrix[0][V_GAS]  = - dF[V_GAS] - a / b * dC[V_GAS];
	matrix[0][P_GAS]  = - dF[P_GAS] - a / b * dC[P_GAS];

        matrix[1][B_DISP] = - dC[B_DISP];
	matrix[1][R_DISP] =   1.0 / sub_dt - dC[R_DISP];
	matrix[1][V_DISP] = - dC[V_DISP];
	matrix[1][P_DISP] = - dC[P_DISP];
	matrix[1][R_GAS]  = - dC[R_GAS];
	matrix[1][V_GAS]  = - dC[V_GAS];
	matrix[1][P_GAS]  = - dC[P_GAS];

	matrix[2][B_DISP] = - dM[B_DISP];
	matrix[2][R_DISP] = - dM[R_DISP];
	matrix[2][V_DISP] =   1.0 / sub_dt - dM[V_DISP];
	matrix[2][P_DISP] = - dM[P_DISP];
	matrix[2][R_GAS]  = - dM[R_GAS];
	matrix[2][V_GAS]  = - dM[V_GAS];
	matrix[2][P_GAS]  = - dM[P_GAS];

	matrix[3][B_DISP] = - dE_disp[B_DISP] + p * dF[B_DISP] + dp[B_DISP] * F + beta * dF[B_DISP] + F * dbeta[B_DISP];
	matrix[3][R_DISP] = - dE_disp[R_DISP] + p * dF[R_DISP] + dp[R_DISP] * F + beta * dF[R_DISP] + F * dbeta[R_DISP];
	matrix[3][V_DISP] = - dE_disp[V_DISP] + p * dF[V_DISP] + dp[V_DISP] * F + beta * dF[V_DISP] + F * dbeta[V_DISP];
	matrix[3][P_DISP] =   1.0 / sub_dt - dE_disp[P_DISP] + p * dF[P_DISP] + dp[P_DISP] * F + beta * dF[P_DISP] + F * dbeta[P_DISP];
	matrix[3][R_GAS]  = - dE_disp[R_GAS]  + p * dF[R_GAS]  + dp[R_GAS]  * F + beta * dF[R_GAS] + F * dbeta[R_GAS];
	matrix[3][V_GAS]  = - dE_disp[V_GAS]  + p * dF[V_GAS]  + dp[V_GAS]  * F + beta * dF[V_GAS] + F * dbeta[V_GAS];
	matrix[3][P_GAS]  = - dE_disp[P_GAS]  + p * dF[P_GAS]  + dp[P_GAS]  * F + beta * dF[P_GAS] + F * dbeta[P_GAS];

	matrix[4][B_DISP] =   dC[B_DISP];
	matrix[4][R_DISP] =   dC[R_DISP];
	matrix[4][V_DISP] =   dC[V_DISP];
	matrix[4][P_DISP] =   dC[P_DISP];
	matrix[4][R_GAS]  =   1.0 / sub_dt + dC[R_GAS];
	matrix[4][V_GAS]  =   dC[V_GAS];
	matrix[4][P_GAS]  =   dC[P_GAS];

	matrix[5][B_DISP] =   dM[B_DISP];
	matrix[5][R_DISP] =   dM[R_DISP];
	matrix[5][V_DISP] =   dM[V_DISP];
	matrix[5][P_DISP] =   dM[P_DISP];
	matrix[5][R_GAS]  =   dM[R_GAS];
	matrix[5][V_GAS]  =   1.0 / sub_dt + dM[V_GAS];
	matrix[5][P_GAS]  =   dM[P_GAS];

	matrix[6][B_DISP] =   dE[B_DISP] - p * dF[B_DISP] - dp[B_DISP] * F;
	matrix[6][R_DISP] =   dE[R_DISP] - p * dF[R_DISP] - dp[R_DISP] * F;
	matrix[6][V_DISP] =   dE[V_DISP] - p * dF[V_DISP] - dp[V_DISP] * F;
	matrix[6][P_DISP] =   dE[P_DISP] - p * dF[P_DISP] - dp[P_DISP] * F;
	matrix[6][R_GAS]  =   dE[R_GAS]  - p * dF[R_GAS]  - dp[R_GAS]  * F;
	matrix[6][V_GAS]  =   dE[V_GAS]  - p * dF[V_GAS]  - dp[V_GAS]  * F;
	matrix[6][P_GAS]  =   1.0 / sub_dt + dE[P_GAS] - p * dF[P_GAS] - dp[P_GAS] * F;

	fi[B_DISP] = (a - u_prev[0]) / sub_dt - F - C * a / b;
	fi[R_DISP] = (b - u_prev[1]) / sub_dt - C;
	fi[V_DISP] = (c - u_prev[2]) / sub_dt - momentum;
	fi[P_DISP] = (d - u_prev[3]) / sub_dt - (E_disp - p * F - beta * F);
	fi[R_GAS]  = (e - u_prev[4]) / sub_dt + C;
	fi[V_GAS]  = (f - u_prev[5]) / sub_dt + momentum;
	fi[P_GAS]  = (g - u_prev[6]) / sub_dt + (E - p * F);

	for (int i = 0; i < n; i++)
	    b_vector[i] = -fi[i];

	u[0] = a;
	u[1] = b;
	u[2] = c;
	u[3] = d;
	u[4] = e;
	u[5] = f;
	u[6] = g;

	for (int i = 0; i < n; i++)
		u_prev[i] = u[i];

	solve_linear_system(matrix, b_vector, n, paramsc->eps_general, du);

	for (int i = 0; i < n; i++)
		u[i] = u[i] + du[i];

	for (int i = 0; i < n; i++)
	    local_epsilon += (du[i] / u[i]) * (du[i] / u[i]);

	local_epsilon = sqrt(local_epsilon);

    }while( local_epsilon > 1.e-5);

    for (int i = 0; i < n; i++){
	v_cons[i] = u[i];
    }

}

// Решение полной системы уравнений (9) по неявной схеме Эйлера методом Ньютона - задача 2d2phc
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params2d - структура с дополнительными параметрами 2d2phc (in)
// v_cons[M] - текущий вектор консервативных переменных в ячейке (in / out)
// sub_dt - шаг по времени (in)
// number_of_block - номер блока (in)
// n - реальный размер векторов
void impicit_euler_full_2d( struct ParametersCommon *paramsc, struct Parameters2d *params2d, double v_cons[M], double sub_dt, int *number_of_block, int n ){

    double a, b, c, k, d, e, f, m, g; // Величины, отвечающие консервативным переменным (по порядку от B_DISP до P_GAS)
    double fi[M]; // вектор fi = 0
    double matrix[M][M]; // матрица частных производных

    double F, dF[M]; // Скорость коспактирования F, ее производная по всем переменным a ... g
    double momentum1, momentum2, dM1[M], dM2[M]; // Residual momentum exchange M, ее производная по всем переменным a ... g
    double E_disp, dE_disp[M]; // Residual energy exchange E дисперсной фазы , ее производная по всем переменным a ... g 
    double E, dE[M]; // Residual energy exchange E газовой фазы , ее производная по всем переменным a ... g
    double p, dp[M]; // давление газа, производная давления по всем переменным a ... g
    double beta, dbeta[M];
    double C, dC[M];
    double H, dH[M];
    double pressure_disp, dp_disp[M];
    double temp, temp_disp;
    double B, dB[M];

    double du[M], b_vector[M], u[M], u_prev[M]; // вектор решение du, ветор правой части b_vector, u - вектор переменных a ... g на текущем шаге, u_prev - вектор переменных a ... g на предыдущем шаге

    // Присваиваем начальные значения векторам u и u_prev
    for (int i = 0; i < n; i++){
	u[i] = v_cons[i];
	u_prev[i] = v_cons[i];
    }

    // Переприсваиваем значения параметров более удобным для использования переменным
    double muc = paramsc->compaction_viscosity; // взякость компактирования
    double gamma_disp = paramsc->g1; // показатель адиабаты дисперсной фазы (УРС)
    double gamma = paramsc->g2; // показатель адиабаты газовой фазы (УРС)
    double P0_disp = paramsc->p01; // параметр P0 дисперсной фазы (УРС)
    double p0 = params2d->block_values[number_of_block[0]][number_of_block[1]][P_GAS_2D]; // начальное значение давления газовой фазы
    double p0_disp = params2d->block_values[number_of_block[0]][number_of_block[1]][P_DISP_2D]; // начальное значение давления дисперсной фазы
    double r0_disp = params2d->block_values[number_of_block[0]][number_of_block[1]][R_DISP_2D]; // начальное значение плотности дисперсной фазы
    double b0_disp = params2d->block_values[number_of_block[0]][number_of_block[1]][B_DISP_2D]; // начальное значение объемной доли дисперсной фазы

    double A_constant = (p0 - p0_disp) / b0_disp / r0_disp * (2.0 - b0_disp) * (2.0 - b0_disp) / log(1.0 - b0_disp); // константа, присутствующая в определении скорости компактирования F, в частности в параметре бета
    double delta = paramsc->interface_drag_coef;
    double b_vir = paramsc->b_virial;
    double q = paramsc->heat_release;
    double sigma = paramsc->reaction_rate_prefactor;

    double local_epsilon;

    do{

	local_epsilon = 0.0;

	a = u[0];
	b = u[1];
	c = u[2];
	k = u[3];
	d = u[4];
	e = u[5];
	f = u[6];
        m = u[7];
        g = u[8];
			
        beta = -A_constant * b * log(1.0 - a) / (2.0 - a) / (2.0 - a);
        p = (g - (f * f + m * m) / 2.0 / e) * (gamma - 1.0) / (1.0 - a) * (1.0 + b_vir * e / (1.0 - a));
        pressure_disp = (d / b - 0.5 * (c * c + k * k) / b / b) * (gamma_disp - 1.0) * b / a - gamma_disp * P0_disp;
        F = a * (1.0 - a) / muc * (pressure_disp - p - beta);
        B = A_constant * log ((2.0 - b0_disp) / (2.0 - a) * pow( (1.0 - a), (1.0 - a) / (2.0 - a)) / pow((1.0 - b0_disp), (1.0 - b0_disp) / (2.0 - b0_disp ) ));

        if (!paramsc->is_compaction){
            B = 0.0;
            beta = 0.0;
            F = 0.0;
        }

        if (paramsc->is_burning && p > paramsc->p_ignition)
            C = -sigma * b * (p - paramsc->p_ignition);
        else
            C = 0.0;

        temp = T_gas(paramsc, p, e / (1.0 - a));
        temp_disp = T_disp(paramsc, pressure_disp, b / a);
        H = paramsc->heat_transfer_coefficient * (temp - temp_disp);
        if (!paramsc->is_heat_transfer)
            H = 0.0;

        momentum1 = (delta + 0.5 * C) * (f/e - c/b) + C * c / b;
        momentum2 = (delta + 0.5 * C) * (m/e - k/b) + C * k / b;
        if (!paramsc->is_friction_force){
            momentum1 = 0.0;
            momentum2 = 0.0;
        }
	E_disp = (momentum1 - C * c / b) * c / b + (momentum2 - C * k / b) * k / b + d / b * C + H;
        E = E_disp + (B + beta * a / b + q) * C;

        dp[B_DISP_2D] =   (g - (f * f + m * m) / 2.0 / e) * (gamma - 1.0) / (1.0 - a) / (1.0 - a) * (1.0 + 2 * b_vir * e / (1.0 - a));
	dp[R_DISP_2D] =   0.0;
	dp[V_DISP_2D] =   0.0;
        dp[U_DISP_2D] =   0.0;
	dp[P_DISP_2D] =   0.0;
	dp[R_GAS_2D]  =   (f * f + m * m) / 2.0 / e / e * (gamma - 1.0) / (1.0 - a) * (1.0 + b_vir * e / (1.0 - a)) + 
            (g - 0.5 * (f * f + m * m) / e ) * (gamma - 1.0) * b_vir / (1.0 - a) / (1.0 - a);
	dp[V_GAS_2D]  = - f / e * (gamma - 1.0) / (1.0 - a) * (1.0 + 1.0 * b_vir * e / (1.0 - a));
        dp[U_GAS_2D]  = - m / e * (gamma - 1.0) / (1.0 - a) * (1.0 + 1.0 * b_vir * e / (1.0 - a));
	dp[P_GAS_2D]  =   (gamma - 1.0) / (1.0 - a) * (1.0 + 1.0 * b_vir * e / (1.0 - a));

        dbeta[B_DISP_2D] = A_constant * b /  (2.0 - a) / (2.0 - a) * ( 1.0 / (a - 1.0) + 2.0 * log(1.0 - a) / (2.0 - a));
	dbeta[R_DISP_2D] = A_constant * log(1.0 - a) / (2.0 - a) / (2.0 - a);
	dbeta[V_DISP_2D] = 0.0;
        dbeta[U_DISP_2D] = 0.0;
	dbeta[P_DISP_2D] = 0.0;
	dbeta[R_GAS_2D] = 0.0;
	dbeta[V_GAS_2D] = 0.0;
        dbeta[U_GAS_2D] = 0.0;
        dbeta[P_GAS_2D] = 0.0;

        dp_disp[B_DISP_2D] = - (d - 0.5 * (c * c + k * k) / b ) * (gamma_disp - 1.0) / a / a;
        dp_disp[R_DISP_2D] = 0.5 * (c * c + k * k) / b / b * (gamma_disp - 1.0) / a;
        dp_disp[V_DISP_2D] = - c / b * (gamma_disp - 1.0) / a;
        dp_disp[U_DISP_2D] = - k / b * (gamma_disp - 1.0) / a;
        dp_disp[P_DISP_2D] = (gamma_disp - 1.0) / a;
        dp_disp[R_GAS_2D]  = 0.0;
        dp_disp[V_GAS_2D]  = 0.0;
        dp_disp[U_GAS_2D]  = 0.0;
        dp_disp[P_GAS_2D]  = 0.0;

        dF[B_DISP_2D] = (1.0 - 2.0 * a) / muc * (pressure_disp - p - beta) + a * (1.0 - a) / muc * (dp_disp[B_DISP_2D] - dp[B_DISP_2D] - dbeta[B_DISP_2D]);
	dF[R_DISP_2D] = a * (1.0 - a) / muc * (dp_disp[R_DISP_2D] - dp[R_DISP_2D] - dbeta[R_DISP_2D]);
	dF[V_DISP_2D] = a * (1.0 - a) / muc * (dp_disp[V_DISP_2D] - dp[V_DISP_2D] - dbeta[V_DISP_2D]);
	dF[U_DISP_2D] = a * (1.0 - a) / muc * (dp_disp[U_DISP_2D] - dp[U_DISP_2D] - dbeta[U_DISP_2D]);
	dF[P_DISP_2D] = a * (1.0 - a) / muc * (dp_disp[P_DISP_2D] - dp[P_DISP_2D] - dbeta[P_DISP_2D]);
	dF[R_GAS_2D]  = a * (1.0 - a) / muc * (dp_disp[R_GAS_2D] - dp[R_GAS_2D] - dbeta[R_GAS_2D]);
	dF[V_GAS_2D]  = a * (1.0 - a) / muc * (dp_disp[V_GAS_2D] - dp[V_GAS_2D] - dbeta[V_GAS_2D]);
        dF[U_GAS_2D]  = a * (1.0 - a) / muc * (dp_disp[U_GAS_2D] - dp[U_GAS_2D] - dbeta[U_GAS_2D]);
	dF[P_GAS_2D]  = a * (1.0 - a) / muc * (dp_disp[P_GAS_2D] - dp[P_GAS_2D] - dbeta[P_GAS_2D]);

        dB[B_DISP_2D] = beta / b;
        dB[R_DISP_2D] = 0.0;
	dB[V_DISP_2D] = 0.0;
        dB[U_DISP_2D] = 0.0;
	dB[P_DISP_2D] = 0.0;
	dB[R_GAS_2D] = 0.0;
	dB[V_GAS_2D] = 0.0;
        dB[U_GAS_2D] = 0.0;
        dB[P_GAS_2D] = 0.0;

        dC[B_DISP_2D] = - sigma * b * dp[B_DISP_2D];
        dC[R_DISP_2D] = - sigma * b * dp[R_DISP_2D] - sigma * (p - paramsc->p_ignition);
	dC[V_DISP_2D] = - sigma * b * dp[V_DISP_2D];
	dC[U_DISP_2D] = - sigma * b * dp[U_DISP_2D];
	dC[P_DISP_2D] = - sigma * b * dp[P_DISP_2D];
	dC[R_GAS_2D]  = - sigma * b * dp[R_GAS_2D];
	dC[V_GAS_2D]  = - sigma * b * dp[V_GAS_2D];
        dC[U_GAS_2D]  = - sigma * b * dp[U_GAS_2D];
	dC[P_GAS_2D]  = - sigma * b * dp[P_GAS_2D];
        /*
        dH[B_DISP_2D] =   paramsc->heat_transfer_coefficient * P0_disp / paramsc->specific_heat_cv1 / b;
        dH[R_DISP_2D] =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv1 / b / b * (d - (c * c + k * k) / b - P0_disp * a);
	dH[V_DISP_2D] =   paramsc->heat_transfer_coefficient * c / paramsc->specific_heat_cv1 / b / b;
	dH[U_DISP_2D] =   paramsc->heat_transfer_coefficient * k / paramsc->specific_heat_cv1 / b / b;
        dH[P_DISP_2D] = - paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv1 / b;
	dH[R_GAS_2D]  = - paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e / e * (g - (f * f + m * m) / e);
	dH[V_GAS_2D]  = - paramsc->heat_transfer_coefficient * f / paramsc->specific_heat_cv2 / e / e;
	dH[U_GAS_2D]  = - paramsc->heat_transfer_coefficient * m / paramsc->specific_heat_cv2 / e / e;
        dH[P_GAS_2D]  =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e;
        */
        dH[B_DISP_2D] = - paramsc->heat_transfer_coefficient * (pressure_disp + P0_disp) / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b 
                        - paramsc->heat_transfer_coefficient * dp_disp[B_DISP_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

        dH[R_DISP_2D] =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv1 / b / b * (pressure_disp + P0_disp) * a
                        - paramsc->heat_transfer_coefficient * dp_disp[R_DISP_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

	dH[V_DISP_2D] = - paramsc->heat_transfer_coefficient * dp_disp[V_DISP_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

	dH[U_DISP_2D] = - paramsc->heat_transfer_coefficient * dp_disp[U_DISP_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

        dH[P_DISP_2D] = - paramsc->heat_transfer_coefficient * dp_disp[P_DISP_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

	dH[R_GAS_2D]  = - paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e / e * (g - (f * f + m * m) / e)
                        - paramsc->heat_transfer_coefficient * dp_disp[R_GAS_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

	dH[V_GAS_2D]  = - paramsc->heat_transfer_coefficient * f / paramsc->specific_heat_cv2 / e / e
                        - paramsc->heat_transfer_coefficient * dp_disp[V_GAS_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

	dH[U_GAS_2D]  = - paramsc->heat_transfer_coefficient * m / paramsc->specific_heat_cv2 / e / e
                        - paramsc->heat_transfer_coefficient * dp_disp[U_GAS_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

        dH[P_GAS_2D]  =   paramsc->heat_transfer_coefficient / paramsc->specific_heat_cv2 / e
                        - paramsc->heat_transfer_coefficient * dp_disp[P_GAS_2D] * a / (gamma_disp - 1.0) / paramsc->specific_heat_cv1 / b ;

        for (int i = 0; i < n; i++){
            if (!paramsc->is_compaction){
                dbeta[i] = 0.0;
                dB[i] = 0.0;
                dF[i] = 0.0;
            }
            if (!paramsc->is_burning || C == 0.0 ){
                dC[i] = 0.0;
            }
            if (!paramsc->is_heat_transfer)
                dH[i] = 0.0;
        }
  
        dM1[B_DISP_2D] =   dC[B_DISP_2D] * c / b + 0.5 * dC[B_DISP_2D] * (f / e - c / b);
        dM1[R_DISP_2D] =   dC[R_DISP_2D] * c / b - C * c / b / b + delta * c / b / b + 0.5 * dC[R_DISP_2D] * (f / e - c / b) + 0.5 * C * c / b / b;
        dM1[V_DISP_2D] =   C / b + dC[V_DISP_2D] * c / b - delta / b + 0.5 * dC[V_DISP_2D] * (f / e - c / b) - 0.5 * C / b;
        dM1[U_DISP_2D] =   dC[U_DISP_2D] * c / b + 0.5 * dC[U_DISP_2D] * (f / e - c / b);
        dM1[P_DISP_2D] =   dC[P_DISP_2D] * c / b + 0.5 * dC[P_DISP_2D] * (f / e - c / b);
        dM1[R_GAS_2D]  =   dC[R_GAS_2D] * c / b - delta * f / e / e + 0.5 * dC[R_GAS_2D] * f / e - 0.5 * C * f / e / e - 0.5 * dC[R_GAS_2D] * c / b;
        dM1[V_GAS_2D]  =   dC[V_GAS_2D] * c / b + delta / e + 0.5 * dC[V_GAS_2D] * f / e + 0.5 * C / e - 0.5 * dC[V_GAS_2D] * c / b;
	dM1[U_GAS_2D]  =   dC[U_GAS_2D] * c / b + 0.5 * dC[U_GAS_2D] * (f / e - c / b);
        dM1[P_GAS_2D]  =   dC[P_GAS_2D] * c / b + 0.5 * dC[P_GAS_2D] * (f / e - c / b);

        dM2[B_DISP_2D] =   dC[B_DISP_2D] * k / b + 0.5 * dC[B_DISP_2D] * (m / e - k / b);
        dM2[R_DISP_2D] =   dC[R_DISP_2D] * k / b - C * k / b / b + delta * k / b / b + 0.5 * dC[R_DISP_2D] * (m / e - k / b) + 0.5 * C * k / b / b;
        dM2[V_DISP_2D] =   dC[V_DISP_2D] * k / b + 0.5 * dC[V_DISP_2D] * (m / e - k / b);
        dM2[U_DISP_2D] =   C / b + dC[U_DISP_2D] * k / b - delta / b + 0.5 * dC[U_DISP_2D] * (m / e - k / b) - 0.5 * C / b;
        dM2[P_DISP_2D] =   dC[P_DISP_2D] * k / b + 0.5 * dC[P_DISP_2D] * (m / e - k / b);
        dM2[R_GAS_2D]  =   dC[R_GAS_2D] * k / b - delta * m / e / e + 0.5 * dC[R_GAS_2D] * m / e - 0.5 * C * m / e / e - 0.5 * dC[R_GAS_2D] * k / b;
        dM2[V_GAS_2D]  =   dC[V_GAS_2D] * k / b + 0.5 * dC[V_GAS_2D] * (m / e - k / b);
	dM2[U_GAS_2D]  =   dC[U_GAS_2D] * k / b + 0.5 * dC[U_GAS_2D] * (m / e - k / b) + 1.0 / e * (delta + 0.5 * C);
        dM2[P_GAS_2D]  =   dC[P_GAS_2D] * k / b + 0.5 * dC[P_GAS_2D] * (m / e - k / b);

        for (int i = 0; i < n; i++){
            if (!paramsc->is_friction_force){
                dM1[i] = 0.0;
                dM2[i] = 0.0;
            }
        }

        dE_disp[B_DISP_2D] =   d / b * dC[B_DISP_2D] + (dM1[B_DISP_2D] - dC[B_DISP_2D] * c / b ) * c / b + (dM2[B_DISP_2D] - dC[B_DISP_2D] * k / b ) * k / b + 
            dH[B_DISP_2D];
        dE_disp[R_DISP_2D] = - d * C / b / b + d / b * dC[R_DISP_2D] + (dM1[R_DISP_2D] - dC[R_DISP_2D] * c / b + C * c / b / b) * c / b - ( momentum1 - C * c / b) * c / b / b +
            (dM2[R_DISP_2D] - dC[R_DISP_2D] * k / b + C * k / b / b) * k / b - ( momentum2 - C * k / b) * k / b / b + dH[R_DISP_2D];
        dE_disp[V_DISP_2D] =   d / b * dC[V_DISP_2D] + (dM1[V_DISP_2D] - dC[V_DISP_2D] * c / b - C / b) * c / b + (momentum1 - C * c / b) / b + (dM2[V_DISP_2D] - dC[V_DISP_2D] * k / b) * k / b + dH[V_DISP_2D];
        dE_disp[U_DISP_2D] =   d / b * dC[U_DISP_2D] + (dM1[U_DISP_2D] - dC[U_DISP_2D] * c / b) * c / b + (momentum2 - C * k / b) / b + (dM2[U_DISP_2D] - dC[U_DISP_2D] * k / b - C / b) * k / b + dH[U_DISP_2D];
        dE_disp[P_DISP_2D] =   C / b + d / b * dC[P_DISP_2D] + (dM1[P_DISP_2D] - dC[P_DISP_2D] * c / b) * c / b + (dM2[P_DISP_2D] - dC[P_DISP_2D] * k / b) * k / b + dH[P_DISP_2D];
        dE_disp[R_GAS_2D] =    d / b * dC[R_GAS_2D] + (dM1[R_GAS_2D] - dC[R_GAS_2D] * c / b) * c / b + (dM2[R_GAS_2D] - dC[R_GAS_2D] * k / b) * k / b + dH[R_GAS_2D];
	dE_disp[V_GAS_2D] =    d / b * dC[V_GAS_2D] + (dM1[V_GAS_2D] - dC[V_GAS_2D] * c / b) * c / b + (dM2[V_GAS_2D] - dC[V_GAS_2D] * k / b) * k / b + dH[V_GAS_2D];
	dE_disp[U_GAS_2D] =    d / b * dC[U_GAS_2D] + (dM1[U_GAS_2D] - dC[U_GAS_2D] * c / b) * c / b + (dM2[U_GAS_2D] - dC[U_GAS_2D] * k / b) * k / b + dH[U_GAS_2D];
        dE_disp[P_GAS_2D] =    d / b * dC[P_GAS_2D] + (dM1[P_GAS_2D] - dC[P_GAS_2D] * c / b) * c / b + (dM2[P_GAS_2D] - dC[P_GAS_2D] * k / b) * k / b + dH[P_GAS_2D];

        dE[B_DISP_2D] = dE_disp[B_DISP_2D] + (dB[B_DISP_2D] + dbeta[B_DISP_2D] * a / b + beta / b) * C + (B + q + beta * a / b) * dC[B_DISP_2D];
        dE[R_DISP_2D] = dE_disp[R_DISP_2D] + (dB[R_DISP_2D] + dbeta[R_DISP_2D] * a / b - beta * a / b / b) * C + (B + q + beta * a / b) * dC[R_DISP_2D];
        dE[V_DISP_2D] = dE_disp[V_DISP_2D] + (dB[V_DISP_2D] + dbeta[V_DISP_2D] * a / b) * C + (B + q + beta * a / b) * dC[V_DISP_2D];
        dE[U_DISP_2D] = dE_disp[U_DISP_2D] + (dB[U_DISP_2D] + dbeta[U_DISP_2D] * a / b) * C + (B + q + beta * a / b) * dC[U_DISP_2D];
        dE[P_DISP_2D] = dE_disp[P_DISP_2D] + (dB[P_DISP_2D] + dbeta[P_DISP_2D] * a / b) * C + (B + q + beta * a / b) * dC[P_DISP_2D];
        dE[R_GAS_2D] = dE_disp[R_GAS_2D] + (dB[R_GAS_2D] + dbeta[R_GAS_2D] * a / b) * C + (B + q + beta * a / b) * dC[R_GAS_2D];
        dE[V_GAS_2D] = dE_disp[V_GAS_2D] + (dB[V_GAS_2D] + dbeta[V_GAS_2D] * a / b) * C + (B + q + beta * a / b) * dC[V_GAS_2D];
        dE[U_GAS_2D] = dE_disp[U_GAS_2D] + (dB[U_GAS_2D] + dbeta[U_GAS_2D] * a / b) * C + (B + q + beta * a / b) * dC[U_GAS_2D];
        dE[P_GAS_2D] = dE_disp[P_GAS_2D] + (dB[P_GAS_2D] + dbeta[P_GAS_2D] * a / b) * C + (B + q + beta * a / b) * dC[P_GAS_2D];

        matrix[0][B_DISP_2D] =   1.0 / sub_dt - dF[B_DISP_2D] - C / b - a / b * dC[B_DISP_2D];
        matrix[0][R_DISP_2D] = - dF[R_DISP_2D] + C * a / b / b - a / b * dC[R_DISP_2D];
	matrix[0][V_DISP_2D] = - dF[V_DISP_2D] - a / b * dC[V_DISP_2D];
        matrix[0][U_DISP_2D] = - dF[U_DISP_2D] - a / b * dC[U_DISP_2D];
	matrix[0][P_DISP_2D] = - dF[P_DISP_2D] - a / b * dC[P_DISP_2D];
	matrix[0][R_GAS_2D]  = - dF[R_GAS_2D] - a / b * dC[R_GAS_2D];
	matrix[0][V_GAS_2D]  = - dF[V_GAS_2D] - a / b * dC[V_GAS_2D];
        matrix[0][U_GAS_2D]  = - dF[U_GAS_2D] - a / b * dC[U_GAS_2D];
	matrix[0][P_GAS_2D]  = - dF[P_GAS_2D] - a / b * dC[P_GAS_2D];

        matrix[1][B_DISP_2D] = - dC[B_DISP_2D];
	matrix[1][R_DISP_2D] =   1.0 / sub_dt - dC[R_DISP_2D];
	matrix[1][V_DISP_2D] = - dC[V_DISP_2D];
        matrix[1][U_DISP_2D] = - dC[U_DISP_2D];
	matrix[1][P_DISP_2D] = - dC[P_DISP_2D];
	matrix[1][R_GAS_2D]  = - dC[R_GAS_2D];
	matrix[1][V_GAS_2D]  = - dC[V_GAS_2D];
        matrix[1][U_GAS_2D]  = - dC[U_GAS_2D];
	matrix[1][P_GAS_2D]  = - dC[P_GAS_2D];

	matrix[2][B_DISP_2D] = - dM1[B_DISP_2D];
	matrix[2][R_DISP_2D] = - dM1[R_DISP_2D];
	matrix[2][V_DISP_2D] =   1.0 / sub_dt - dM1[V_DISP_2D];
        matrix[2][U_DISP_2D] = - dM1[U_DISP_2D];
	matrix[2][P_DISP_2D] = - dM1[P_DISP_2D];
	matrix[2][R_GAS_2D]  = - dM1[R_GAS_2D];
	matrix[2][V_GAS_2D]  = - dM1[V_GAS_2D];
        matrix[2][U_GAS_2D]  = - dM1[U_GAS_2D];
	matrix[2][P_GAS_2D]  = - dM1[P_GAS_2D];

        matrix[3][B_DISP_2D] = - dM2[B_DISP_2D];
	matrix[3][R_DISP_2D] = - dM2[R_DISP_2D];
	matrix[3][V_DISP_2D] = - dM2[V_DISP_2D];
        matrix[3][U_DISP_2D] =   1.0 / sub_dt- dM2[U_DISP_2D];
	matrix[3][P_DISP_2D] = - dM2[P_DISP_2D];
	matrix[3][R_GAS_2D]  = - dM2[R_GAS_2D];
	matrix[3][V_GAS_2D]  = - dM2[V_GAS_2D];
        matrix[3][U_GAS_2D]  = - dM2[U_GAS_2D];
	matrix[3][P_GAS_2D]  = - dM2[P_GAS_2D];

        matrix[4][B_DISP_2D] =  (dp[B_DISP_2D] + dbeta[B_DISP_2D]) * F + (p + beta) * dF[B_DISP_2D] - dE_disp[B_DISP_2D];
	matrix[4][R_DISP_2D] =  (dp[R_DISP_2D] + dbeta[R_DISP_2D]) * F + (p + beta) * dF[R_DISP_2D] - dE_disp[R_DISP_2D];
	matrix[4][V_DISP_2D] =  (dp[V_DISP_2D] + dbeta[V_DISP_2D]) * F + (p + beta) * dF[V_DISP_2D] - dE_disp[V_DISP_2D];
	matrix[4][U_DISP_2D] =  (dp[U_DISP_2D] + dbeta[U_DISP_2D]) * F + (p + beta) * dF[U_DISP_2D] - dE_disp[U_DISP_2D];
        matrix[4][P_DISP_2D] =   1.0 / sub_dt + (dp[P_DISP_2D] + dbeta[P_DISP_2D]) * F + (p + beta) * dF[P_DISP_2D] - dE_disp[P_DISP_2D];
	matrix[4][R_GAS_2D]  =  (dp[R_GAS_2D] + dbeta[R_GAS_2D]) * F + (p + beta) * dF[R_GAS_2D] - dE_disp[R_GAS_2D];
	matrix[4][V_GAS_2D]  =  (dp[V_GAS_2D] + dbeta[V_GAS_2D]) * F + (p + beta) * dF[V_GAS_2D] - dE_disp[V_GAS_2D];
	matrix[4][U_GAS_2D]  =  (dp[U_GAS_2D] + dbeta[U_GAS_2D]) * F + (p + beta) * dF[U_GAS_2D] - dE_disp[U_GAS_2D];
        matrix[4][P_GAS_2D]  =  (dp[P_GAS_2D] + dbeta[P_GAS_2D]) * F + (p + beta) * dF[P_GAS_2D] - dE_disp[P_GAS_2D];

	matrix[5][B_DISP_2D] =   dC[B_DISP_2D];
	matrix[5][R_DISP_2D] =   dC[R_DISP_2D];
	matrix[5][V_DISP_2D] =   dC[V_DISP_2D];
        matrix[5][U_DISP_2D] =   dC[U_DISP_2D];
	matrix[5][P_DISP_2D] =   dC[P_DISP_2D];
	matrix[5][R_GAS_2D]  =   1.0 / sub_dt + dC[R_GAS_2D];
	matrix[5][V_GAS_2D]  =   dC[V_GAS_2D];
        matrix[5][U_GAS_2D] =    dC[U_GAS_2D];
	matrix[5][P_GAS_2D]  =   dC[P_GAS_2D];

	matrix[6][B_DISP_2D] =   dM1[B_DISP_2D];
	matrix[6][R_DISP_2D] =   dM1[R_DISP_2D];
	matrix[6][V_DISP_2D] =   dM1[V_DISP_2D];
        matrix[6][U_DISP_2D] =   dM1[U_DISP_2D];
	matrix[6][P_DISP_2D] =   dM1[P_DISP_2D];
	matrix[6][R_GAS_2D]  =   dM1[R_GAS_2D];
	matrix[6][V_GAS_2D]  =   1.0 / sub_dt + dM1[V_GAS_2D];
        matrix[6][U_GAS_2D] =    dM1[U_GAS_2D];
	matrix[6][P_GAS_2D]  =   dM1[P_GAS_2D];

        matrix[7][B_DISP_2D] =   dM2[B_DISP_2D];
	matrix[7][R_DISP_2D] =   dM2[R_DISP_2D];
	matrix[7][V_DISP_2D] =   dM2[V_DISP_2D];
        matrix[7][U_DISP_2D] =   dM2[U_DISP_2D];
	matrix[7][P_DISP_2D] =   dM2[P_DISP_2D];
	matrix[7][R_GAS_2D]  =   dM2[R_GAS_2D];
	matrix[7][V_GAS_2D]  =   dM2[V_GAS_2D];
        matrix[7][U_GAS_2D] =    1.0 / sub_dt + dM2[U_GAS_2D];
	matrix[7][P_GAS_2D]  =   dM2[P_GAS_2D];

	matrix[8][B_DISP_2D] =   dE[B_DISP_2D] - p * dF[B_DISP_2D] - dp[B_DISP_2D] * F;
	matrix[8][R_DISP_2D] =   dE[R_DISP_2D] - p * dF[R_DISP_2D] - dp[R_DISP_2D] * F;
	matrix[8][V_DISP_2D] =   dE[V_DISP_2D] - p * dF[V_DISP_2D] - dp[V_DISP_2D] * F;
        matrix[8][U_DISP_2D] =   dE[U_DISP_2D] - p * dF[U_DISP_2D] - dp[U_DISP_2D] * F;
	matrix[8][P_DISP_2D] =   dE[P_DISP_2D] - p * dF[P_DISP_2D] - dp[P_DISP_2D] * F;
	matrix[8][R_GAS_2D]  =   dE[R_GAS_2D]  - p * dF[R_GAS_2D]  - dp[R_GAS_2D]  * F;
	matrix[8][V_GAS_2D]  =   dE[V_GAS_2D]  - p * dF[V_GAS_2D]  - dp[V_GAS_2D]  * F;
        matrix[8][U_GAS_2D] =    dE[U_GAS_2D]  - p * dF[U_GAS_2D]  - dp[U_GAS_2D] * F;
	matrix[8][P_GAS_2D]  =   1.0 / sub_dt + dE[P_GAS_2D] - p * dF[P_GAS_2D] - dp[P_GAS_2D] * F;

	fi[B_DISP_2D] = (a - u_prev[0]) / sub_dt - F - C * a / b;
	fi[R_DISP_2D] = (b - u_prev[1]) / sub_dt - C;
	fi[V_DISP_2D] = (c - u_prev[2]) / sub_dt - momentum1;
        fi[U_DISP_2D] = (k - u_prev[3]) / sub_dt - momentum2;
	fi[P_DISP_2D] = (d - u_prev[4]) / sub_dt - (E_disp - p * F - beta * F);
	fi[R_GAS_2D]  = (e - u_prev[5]) / sub_dt + C;
	fi[V_GAS_2D]  = (f - u_prev[6]) / sub_dt + momentum1;
        fi[U_GAS_2D]  = (m - u_prev[7]) / sub_dt + momentum2;
	fi[P_GAS_2D]  = (g - u_prev[8]) / sub_dt + (E - p * F);

	for (int i = 0; i < n; i++)
	    b_vector[i] = -fi[i];

	u[0] = a;
	u[1] = b;
	u[2] = c;
        u[3] = k;
	u[4] = d;
	u[5] = e;
	u[6] = f;
        u[7] = m;
	u[8] = g;

	for (int i = 0; i < n; i++)
	    u_prev[i] = u[i];

	solve_linear_system(matrix, b_vector, n, paramsc->eps_general, du);

	for (int i = 0; i < n; i++)
	    u[i] = u[i] + du[i];

	for (int i = 0; i < n; i++){

            if ( fabs(u[i] ) > 1.e-14 )
	        local_epsilon += (du[i] / u[i]) * (du[i] / u[i]);

        }
	local_epsilon = sqrt(local_epsilon);
    //    printf("local epsilon = %.25lf\n", local_epsilon);
    }while( local_epsilon > 1.e-5);
    
    for (int i = 0; i < n; i++){
	v_cons[i] = u[i];
    }
}