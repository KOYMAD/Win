// hll.cc
// Метод HLL из Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) Уткин Павел, 2015
// Создан: 8 июля 2015 г.

#include "hll.h"
#include "Moving\Constants.h"
#include "godunov_bn.h"
// Метод HLL численного интегрирования уравнений Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией (in)
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге
// step_number - номер текущего шага по пространству
// n - реальный размер векторов
// is_pressure_relaxation_after_this_step - true, если после решения одномерной задачи по данному направлению нужно применять релаксацию давлений; false - иначе
//                              если в данной задачу вообще не нужно применять релаксацию давлений, то true не повлияет на результат, так как данная опция контролируется в parameters.dat
// number_of_scalars - число лагранжевых скаляров
// curr_time - текущий момент времени
// configuration_presure - конфигурационное давление
void hll( struct ParametersCommon* paramsc, struct Parameters1d* params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
          double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure, double body_velocity, int i, int *status, double cont_left[M], double cont_right[M]  ) {
    double solution_ncons_Lh[M]; // вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
    if (status[i+1] == BOUNDARY)
        printf("\n before hll %lf",center_ncons[P_GAS]);
          Lh_HLL(debug_info, params1d, paramsc, left_ncons, center_ncons, right_ncons, slopes_left, slopes_center, slopes_right,
          dt, h, solution_ncons_Lh, n, number_of_scalars, i, status, body_velocity,  cont_left, cont_right ); 
              if (status[i+1] == BOUNDARY)
        printf("\n after hll %lf",solution_ncons_Lh[P_GAS]);// действие гиперболического оператора
          //printf("\n   hll %lf %lf %lf %lf %d %d", center_ncons[P_GAS], center_ncons[P_DISP], solution_ncons_Lh[P_GAS], solution_ncons_Lh[P_DISP], i, status[i]);
    if ( paramsc->pressure_relaxation == true && is_pressure_relaxation_after_this_step == true ){
        if (status[i] == INNER){

        Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // действие релаксационного оператора
        if (status[i+1] == BOUNDARY)
        printf("\n after rel %lf",solution_ncons[P_GAS]);
        }
        else{}

        if ( status[i] == BOUNDARY){
            if (status[i+1] == GHOST){
                //printf("\n %lf",solution_ncons[P_GAS]);
                //double right_ghost[M];
                //for (int l = 0; l < n; l++)
                //    right_ghost[l] = center_ncons[l];
                //right_ghost[V_GAS] = 2 * body_velocity -  right_ghost[V_GAS];
                //right_ghost[V_DISP] = 2 * body_velocity - right_ghost[V_DISP];
                Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh,cont_left , dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // действие релаксационного оператора
            }
        
            else {
                //printf("\n %lf",solution_ncons[P_GAS]);
                //double left_ghost[M];
                //for (int l = 0; l < n; l++)
                //    left_ghost[l] = center_ncons[l];
                //left_ghost[V_GAS] = 2 * body_velocity - left_ghost[V_GAS];
                //left_ghost[V_DISP] = 2 * body_velocity - left_ghost[V_DISP];
                Lr( paramsc, params1d, debug_info, cont_right, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // действие релаксационного оператора
            }
        }
        else{}
    if (status[i] == GHOST || status[i] == OUTER){
        
       for ( int j = 0; j < n; j++ )
            solution_ncons[j] = solution_ncons_Lh[j];
       //printf(" g 0");
        }
    
    else{}
    }
    else {
        for ( int j = 0; j < n; j++ )
            solution_ncons[j] = solution_ncons_Lh[j];
    }
                if (solution_ncons[B_DISP] > 0.9)
                solution_ncons[B_DISP] = 0.9;
    //printf(" rel %lf", solution_ncons[P_GAS]);

}

// Гиперболический оператор метода HLL
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
// n - реальный размер векторов без учета лагранжевых скаляров
// number_of_scalars - количество дополнительных уравнений, оно же количество лагранжевых скаляров
// cont_ncons - вектор переменных на контактном разрыве для потока через граничную ячейку подвижного тела
void Lh_HLL(struct DebugInfo *debug_info, struct Parameters1d* params1d,  struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars,  int i, int *status, double body_velocity, double cont_left[M], double cont_right[M] ) {
                 //printf("\n center %lf %lf %lf", center_ncons[P_GAS], left_ncons[P_GAS], right_ncons[P_GAS]);
    // реконструкция сеточных функций
    double left_minus_ncons[M];
    double left_minus_cons[M];
    double left_plus_ncons[M];
    double left_plus_cons[M];
    double right_minus_ncons[M];
    double right_minus_cons[M];
    double right_plus_ncons[M];
    double right_plus_cons[M];
    
    for ( int j = 0; j < n; j++ ) {
        // реконструкция значений на левой границе слева
        convert_noncons_to_cons( paramsc, left_ncons, left_minus_cons, number_of_scalars );
        left_minus_cons[j] += 0.5 * h * slopes_left[j];
        convert_cons_to_noncons( paramsc, left_minus_cons, left_minus_ncons, number_of_scalars );
        // реконструкция значений на левой границе справа
        convert_noncons_to_cons( paramsc, center_ncons, left_plus_cons, number_of_scalars );
        left_plus_cons[j] -= 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, left_plus_cons, left_plus_ncons, number_of_scalars );
        // реконструкция значений на правой границе слева
        convert_noncons_to_cons( paramsc, center_ncons, right_minus_cons, number_of_scalars );
        right_minus_cons[j] += 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, right_minus_cons, right_minus_ncons, number_of_scalars );
        // реконструкция значений на правой границе справа
        convert_noncons_to_cons( paramsc, right_ncons, right_plus_cons, number_of_scalars );
        right_plus_cons[j] -= 0.5 * h * slopes_right[j];
        convert_cons_to_noncons( paramsc, right_plus_cons, right_plus_ncons, number_of_scalars ); 
        
    } 
    debug_info->current_cell = i;
    debug_info->neighbour_cell = debug_info->current_cell - 1;
    for ( int j = 0; j < n; j++ )
        debug_info->current_cell_vncons[j] = left_plus_ncons[j];
    for ( int j = 0; j < n; j++ )
        debug_info->neighbour_cell_vncons[j] = left_minus_ncons[j];
    double splus_l, sminus_l;
    double splus_r, sminus_r;
    array1D left_flux( M );
    array1D right_flux( M );
    if (status[i] == GHOST){
        if ( BOUNDARY == status[i-1]){
            // расчет потока через левую границу ячейки

            for (int j = 0;  j < n; j++)
                left_plus_ncons[j] = left_minus_ncons[j];

            left_plus_ncons[V_GAS] = 2 * body_velocity - left_plus_ncons[V_GAS];
            left_plus_ncons[V_DISP] = 2 * body_velocity - left_plus_ncons[V_DISP];
            hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars );

            // расчет потока через правую границу ячейки
            hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
            for (int j = 0; j < n; j++)
                right_flux[j] = 0;

        }
        else{}
        
        if ( BOUNDARY == status[i + 1]){            
            hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars );
            for (int j = 0; j < n; j++)
                left_flux[j] = 0;
            // расчет потока через правую границу ячейки
            for (int j = 0;  j < n; j++)
                right_minus_ncons[j] = right_plus_ncons[j];
            right_minus_ncons[V_GAS] = 2 * body_velocity - right_minus_ncons[V_GAS];
            right_minus_ncons[V_DISP] = 2 * body_velocity - right_minus_ncons[V_DISP];

            hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
        }
        else{}
    }
    
 
    if (status[i] == OUTER){
        for (int j = 0; j < n; j++){
            right_flux[j] = 0;
            left_flux[j] = 0;
        }
    }
     
    if ( status[i] == BOUNDARY){
        if ( GHOST == status[i+1]){

            for (int j = 0;  j < n; j++)
                right_plus_ncons[j] = right_minus_ncons[j];
            //printf("\n before %lf", right_minus_ncons[P_GAS]);
            right_plus_ncons[V_GAS] = 2 * body_velocity - right_plus_ncons[V_GAS]; 
            right_plus_ncons[V_DISP] = 2 * body_velocity - right_plus_ncons[V_DISP];
            //printf("\n rg %lf %lf", right_plus_ncons[V_GAS], right_plus_ncons[V_DISP]);
            //convert_noncons_to_cons(paramsc,right_plus_ncons,right_plus_cons, number_of_scalars);
            double cont_ncons[M];
            debug_info->neighbour_cell = debug_info->current_cell +1;
            for ( int j = 0; j < n; j++ )
                debug_info->current_cell_vncons[j] = right_minus_ncons[j];
            for ( int j = 0; j < n; j++ )
                debug_info->neighbour_cell_vncons[j] = right_plus_ncons[j]; 
            for (int counter = 0; counter < M; counter++)
                printf("\n left_ncons = %lf", right_minus_ncons[counter]);
            printf("\n Left begin");
            printf("\n gas phase beginning %d", i);
            // решение задачи о распаде разрыва
            double c1 = calc_sound_velocity_one_phase(paramsc,right_minus_ncons, GAS_PHASE ); // скорость звука газовой фазы в граничной ячейке
            double c2 = calc_sound_velocity_one_phase(paramsc,right_minus_ncons, DISPERSED_PHASE ); // скорость звука дисперсной фазы в граничной ячейке
            if (right_plus_ncons[V_GAS] < right_minus_ncons[V_GAS])
            {
                
                cont_ncons[P_GAS] = right_minus_ncons[P_GAS] + (paramsc->g2 + 1)*right_minus_ncons[R_GAS] * (body_velocity - right_minus_ncons[V_GAS]) * (body_velocity - right_minus_ncons[V_GAS]) / 4 *( 1 + sqrt( 1 + (4 * c1 /( (paramsc->g2 + 1) * (body_velocity - right_minus_ncons[V_GAS]))) * (4 * c1 /( (paramsc->g2 + 1) * (body_velocity - right_minus_ncons[V_GAS])))));
                cont_ncons[V_GAS] = body_velocity;
                cont_ncons[R_GAS] = right_minus_ncons[R_GAS] * ( (paramsc->g2 - 1) * right_minus_ncons[P_GAS] + (paramsc->g2 + 1) * cont_ncons[P_GAS]) / ( (paramsc->g2 + 1) * right_minus_ncons[P_GAS] + (paramsc->g2 - 1) * cont_ncons[P_GAS]);
            }
            else
            {
                cont_ncons[P_GAS] = right_minus_ncons[P_GAS] * pow( 1 - (paramsc->g2 - 1) * (body_velocity - right_minus_ncons[V_GAS]) / ( 2 * c1 ), 2 * paramsc->g2/(paramsc->g2 - 1));
                cont_ncons[V_GAS] = body_velocity;
                cont_ncons[R_GAS] = right_minus_ncons[R_GAS] * pow(cont_ncons[P_GAS] / right_minus_ncons[P_GAS], 1/paramsc->g2);
            }

            printf("\n gas phase end %d", i);
            printf("\n Disperced phase beginning %d", i);
            if (right_plus_ncons[V_DISP] < right_minus_ncons[V_DISP])
            {
               printf("\n shock"); 
                cont_ncons[P_DISP] = right_minus_ncons[P_DISP] + (paramsc->g1 + 1)*right_minus_ncons[R_DISP] * (body_velocity - right_minus_ncons[V_DISP]) * (body_velocity - right_minus_ncons[V_DISP]) / 4 *( 1 + sqrt( 1 + (4 * c2 /( (paramsc->g1 + 1) * (body_velocity - right_minus_ncons[V_DISP]))) * (4 * c2 /( (paramsc->g1 + 1) * (body_velocity - right_minus_ncons[V_DISP])))));
                cont_ncons[V_DISP] = body_velocity;
                cont_ncons[R_DISP] = right_minus_ncons[R_DISP] * ( (paramsc->g1 - 1) * (right_minus_ncons[P_DISP] + paramsc-> p01) + (paramsc->g1 + 1) * (cont_ncons[P_DISP] + paramsc-> p01)) / ( (paramsc->g1 + 1) * (right_minus_ncons[P_DISP] + paramsc-> p01) + (paramsc->g1 - 1) * (cont_ncons[P_DISP] + paramsc-> p01));
            }
            else
            {
                printf("\n fan");
                cont_ncons[P_DISP] = (right_minus_ncons[P_DISP] + paramsc->p01) * pow( 1 - (paramsc->g1 - 1) * (body_velocity - right_minus_ncons[V_DISP])* (body_velocity - right_minus_ncons[V_DISP]) / ( 2 * c2 ), 2 * paramsc->g1/(paramsc->g1 - 1)) - paramsc->p01;
                cont_ncons[V_DISP] = body_velocity;
                cont_ncons[R_DISP] = right_minus_ncons[R_DISP] * pow((cont_ncons[P_DISP] + paramsc-> p01) / ( right_minus_ncons[P_DISP] + paramsc-> p01), 1/paramsc->g1);
            }

            //cont_ncons[P_DISP] = center_ncons[P_DISP];
            //cont_ncons[R_DISP] = center_ncons[R_DISP];
            //cont_ncons[V_DISP] = -center_ncons[V_DISP];
            cont_ncons[B_DISP] =  center_ncons[B_DISP];
            printf("\n Disperced phase end %d", i);
            for (int counter = 0; counter < M; counter++)
                cont_left[counter] = cont_ncons[counter];
            double contact_discontinuity_sound_velocity_gas, contact_discontinuity_sound_velocity_disp;
            calc_sound_velocity(paramsc, cont_ncons, &contact_discontinuity_sound_velocity_gas, &contact_discontinuity_sound_velocity_disp);            
 
            for(int l = 0; l < n; l++)
                printf("\n cont_ncons %lf", cont_ncons[l]);
            //выбор того, в какую зону попадает граница ячейки для газовой фазы
            // Формулы из Chertock, Kurganov
            if ( body_velocity > 0 )
            {
                if (body_velocity < contact_discontinuity_sound_velocity_gas)
                {
                    right_plus_ncons[P_GAS] = cont_ncons[P_GAS];
                    right_plus_ncons[R_GAS] = cont_ncons[R_GAS];
                    right_plus_ncons[V_GAS] = cont_ncons[V_GAS];
                }
                else
                {
                    right_plus_ncons[P_GAS] = right_minus_ncons[P_GAS];
                    right_plus_ncons[R_GAS] = right_minus_ncons[R_GAS];
                    right_plus_ncons[V_GAS] = right_minus_ncons[V_GAS];
                }
            }
            else
            {
                if (body_velocity > -contact_discontinuity_sound_velocity_gas)
                {
                    right_plus_ncons[P_GAS] = cont_ncons[P_GAS];
                    right_plus_ncons[R_GAS] = cont_ncons[R_GAS];
                    right_plus_ncons[V_GAS] = cont_ncons[V_GAS];
                }
                else
                {
                }
            }

            //выбор того, в какую зону попадает граница ячейки для дисперсной фазы
            if ( body_velocity > 0 )
            {
                if (body_velocity < contact_discontinuity_sound_velocity_disp)
                {
                    right_plus_ncons[P_DISP] = cont_ncons[P_DISP];
                    right_plus_ncons[R_DISP] = cont_ncons[R_DISP];
                    right_plus_ncons[V_DISP] = cont_ncons[V_DISP];
                }
                else
                {
                    right_plus_ncons[P_DISP] = right_minus_ncons[P_DISP];
                    right_plus_ncons[R_DISP] = right_minus_ncons[R_DISP];
                    right_plus_ncons[V_DISP] = right_minus_ncons[V_DISP];
                }
            }
            else
            {
                if (body_velocity > -contact_discontinuity_sound_velocity_disp)
                {
                    right_plus_ncons[P_DISP] = cont_ncons[P_DISP];
                    right_plus_ncons[R_DISP] = cont_ncons[R_DISP];
                    right_plus_ncons[V_DISP] = cont_ncons[V_DISP];
                }
                else
                {
                }
            }
            right_plus_ncons[B_DISP] = right_minus_ncons[B_DISP];
            if (body_velocity ==0 )
            {
                    right_plus_ncons[P_DISP] = right_minus_ncons[P_DISP];
                    right_plus_ncons[R_DISP] = right_minus_ncons[R_DISP];
                    right_plus_ncons[V_DISP] = -right_minus_ncons[V_DISP];
                    right_plus_ncons[P_GAS] = right_minus_ncons[P_GAS];
                    right_plus_ncons[R_GAS] = right_minus_ncons[R_GAS];
                    right_plus_ncons[V_GAS] = -right_minus_ncons[V_GAS];

            }
            hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
            hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars );
            
            printf("\n Left end");
            
        }
    else{}
        if ( GHOST == status[i-1]){

             for (int j = 0;  j < n; j++)
                left_minus_ncons[j] = left_plus_ncons[j];

            left_minus_ncons[V_GAS] = 2 * body_velocity - left_minus_ncons[V_GAS]; 
            left_minus_ncons[V_DISP] = 2 * body_velocity - left_minus_ncons[V_DISP];
            //convert_noncons_to_cons(paramsc, left_minus_ncons, left_minus_cons, number_of_scalars);
            double cont_ncons[M];
            cont_ncons[7] = left_plus_ncons[7];
            debug_info->neighbour_cell = debug_info->current_cell -11;
            for ( int j = 0; j < n; j++ )
                debug_info->current_cell_vncons[j] = left_plus_ncons[j];
            for ( int j = 0; j < n; j++ )
                debug_info->neighbour_cell_vncons[j] = left_minus_ncons[j]; 
            printf("\n Right begin");
            printf("\n gas phase beginning %d", i);
            // решение задачи о распаде разрыва
            double c1 = calc_sound_velocity_one_phase(paramsc,left_plus_ncons, GAS_PHASE ); // скорость звука газовой фазы в граничной ячейке
            double c2 = calc_sound_velocity_one_phase(paramsc,left_plus_ncons, DISPERSED_PHASE ); // скорость звука дисперсной фазы в граничной ячейке
            if (left_minus_ncons[V_GAS] > left_plus_ncons[V_GAS])
            {
                
                cont_ncons[P_GAS] = left_plus_ncons[P_GAS] + (paramsc->g2 + 1)*left_plus_ncons[R_GAS] * (body_velocity - left_plus_ncons[V_GAS]) * (body_velocity - left_plus_ncons[V_GAS]) / 4 *( 1 + sqrt( 1 + (4 * c1 /( (paramsc->g2 + 1) * (body_velocity - left_plus_ncons[V_GAS]))) * (4 * c1 /( (paramsc->g2 + 1) * (body_velocity - left_plus_ncons[V_GAS])))));
                cont_ncons[V_GAS] = body_velocity;
                cont_ncons[R_GAS] = left_plus_ncons[R_GAS] * ( (paramsc->g2 - 1) * left_plus_ncons[P_GAS] + (paramsc->g2 + 1) * cont_ncons[P_GAS]) / ( (paramsc->g2 + 1) * left_plus_ncons[P_GAS] + (paramsc->g2 - 1) * cont_ncons[P_GAS]);
            }
            else
            {
                cont_ncons[P_GAS] = left_plus_ncons[P_GAS] * pow( 1 + (paramsc->g2 - 1) * (body_velocity - left_plus_ncons[V_GAS]) / ( 2 * c1 ), 2 * paramsc->g2/(paramsc->g2 - 1));
                cont_ncons[V_GAS] = body_velocity;
                
                cont_ncons[R_GAS] = left_plus_ncons[R_GAS] * pow(cont_ncons[P_GAS] / left_plus_ncons[P_GAS], 1/paramsc->g2);
            }

            printf("\n gas phase end %d", i);
            printf("\n Disperced phase beginning %d", i);
            if (left_minus_ncons[V_GAS] > left_plus_ncons[V_GAS])
            {
                
                cont_ncons[P_DISP] = left_plus_ncons[P_DISP] + (paramsc->g1 + 1)*left_plus_ncons[R_DISP] * (body_velocity - left_plus_ncons[V_DISP]) * (body_velocity - left_plus_ncons[V_DISP]) / 4 *( 1 + sqrt( 1 + (4 * c2 /( (paramsc->g1 + 1) * (body_velocity - left_plus_ncons[V_DISP]))) * (4 * c2 /( (paramsc->g1 + 1) * (body_velocity - left_plus_ncons[V_DISP])))));
                cont_ncons[V_DISP] = body_velocity;
                cont_ncons[R_DISP] = left_plus_ncons[R_DISP] * ( (paramsc->g1 - 1) * (left_plus_ncons[P_DISP] + paramsc-> p01) + (paramsc->g1 + 1) * (cont_ncons[P_DISP] + paramsc-> p01)) / ( (paramsc->g1 + 1) * (left_plus_ncons[P_DISP]  + paramsc-> p01) + (paramsc->g1 - 1) * (cont_ncons[P_DISP] + paramsc-> p01));
            }
            else
            {
                cont_ncons[P_DISP] = (left_plus_ncons[P_DISP] + paramsc->p01) * pow( 1 + (paramsc->g1 - 1) * (body_velocity - left_plus_ncons[V_DISP]) / ( 2 * c2 ), 2 * paramsc->g1/(paramsc->g1 - 1)) - paramsc->p01;
                cont_ncons[V_DISP] = body_velocity;
                cont_ncons[R_DISP] = left_plus_ncons[R_DISP] * pow((cont_ncons[P_DISP] + paramsc-> p01) / ( left_plus_ncons[P_DISP] + paramsc-> p01), 1/paramsc->g1);
            }
           
            //cont_ncons[P_DISP] = center_ncons[P_DISP];
            //cont_ncons[R_DISP] = center_ncons[R_DISP];
            //cont_ncons[V_DISP] = -center_ncons[V_DISP];
            cont_ncons[B_DISP] =  center_ncons[B_DISP];


            printf("\n Disperced phase end %d", i);
            for (int counter = 0; counter < M; counter++)
                cont_right[counter] = cont_ncons[counter];
            double contact_discontinuity_sound_velocity_gas, contact_discontinuity_sound_velocity_disp;
            calc_sound_velocity(paramsc, cont_ncons, &contact_discontinuity_sound_velocity_gas, &contact_discontinuity_sound_velocity_disp);            
 
            for(int l = 0; l < n; l++)
                printf("\n cont_ncons %lf", cont_ncons[l]);
            //выбор того, в какую зону попадает граница ячейки для газовой фазы
            // если тело движется влево быстрее скорости звука на контактном разрыве
            if ( body_velocity < -contact_discontinuity_sound_velocity_gas )
            {
                left_minus_ncons[P_GAS] = left_plus_ncons[P_GAS];
                left_minus_ncons[R_GAS] = left_plus_ncons[R_GAS];
                left_minus_ncons[V_GAS] = left_plus_ncons[V_GAS];
            }
            // если тело движется медленнее скорости звука на контактном разрыве
            else if ( body_velocity < contact_discontinuity_sound_velocity_gas )
            {
                left_minus_ncons[P_GAS] = cont_ncons[P_GAS];
                left_minus_ncons[R_GAS] = cont_ncons[R_GAS];
                left_minus_ncons[V_GAS] = cont_ncons[V_GAS];
            }
          
            // если тело движется вправо быстрее скорости звука на контактном разрыве
            else
            {


            }
            //выбор того, в какую зону попадает граница ячейки для дисперсной фазы
            if ( body_velocity < -contact_discontinuity_sound_velocity_disp )
            {
                left_minus_ncons[P_DISP] = left_plus_ncons[P_DISP];
                left_minus_ncons[R_DISP] = left_plus_ncons[R_DISP];
                left_minus_ncons[V_DISP] = left_plus_ncons[V_DISP];
            }
            // если тело движется медленнее скорости звука на контактном разрыве
            else if ( body_velocity < contact_discontinuity_sound_velocity_gas )
            {
                left_minus_ncons[P_DISP] = cont_ncons[P_DISP];
                left_minus_ncons[R_DISP] = cont_ncons[R_DISP];
                left_minus_ncons[V_DISP] = cont_ncons[V_DISP];
            }
          
            // если тело движется вправо быстрее скорости звука на контактном разрыве
            else
            {


            }
            hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
            hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars );
            printf("\n right end");

        }
    else{}
     }
            
                if (status[i] == INNER){

                    hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars );
                    // расчет потока через правую границу ячейки

                    hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
                }
                else{}
                
               

       
        if (status[i] != OUTER && status[i] != GHOST){
        // аппроксимация неконсервативных составляющих в правых частях
        array1D rhst( M ); // неконсервативный вектор правых частей
        rhst_ncons( paramsc, center_ncons, &rhst, number_of_scalars );
        double dB2dx; // аппроксимация градиента объемной доли газовой фазы
        calc_dB2( left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, splus_l, sminus_l, splus_r, sminus_r, &dB2dx );
        dB2dx /= h;
        if (status[i] == BOUNDARY){
            dB2dx = 0;
        }
        
        // рассчитываем компоненты вектора-решения на следующем шаге по времени без первой компоненты (объемная доля дисперсной фазы)
        double center_cons[M]; // вектор консервативных переменных в текущей рассчитываемой ячейке
        convert_noncons_to_cons( paramsc, center_ncons, center_cons, number_of_scalars );
        double solution_cons[M]; // вектор-решение консервативных переменных на следующем шаге по времени
        for ( int l = 1; l < n; l++ ){
            solution_cons[l] = center_cons[l] - dt * ( right_flux[l] - left_flux[l] ) / h + dt * dB2dx * rhst[l];
            //printf("\n fluxes %lf %lf %lf", right_flux[l], left_flux[l], rhst[l]);
        }
        //printf("\n cons %lf %lf", solution_cons[P_GAS], center_cons[P_GAS]);
        // расчет переноса объемной доли дисперсной фазы


            double u_i = calc_u_i( paramsc, center_ncons );
            double t1 = ( u_i * ( splus_r * right_minus_ncons[B_DISP] - sminus_r * right_plus_ncons[B_DISP] ) + splus_r * sminus_r * ( right_plus_ncons[B_DISP] - right_minus_ncons[B_DISP] ) )
                / ( splus_r - sminus_r );
            double t2 = ( u_i * ( splus_l * left_minus_ncons[B_DISP] - sminus_l * left_plus_ncons[B_DISP] ) + splus_l * sminus_l * ( left_plus_ncons[B_DISP] - left_minus_ncons[B_DISP] ) )
                / ( splus_l - sminus_l );
            //printf("\n t %lf %lf", t1, t2);

            solution_cons[B_DISP] = center_ncons[B_DISP] - dt * ( t1 - t2 ) / h;
            if (solution_cons[B_DISP] > 0.9)
                solution_cons[B_DISP] = 0.9;

        convert_cons_to_noncons( paramsc, solution_cons, solution_ncons, number_of_scalars );

    }
        else
            for ( int l = 0; l < n; l++)
                solution_ncons[l] = center_ncons[l];
        //printf(" \n final sol %lf %lf ",center_ncons[B_DISP], solution_ncons[B_DISP]);
}

// Метод Harten - Lax - van Leer (HLL) расчета потоков в двухфазной среде
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных слева от разрыва
// right_ncons - вектор примитивных переменных справа от разрыва
// flux - рассчитанный вектор потока
// splus - оценка для скорости волны S+
// sminus - оценка для скорости волны S-
// n - реальный размер векторов
void hll_flux( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M],
               array1D* flux, double* splus, double* sminus, int n, int number_of_scalars ) {

    // расчет векторов дифференциального потока по параметрам слева и справа от разрыва
    array1D left_diff_flux( M );
    diff_flux_ncons( paramsc, left_ncons, &left_diff_flux, number_of_scalars );
    array1D right_diff_flux( M );
    diff_flux_ncons( paramsc, right_ncons, &right_diff_flux, number_of_scalars );
    // расчет оценок для скоростей волн S+ и S- для газовой фазы
    double splus_g, sminus_g;
    calc_splus_sminus( paramsc, left_ncons, right_ncons, GAS_PHASE, &splus_g, &sminus_g );

    // расчет оценок для скоростей волн S+ и S- для дисперсной фазы
    double splus_s, sminus_s;
    calc_splus_sminus( paramsc, left_ncons, right_ncons, DISPERSED_PHASE, &splus_s, &sminus_s );

    // расчет оценок для скоростей волн S+ и S- для полной системы
    *splus = max( splus_g, splus_s );
    *sminus = min( sminus_g, sminus_s );

    // расчет векторов консервативных переменных слева и справа от разрыва
    double left_cons[M];
    convert_noncons_to_cons( paramsc, left_ncons, left_cons, number_of_scalars );
    double right_cons[M];
    convert_noncons_to_cons( paramsc, right_ncons, right_cons, number_of_scalars );
	
    for ( int i = 1; i < n; i++ )
        (*flux)[i] = ( (*splus) * left_diff_flux[i] - (*sminus) * right_diff_flux[i] +
            (*splus) * (*sminus) * ( right_cons[i] - left_cons[i] ) ) / ( (*splus) - (*sminus) );

}

// Расчет оценок для скоростей волн S+ и S-
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons[M] - вектор примитивных переменных слева от разрыва
// right_ncons[M] - вектор примитивных переменных справа от разрыва
// phase - фаза, для которой оцениваются скорости
// splus - оценка для скорости волны S+
// sminus - оценка для скорости волны S-
void calc_splus_sminus( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M],
                        const Phase phase, double* splus, double* sminus ) {
                           
    double c_l = calc_sound_velocity_one_phase( paramsc, left_ncons, phase ); // скорость звука слева от разрыва
    double c_r = calc_sound_velocity_one_phase( paramsc, right_ncons, phase ); // скорость звука справа от разрыва

    double v_l, v_r;
    if ( phase == GAS_PHASE ) {
        v_l = left_ncons[V_GAS];
        v_r = right_ncons[V_GAS];
    }
    else {
        v_l = left_ncons[V_DISP];
        v_r = right_ncons[V_DISP];
    }

    // В статье:
    // Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
    // and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467,
    // написано так:
    // *splus = max( 0.0, max( v_l + c_l, v_r + c_r ) );
    // *sminus = min( 0.0, min( v_l - c_l, v_r - c_r ) );
    // Но в оригинальной статье, на которую эта статья ссылается:
    // Davis S.F. Simplified second-order Godunov-type methods // SIAM J. Sci. Stat. Comput. - 1988.
    // - V. 9, No. 3. - P. 445 - 473,
    // как и в книге Торо (издание 2009 г., стр. 328), написано так:
    *splus = max( v_l + c_l, v_r + c_r );
    *sminus = min( v_l - c_l, v_r - c_r );
    // В дозвуковом случае разницы быть не должно, иначе - возможна

}

// Аппроксимация дифференциала dB2 в правой части
// left_minus - левая грань, вектор консервативных переменных слева от разрыва
// left_plus - левая грань, вектор консервативных переменных справа от разрыва
// right_minus - правая грань, вектор консервативных переменных слева от разрыва
// right_plus - правая грань, вектор консервативных переменных справа от разрыва
// splus_l - оценка для скорости волны S+ на левом ребре
// sminus_l - оценка для скорости волны S- на левом ребре
// splus_r - оценка для скорости волны S+ на правом ребре
// sminus_l - оценка для скорости волны S- на правом ребре
// dB2dx - искомая аппроксимация
void calc_dB2( const double left_minus[M], const double left_plus[M], const double right_minus[M], const double right_plus[M],
               double splus_l, double sminus_l, double splus_r, double sminus_r, double* dB2dx ) {

    *dB2dx = ( ( splus_r * right_minus[B_DISP] - sminus_r * right_plus[B_DISP] ) / ( splus_r - sminus_r )
        - ( splus_l * left_minus[B_DISP] - sminus_l * left_plus[B_DISP] ) / ( splus_l - sminus_l ) );

}

