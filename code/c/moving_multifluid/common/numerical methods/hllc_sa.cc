// hllc_sa.cc
// Метод HLLC численного интегрирования уравнений Saurel-Abgrall
// Реализовано по: Li Q. et al. Difference scheme for two-phase flow // Applied Mathematics and Mechanics. - 2004. - V. 25, No. 5. - P. 536 - 545
// и
// Liang S. et al. Solving seven-equation model for compressible two-phase flow using multiple GPUs // Computers & Fluids. - 2014. - V. 99. - P. 156 - 171
// (c) Уткин Павел, 2017
// Создан: 29 марта 2017 г.

#include "hllc_sa.h"


// Метод HLLC численного интегрирования уравнений Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге
// n - реальный размер векторов
void hllc_1d( struct ParametersCommon* paramsc, struct Parameters1d* params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
          double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure, double body_velocity, int i, int *status, double cont_left[M], double cont_right[M]  ) {
    double solution_ncons_Lh[M]; // вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
    if (status[i+1] == BOUNDARY)
        printf("\n before hll %lf",center_ncons[P_GAS]);
          Lh_HLL(debug_info, params1d, paramsc, left_ncons, center_ncons, right_ncons, slopes_left, slopes_center, slopes_right,
          dt, h, solution_ncons_Lh, n, number_of_scalars, i, status, body_velocity,  cont_left, cont_right ); 
              if (status[i+1] == BOUNDARY)
        printf("\n after hll %lf",solution_ncons_Lh[P_GAS]);// действие гиперболического оператора
    if ( paramsc->pressure_relaxation == true && is_pressure_relaxation_after_this_step == true ){
        if (status[i] == INNER){

        Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // действие релаксационного оператора
        if (status[i+1] == BOUNDARY)
        printf("\n after rel %lf",solution_ncons[P_GAS]);
        }
        else{}

        if ( status[i] == BOUNDARY){
            if (status[i+1] == GHOST){
                Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh,cont_left , dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // действие релаксационного оператора
            }
        
            else {
                Lr( paramsc, params1d, debug_info, cont_right, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // действие релаксационного оператора
            }
        }
        else{}
    if (status[i] == GHOST || status[i] == OUTER){
        
       for ( int j = 0; j < n; j++ )
            solution_ncons[j] = solution_ncons_Lh[j];
        }
    
    else{}
    }
    else {
        for ( int j = 0; j < n; j++ )
            solution_ncons[j] = solution_ncons_Lh[j];
    }
                if (solution_ncons[B_DISP] > 0.9)
                solution_ncons[B_DISP] = 0.9;

}

// Гиперболический оператор метода HLLC
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
void Lh_HLLC(struct DebugInfo *debug_info, struct Parameters1d* params1d,  struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars,  int i, int *status, double body_velocity, double cont_left[M], double cont_right[M] ) {

    // реконструкция сеточных функций
    double left_minus_ncons[M];
    double left_minus_cons[M];
    double left_plus_ncons[M];
    double left_plus_cons[M];
    double right_minus_ncons[M];
    double right_minus_cons[M];
    double right_plus_ncons[M];
    double right_plus_cons[M];
    for ( int j = 0; j < M; j++ ) {
        // реконструкция значений на левой границе слева
        left_minus_cons[j] += 0.5 * h * slopes_left[j];
        // реконструкция значений на левой границе справа
        left_plus_cons[j] -= 0.5 * h * slopes_center[j];
        // реконструкция значений на правой границе слева
        convert_noncons_to_cons( paramsc, center_ncons, right_minus_cons, 0 );
        right_minus_cons[j] += 0.5 * h * slopes_center[j];
        // реконструкция значений на правой границе справа
        right_plus_cons[j] -= 0.5 * h * slopes_right[j];
        
    }
convert_noncons_to_cons( paramsc, right_ncons, right_plus_cons, 0 );
convert_cons_to_noncons( paramsc, right_plus_cons, right_plus_ncons, 0 );
convert_cons_to_noncons( paramsc, right_minus_cons, right_minus_ncons, 0 );
convert_cons_to_noncons( paramsc, left_plus_cons, left_plus_ncons, 0 );
convert_noncons_to_cons( paramsc, center_ncons, left_plus_cons, 0 );
convert_cons_to_noncons( paramsc, left_minus_cons, left_minus_ncons, 0 );
convert_noncons_to_cons( paramsc, left_ncons, left_minus_cons, 0 );
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
    //array1D left_flux( M );
    double left_phi; // объемная доля дисперсной фазы на левом ребре для аппроксимации неконсервативного члена
    double s_left; // оценка для скорости контактного разрыва

    //array1D right_flux( M );
    double right_phi; // объемная доля дисперсной фазы на правом ребре для аппроксимации неконсервативного члена
    double s_right; // оценка для скорости контактного разрыва
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

                    //hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars );
                    hllc_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &left_phi, &s_left );
                    // расчет потока через правую границу ячейки

                    //hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
                    hllc_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &right_phi, &s_right );
                }
                else{}
                
               

       

    
if (status[i] != OUTER && status[i] != GHOST){
    // аппроксимация неконсервативных составляющих в правых частях
    array1D rhst( M ); // неконсервативный вектор правых частей
    rhst_ncons( paramsc, center_ncons, &rhst, 0 );
    double theta = ( right_phi - left_phi ) / h;

    // рассчитываем компоненты вектора-решения на следующем шаге по времени без первой компоненты (объемная доля дисперсной фазы)
    double center_cons[M]; // вектор консервативных переменных в текущей рассчитываемой ячейке
    convert_noncons_to_cons( paramsc, center_ncons, center_cons, 0 );
    double solution_cons[M]; // вектор-решение консервативных переменных на следующем шаге по времени
    for ( int i = 1; i < M; i++ )
        solution_cons[i] = center_cons[i] - dt * ( right_flux[i] - left_flux[i] ) / h + dt * theta * rhst[i];

    // расчет переноса объемной доли дисперсной фазы
    
    // подход S. Liang, W. Liu, L. Yuan
    double u_i = calc_u_i( paramsc, center_ncons );

    // возможно, в оригинальной работе имелось в виду и так
    /* if ( u_i >= 0.0 ) {
        right_phi = center_ncons[B_DISP];
        left_phi = left_ncons[B_DISP];
    }
    else {
        right_phi = right_ncons[B_DISP];
        left_phi = center_ncons[B_DISP];
    } */
    
    solution_cons[B_DISP] = center_ncons[B_DISP] - dt * u_i * ( right_phi - left_phi ) / h;
    
    // подход Q. Li et al.
    /* double t1 = - u_i * ( right_phi - left_phi );
    double right_phi_star = ( sqrt( center_ncons[B_DISP] * center_ncons[R_DISP] ) * center_ncons[B_DISP] +
        sqrt( right_ncons[B_DISP] * right_ncons[R_DISP] ) * right_ncons[B_DISP] ) / ( sqrt( center_ncons[B_DISP] * center_ncons[R_DISP] ) + sqrt( right_ncons[B_DISP] * right_ncons[R_DISP] ) );
    double left_phi_star = ( sqrt( center_ncons[B_DISP] * center_ncons[R_DISP] ) * center_ncons[B_DISP] +
        sqrt( left_ncons[B_DISP] * left_ncons[R_DISP] ) * left_ncons[B_DISP] ) / ( sqrt( center_ncons[B_DISP] * center_ncons[R_DISP] ) + sqrt( left_ncons[B_DISP] * left_ncons[R_DISP] ) );
    double t2 = - ( s_right * ( right_phi_star - right_phi ) - s_left * ( left_phi_star - left_phi ) );
    solution_cons[B_DISP] = center_ncons[B_DISP] + dt * ( t1 + t2 ) / h; */

    convert_cons_to_noncons( paramsc, solution_cons, solution_ncons, 0 );
        }
        else
            for ( int l = 0; l < n; l++)
                solution_ncons[l] = center_ncons[l];

}

// Метод Harten - Lax - van Leer - Contact (HLLC) расчета потоков в двухфазной среде
// Liang S. et al. Solving seven-equation model for compressible two-phase flow using multiple GPUs //
// Computers & Fluids. - 2014. - V. 99. - P. 156 - 171.
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных слева от разрыва
// right_ncons - вектор примитивных переменных справа от разрыва
// flux - рассчитанный вектор потока
// phi - объемная доля дисперсной фазы на ребре
// s_cont - оценка для скорости контактного разрыва
void hllc_flux( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M], array1D* flux, double* phi, double* s_cont ) {

    // расчет векторов дифференциального потока по параметрам слева и справа от разрыва
    array1D left_diff_flux( M );
    diff_flux_ncons( paramsc, left_ncons, &left_diff_flux, 0 );
    array1D right_diff_flux( M );
    diff_flux_ncons( paramsc, right_ncons, &right_diff_flux, 0 );

    // расчет оценок для скоростей волн s_L и s_R для газовой фазы
    double s_L_g, s_R_g;
    calc_splus_sminus( paramsc, left_ncons, right_ncons, GAS_PHASE, &s_R_g, &s_L_g );

    // расчет оценок для скоростей волн s_L и s_R для дисперсной фазы
    double s_L_s, s_R_s;
    calc_splus_sminus( paramsc, left_ncons, right_ncons, DISPERSED_PHASE, &s_R_s, &s_L_s );

    // расчет оценок для скоростей волн s_L и s_R для полной системы
    double s_L = min( s_L_g, s_L_s );
    double s_R = max( s_R_g, s_R_s );
    
    // полученные оценки для s_L и s_R аналогичны методу HLL

    // расчет оценки для s*

    // расчет средних значений на ребре слева
    double r_m_L = calc_r_i( paramsc, left_ncons );
    double p_m_L = calc_p_i( paramsc, left_ncons );
    double u_m_L = calc_u_i( paramsc, left_ncons );

    // расчет средних значений на ребре справа
    double r_m_R = calc_r_i( paramsc, right_ncons );
    double p_m_R = calc_p_i( paramsc, right_ncons );
    double u_m_R = calc_u_i( paramsc, right_ncons );

    double s_star = ( p_m_R - p_m_L + r_m_L * u_m_L * ( s_L - u_m_L ) - r_m_R * u_m_R * ( s_R - u_m_R ) ) /
        ( r_m_L * ( s_L - u_m_L ) - r_m_R * ( s_R - u_m_R ) );

    // оценка параметров слева и справа от контактного разрыва
    array1D q_star_L( M );
    array1D q_star_R( M );

    calc_contact_vector( paramsc, left_ncons, s_L, s_star, &q_star_L );
    calc_contact_vector( paramsc, right_ncons, s_R, s_star, &q_star_R );

    // собственно расчет потока
    
    double left_cons[M];
    convert_noncons_to_cons( paramsc, left_ncons, left_cons, 0 );
    double right_cons[M];
    convert_noncons_to_cons( paramsc, right_ncons, right_cons, 0 );
    
    if ( s_L >= 0.0 ) {
        for ( int i = 1; i < M; i++ )
            (*flux)[i] = left_diff_flux[i];
        *s_cont = 0.0;
    }
    else if ( s_L <= 0.0 && s_star >= 0.0 ) {
        for ( int i = 1; i < M; i++ )
            (*flux)[i] = left_diff_flux[i] + s_L * ( q_star_L[i] - left_cons[i] );
        *s_cont = s_L;
    }
    else if ( s_star <= 0.0 && s_R >=0 ) {
        for ( int i = 1; i < M; i++ )
            (*flux)[i] = right_diff_flux[i] + s_R * ( q_star_R[i] - right_cons[i] );
        *s_cont = s_R;
    }
    else if ( s_R <= 0 ) {
        for ( int i = 1; i < M; i++ )
            (*flux)[i] = right_diff_flux[i];
        *s_cont = 0.0;
    }

    // значение объемной доли дисперсной фазы на ребре
    if ( s_star >= 0.0 )
        *phi = left_ncons[B_DISP];
    else
        *phi = right_ncons[B_DISP];

}

// Расчет вектора - оценки состояния слева или справа от контактного разрыва в методе HLLC
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных слева или справа от разрыва
// s - оценка скорости левой или правой волны
// s_star - оценка для s со звездочкой
// q_star - оценки состояния слева или справа от контактного разрыва
void calc_contact_vector( const struct ParametersCommon *paramsc, const double v_ncons[M], const double s, const double s_star, array1D* q_star ) {

    if ( fabs( s - s_star ) < paramsc->eps_general || fabs( s - v_ncons[V_DISP] ) < paramsc->eps_general ) {
        printf( "calc_contact_vector -> something is wrong with hllc method, s is equal to s_star\n" );
        exit( EXIT_FAILURE );
    }

    double v_cons[M];
    convert_noncons_to_cons( paramsc, v_ncons, v_cons, 0 );

    // дисперсная фаза

    double mult = v_cons[R_DISP] * ( s - v_ncons[V_DISP] ) / ( s - s_star );

    (*q_star)[R_DISP] = mult;
    (*q_star)[V_DISP] = mult * s_star;
    (*q_star)[P_DISP] = mult * (  v_cons[P_DISP] /  v_cons[R_DISP] + ( s_star - v_ncons[V_DISP] ) * ( s_star + v_ncons[P_DISP] / v_ncons[R_DISP] / ( s - v_ncons[V_DISP] ) ) );

    // газовая фаза

    mult = v_cons[R_GAS] * ( s - v_ncons[V_GAS] ) / ( s - s_star );

    (*q_star)[R_GAS] = mult;
    (*q_star)[V_GAS] = mult * s_star;
    (*q_star)[P_GAS] = mult * (  v_cons[P_GAS] /  v_cons[R_GAS] + ( s_star - v_ncons[V_GAS] ) * ( s_star + v_ncons[P_GAS] / v_ncons[R_GAS] / ( s - v_ncons[V_GAS] ) ) );

}