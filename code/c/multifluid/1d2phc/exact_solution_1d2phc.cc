// exact_solution_1d2phc.cc
// Построение точного решения задачи о распаде разрыва для системы уравнений Баера - Нунциато.
// (c) Уткин Павел, 2013
// Создан: 24 августа 2013 г.

#include "exact_solution_1d2phc.h"

// Построение точного решения задачи о распаде разрыва для системы уравнений Баера-Нунзиато
// output_directory - директория, куда будет записан файл с точным решением (in) 
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с параметрами одномерной задачи (in)
// debug_info - структура с отладочной информацией (in)
void build_exact_sol( char *output_directory, struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, int n ) {

    char output_filename[MAX_STRING_SIZE]; // имя выходного файла файла
    double *xc; // массив координат центров ячеек сетки для точного решения
    double *x; // массив координат граней ячеек сетки для точного решения
    double s; // текущее значение автомодельной переменной
    FILE *ex_sol_out;
    double flux[M]; // вектор потока - фактически не используется, но нужно для использования функции godunov_cons_flux
    
    double p_solid_cont_l, p_solid_cont_r; // давления в дисперсной фазе слева и справа от контактного разрыва в дисперсной фазе
    double v_solid_cont; // скорость контактного разрыва в дисперсной фазе
   
    double v_ncons_res[M]; // вектор-решение задачи о распаде разрыва

    double time_moment; // момент времени, для которого строится точное решение

    Disp_phase_cases solver_part; // константа, которая определяет, какую из четырех частей солвера использовать для расчета "потока"
    
    // определение момента времени, для которого строится точное решение, исходя из опций задачи
    if ( paramsc->is_output_on_time )
        time_moment = paramsc->stop_time;
    else
        time_moment = paramsc->stop_steps_number * paramsc->dt;

    // формирование имени выходного файла с точным решением
    strcpy_s( output_filename, output_directory );
    strcat_s( output_filename, "\\exact_solution.dat" );
    fopen_s( &ex_sol_out, output_filename, "wt" );
    if ( NULL == ex_sol_out ) {
        printf( "build_exact_sol -> Can't open file %s for writing.\n", output_filename );
    }

    // выделение памяти под сетку для точного решения
    get_memory_for_1D_double_array( params1d->cells_number_for_exact_solution, &xc );
    get_memory_for_1D_double_array( params1d->cells_number_for_exact_solution + 1, &x );

    // построение расчетной сетки
    build_grid( LEFT_BOUN_EX_SOL, RIGHT_BOUN_EX_SOL, params1d->cells_number_for_exact_solution, xc, x );

    // анализ соотношения между величинами объемной доли дисперсной фазы слева и справа от разрыва
    solver_part = what_case( paramsc, (params1d->block_values[0])[B_DISP], (params1d->block_values[1])[B_DISP] );
    double gas_left_ncons_reduced[M_REDUCTION]; // сокращенный вектор для параметров газовой фазы слева от разрыва
    double gas_right_ncons_reduced[M_REDUCTION]; // сокращенный вектор для параметров газовой фазы справа от разрыва
    double disp_left_ncons_reduced[M_REDUCTION]; // сокращенный вектор для параметров дисперсной фазы слева от разрыва
    double disp_right_ncons_reduced[M_REDUCTION]; // сокращенный вектор для параметров дисперсной фазы справа от разрыва
    if ( solver_part == NO_GRAD ) {
        // формируем сокращенные векторы
        convert_full_to_reduced( params1d->block_values[0], GAS_PHASE, gas_left_ncons_reduced );
        convert_full_to_reduced( params1d->block_values[1], GAS_PHASE, gas_right_ncons_reduced );
        convert_full_to_reduced( params1d->block_values[0], DISPERSED_PHASE, disp_left_ncons_reduced );
        convert_full_to_reduced( params1d->block_values[1], DISPERSED_PHASE, disp_right_ncons_reduced );
    }

    // цикл по ячейкам
    for ( int i_cell = 0; i_cell < params1d->cells_number_for_exact_solution; i_cell++ ) {
        s = xc[i_cell] / time_moment; // текущее значение автомодельной переменной
        // запись координаты в файл
        fprintf( ex_sol_out, "%e ", xc[i_cell] + ( params1d->left_boundary_x - LEFT_BOUN_EX_SOL ) ); // преобразуем в соответствии с сновной расчетной областью
        if ( solver_part == NO_GRAD ) { // случай полного расщепления фаз
                        
            // строим решение для дисперсной фазы
            double cl = calc_sound_velocity_reduced( paramsc, disp_left_ncons_reduced, DISPERSED_PHASE );
            double cr = calc_sound_velocity_reduced( paramsc, disp_right_ncons_reduced, DISPERSED_PHASE );
            // итерационная процедура расчета давления и скорости дисперсной фазы на контактном разрыве
            double p_cont; // давление на контакном разрыве
            double v_cont; // скорость на контактном разрыве
            calc_contact_pressure_velocity( paramsc, debug_info, disp_left_ncons_reduced, disp_right_ncons_reduced, M_REDUCTION,
                cl, cr, DISPERSED_PHASE, &p_cont, &v_cont );
            // отбор решения
            double v_ncons_disp_reduced[M_REDUCTION]; // сокращенный вектор
            sample_reduced( paramsc, disp_left_ncons_reduced, disp_right_ncons_reduced, M_REDUCTION, cl, cr, DISPERSED_PHASE,
                p_cont, v_cont, s, v_ncons_disp_reduced );
            
            // строим решение для газовой фазы
            cl = calc_sound_velocity_reduced( paramsc, gas_left_ncons_reduced, GAS_PHASE );
            cr = calc_sound_velocity_reduced( paramsc, gas_right_ncons_reduced, GAS_PHASE );
            // итерационная процедура расчета давления и скорости газа на контактном разрыве
            calc_contact_pressure_velocity( paramsc, debug_info, gas_left_ncons_reduced, gas_right_ncons_reduced, M_REDUCTION,
                cl, cr, GAS_PHASE, &p_cont, &v_cont );
            // отбор решения
            double v_ncons_gas_reduced[M_REDUCTION]; // сокращенный вектор
            sample_reduced( paramsc, gas_left_ncons_reduced, gas_right_ncons_reduced, M_REDUCTION, cl, cr, GAS_PHASE,
                p_cont, v_cont, s, v_ncons_gas_reduced );

            // формируем полный вектор-решение
            v_ncons_res[B_DISP] = (params1d->block_values[0])[B_DISP];
            convert_reduced_to_full( v_ncons_disp_reduced, DISPERSED_PHASE, v_ncons_res );
            convert_reduced_to_full( v_ncons_gas_reduced, GAS_PHASE, v_ncons_res );

        }
        else { // есть разрыв объемной доли дисперсной фазы
            ReturnCodes code = godunov_cons_flux( paramsc, debug_info, params1d->block_values[0], params1d->block_values[1], s, solver_part,
                flux, v_ncons_res, &v_solid_cont, &p_solid_cont_l, &p_solid_cont_r, n );
        }
        // записываем полный вектор-решение
        for ( int i_component = 0; i_component < n; i_component++ )
            fprintf( ex_sol_out, "%e ", v_ncons_res[i_component] );
        fprintf( ex_sol_out, "\n" );
    }

    // свобождение памяти
    free( xc );
    free( x );

    fclose( ex_sol_out );

}