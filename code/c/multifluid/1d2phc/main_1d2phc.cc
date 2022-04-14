// main.cc
// Одномерный (1D) двухфазный (2PH) сжимаемый (C) код
// (c) Уткин Павел, 2013 - 2018
// Создан: 16 февраля 2013 г.

#include "main_1d2phc.h"

// общие для всех методов функции
#include "grid.h"
#include "utils.h"
#include "utils_1d2phc.h"
#include "io_1d.h"
#include "exact_solution_1d2phc.h"
#include "memory.h"

// "физика"
#include "physics_solver.h"
#include "source_terms.h"

// повышение порядка аппроксимации
#include "minmod.h"

int main( int argc, char *argv[] ) {

    double initial_total_mass; // исходная суммарная масса вещества в системе
    double mass_diff; // относительная погрешность суммарной массы вещества в системе относительно исходной
    double time_gap; // время между записями файлов с промежуточными результатами, имеет смысл при params.is_output_on_time = 1
    int steps_gap; // количество шагов по времени между записями файлов с промежуточными результатами, имеет смысл при params.is_output_on_time = 0
    int files_counter = 0; // счетчик файлов с промежуточными результатами
    char output_filename[MAX_STRING_SIZE]; // имя текущего файла с промежуточными результатами
    char tmp_str[MAX_STRING_SIZE]; // строковая переменная для перевода числовых данных

    printf( "\n1D two-phase compressible solver\n(c) Pavel Utkin, ICAD RAS, MIPT, 2013-2018\ne-mail: pavel_utk@mail.ru\n" );

    printf( "\nPreparing:\n" );

    // обработка параметров командной строки
    if ( argc != 3 ) {
        // argv[1] - путь к директории, где находится файл с параметрами задачи
        // argv[2] - путь к директории, где будут записаны выходные файлы
        printf( "\nmain -> wrong command line arguments number\n" );
        exit( EXIT_FAILURE );
    }

    // считывание файла с параметрами задачи и заполнение структуры params
    int file_num = 0;
    struct ParametersCommon paramsc; // структура с основными параметрами вычислительного эксперимента
    struct Parameters1d params1d; // структура с параметрами одномерной задачи
    struct Parameters2d empty_structure; // структура с параметрами 2d задачи, не заполняемая в этой задаче
    fill_parameters_struct_1d( &paramsc, &params1d, argv[1], argv[2], file_num );
    printf( "\n> File Parameters1d.dat is processed successfully. Correspondent structure is filled.\n" );

    // выделение памяти под массивы
    double *xc; // массив координат центров ячеек сетки
    get_memory_for_1D_double_array( params1d.cells_number, &xc );
    double *x; // массив координат узлов сетки
    get_memory_for_1D_double_array( params1d.cells_number + 1, &x );
    double **u_prev; // вектора примитивных переменных на n-ом слое
    get_memory_for_2D_double_array( params1d.cells_number, M, &u_prev );
    double **u_next; // вектора примитивных переменных на (n+1)-ом слое
    get_memory_for_2D_double_array( params1d.cells_number, M, &u_next );
    double **slopes; // наклоны всех компонент векторов консервативных переменных во всех ячейках расчетной области для повышения порядка аппроксимации
    get_memory_for_2D_double_array( params1d.cells_number, M, &slopes );
    double *beta; // configuration pressure
    get_memory_for_1D_double_array( params1d.cells_number, &beta );
    double *B; // энергия компатирования
    get_memory_for_1D_double_array( params1d.cells_number, &B );
    double *C; // влияние горения
    get_memory_for_1D_double_array( params1d.cells_number, &C );
    printf( "\n> All the necessary memory is allocated successfully.\n" );
    int *number_of_block;
    get_memory_for_1D_int_array( 1, &number_of_block );
    int *ignition_flag; // индикатор начатого горения в ячейке
    get_memory_for_1D_int_array( params1d.cells_number, &ignition_flag );
    int *initial_ignition; // индикатор начального поджига
    get_memory_for_1D_int_array( params1d.cells_number, &initial_ignition );
    double *S; // Площадь канала
    get_memory_for_1D_double_array( params1d.cells_number, &S );
    // приведение размерных параметров задачи к безразмерному виду
    if ( paramsc.use_dimensions ) {
        dimensionalization( &paramsc, &params1d, &empty_structure );
        printf( "\n> Dimensionalization is done.\n" );
    }

    // определение координат центров ячеек сетки
    build_grid( params1d.left_boundary_x, params1d.right_boundary_x, params1d.cells_number, xc, x );
    printf( "\n> Computational grid is prepared successfully.\n" );

    // инициализация вектора-решения
    struct TimeMoment time_mom = { 0, 0.0 }; // структура, определяющая текущий момент времени
    init_solution_1d( &paramsc, &params1d, argv[1], &time_mom, u_prev );
    printf( "\n> The solution is initiated successfully.\n" );

    for (int i = 0; i < params1d.cells_number; i++){
	// определение блока, в котором сейчас находимся
	number_of_block[0] = 0;
	for (int j = 0; j < params1d.ic_blocks_number; j++){
	    if (i <= params1d.cell_end[j] && i >= params1d.cell_begin[j])
                number_of_block[0] = j;
	}
	ignition_flag[i] = 0;
	initial_ignition[i] = params1d.initial_burning[number_of_block[0]];
        beta[i] = calc_configuration_pressure(&paramsc, &params1d, &empty_structure, X_DIRECTION, u_prev[i], number_of_block);
	B[i] = compaction_energy(&paramsc, &params1d, &empty_structure, X_DIRECTION, u_prev[i], number_of_block);
        C[i] = calc_chemical_reaction(&paramsc, u_prev[i],ignition_flag[i]);
        S[i] = params1d.S[number_of_block[0]];
    }
    
    double curr_inflow_params[M]; // параметры втекания в данный момент в данной ячейке

    // Запись файла с начальными распределениями
    write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, 0, 0, output_filename, beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );

    // подготовительные процедуры для работы с датчиками в заданных местах
    int sensors_cells[MAX_SENSORS_NUM]; // для каждого датчика хранится номер ячейки, в которой он находится
    FILE *sensors_files[MAX_SENSORS_NUM]; // массив файловых дескрипторов для записи данных с датчиков
    if ( paramsc.are_sensors ) {
        prepare_sensors( &paramsc, &params1d, x, sensors_cells, sensors_files );
    }

    // инициализация полей структуры debug_info
    struct DebugInfo debug_info; // структура с отладочной информацией
    init_debug_info( argv[2], &debug_info );

    // построение точного решения
    if ( params1d.build_exact_solution ) {   
        build_exact_sol( argv[2], &paramsc, &params1d, &debug_info, M1D );
        printf( "\n> The exact solution is constructed successfully.\n" );
    }

    printf( "\nCalculation:\n" );

    // расчет исходного значения суммарной массы вещества
    initial_total_mass = get_total_mass( &params1d, x, u_prev );

    // промежуток времени или количество шагов между записями файлов с промежуточными результатами
    if ( paramsc.is_output_on_time )
        time_gap = paramsc.stop_time / paramsc.output_number;
    else
        steps_gap = paramsc.stop_steps_number / paramsc.output_number;

    // обнулили наклоны  перед началом расчета
    for ( int i = 0; i < params1d.cells_number; i++ ){
        for ( int j = 0; j < M1D + params1d.number_of_scalars; j++ )
            slopes[i][j] = 0.0;
	}
    double zero_slopes[M]; // массив нулевых наклонов для потоков на границах
    for ( int i = 0; i < M; i++ )
        zero_slopes[i] = 0.0;

    double boun_v[M]; // вектор примитивных переменных для реализации граничного условия
    double dt; // шаг по времени

    // основной цикл по времени
    while ( ( ( paramsc.stop_steps_number - time_mom.steps_num > 0 ) && ( !paramsc.is_output_on_time ) ) ||
          ( ( paramsc.stop_time - time_mom.curr_t > 0 ) && ( paramsc.is_output_on_time ) ) ) {
	
        // подготовка полей отладочной структуры для анализа процедуры релаксации
        for ( int i_relax_case = 0; i_relax_case < RELAX_CASES_NUM; i_relax_case++ )
            debug_info.relaxation_cases[i_relax_case] = 0;
        // подготовка полей отладочной структуры для анализа работы метода Годунова
        for ( int i_godunov_case = 0; i_godunov_case < GODUNOV_CASES_NUM; i_godunov_case++ )
            debug_info.godunov_cases[i_godunov_case] = 0;

        // расчет текущего шага по времени
        if ( time_mom.steps_num >= 188 ){
            dt = get_time_step(&debug_info, &paramsc, &params1d, x, u_prev, 1 );
        }
        else {
            dt = get_time_step(&debug_info, &paramsc, &params1d, x, u_prev, 0 );
        }

        // расчет наклонов сеточных функций
        // дальнейшая реконструкция сеточных функций проводится в предположении, что сетка равномерная
        if ( paramsc.approximation_order == 2 )
            reconstruction( &paramsc, &params1d, u_prev, xc, slopes, params1d.number_of_scalars );

        // цикл по ячейкам
        for ( int i = 0; i < params1d.cells_number; i++ ) {
            //printf("\n%lf", u_next[i][P_GAS]);

            // заполнение части полей структуры для отладки кода, касающихся данных о текущей ячейке
           /* debug_info.current_cell = i;
            debug_info.current_cell_x = xc[i];
            for ( int j = 0; j < M1D + params1d.number_of_scalars; j++ )
                debug_info.current_cell_vncons[j] = u_prev[i][j];*/

            // расчет параметров в i-ой ячейке на следующем временном слое
            if ( i == 0 ) {
                // обработка левого граничного условия
                current_inflow_parameters(&paramsc, &params1d, &empty_structure, LEFT_BOUNDARY, curr_inflow_params); 
                boundary( &params1d, &empty_structure, u_prev[0], params1d.left_bc, boun_v, time_mom.curr_t, curr_inflow_params );
                calc_step_in_cell( &paramsc, &params1d, &debug_info, boun_v, u_prev[0], u_prev[1], zero_slopes, zero_slopes, slopes[1], dt, x[1]-x[0], u_next[0], i, M1D + params1d.number_of_scalars, true, params1d.number_of_scalars, time_mom.curr_t, &beta[i], S[i], S[i], S[i+1], i );
            }
            else if ( i == params1d.cells_number - 1 ) {
                // обработка правого граничного условия
                current_inflow_parameters(&paramsc, &params1d, &empty_structure, RIGHT_BOUNDARY, curr_inflow_params); 
                boundary( &params1d, &empty_structure, u_prev[params1d.cells_number-1], params1d.right_bc, boun_v, time_mom.curr_t, curr_inflow_params );

                calc_step_in_cell( &paramsc, &params1d,&debug_info, u_prev[params1d.cells_number-2], u_prev[params1d.cells_number-1], boun_v,
                    slopes[params1d.cells_number-2], zero_slopes, zero_slopes, dt, x[params1d.cells_number] - x[params1d.cells_number-1], u_next[params1d.cells_number-1], i, M1D + params1d.number_of_scalars, true, params1d.number_of_scalars, time_mom.curr_t, &beta[i], S[i], S[i-1], S[i], i );

            }
            else {
                // расчет параметров во внутренней ячейке
                calc_step_in_cell( &paramsc, &params1d, &debug_info, u_prev[i-1], u_prev[i], u_prev[i+1], slopes[i-1], slopes[i], slopes[i+1], dt, x[i+1] - x[i], u_next[i], i, M1D + params1d.number_of_scalars, true, params1d.number_of_scalars, time_mom.curr_t, &beta[i], S[i] , S[i-1], S[i+1], i);
            }


            if ( paramsc.is_physics ) {
                // учет "физики" задачи	

		if( initial_ignition[i] == 1 && time_mom.curr_t < params1d.ignition_time){ // Условие поджига
		    initial_ignition[i] = 1;
		}
		else{
			if (initial_ignition[i] == 2)
			{
			}
			else
				initial_ignition[i] = 0;
		}


                physics_solver( &paramsc, &params1d, &empty_structure, X_DIRECTION, dt, u_next[i], u_next[i+1], u_next[i-1], number_of_block, &beta[i], &B[i], &C[i], M1D + params1d.number_of_scalars, ignition_flag[i], initial_ignition[i], S[i], time_mom.curr_t );
            }

        } // конец цикла по ячейкам

        // запись информации в файлы статистики на каждом шаге
        write_statistics( &paramsc, &params1d, &time_mom, &debug_info, u_prev );

        // Запись в файл информации о параметрах в ячейке 
        write_time_dependent_information_in_a_cell(&paramsc, &params1d, time_mom.curr_t, argv[2], "last_cell_infromation.dat", u_prev[0], M1D + params1d.number_of_scalars);

        for ( int i = 0; i < params1d.cells_number; i++ )
            for ( int j = 0; j < M1D + params1d.number_of_scalars; j++ )
                u_prev[i][j] = u_next[i][j];

        if (time_mom.steps_num == 188){
            printf("time_mom.steps_num == 188");
            write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, files_counter + 1 ,time_mom.curr_t, output_filename,  beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );
        }

        
        time_mom.curr_t += dt;
        time_mom.steps_num++;

        // проверка сохранения массы в системе
        check_mass( &params1d, x, u_next, initial_total_mass, &mass_diff );

        // печать информационного сообщения в конце шага
        printf( "\nStep %d, time %f: dt = %.2e\n", time_mom.steps_num - 1, time_mom.curr_t - dt, dt );
        printf( "Pressure relaxation info - case 1: %d  case 2: %d  case 3: %d  case 4: %d  case 5: %d\n", debug_info.relaxation_cases[0], debug_info.relaxation_cases[1],
        debug_info.relaxation_cases[2], debug_info.relaxation_cases[3], debug_info.relaxation_cases[4] );
        printf( "Godunov algorithm info - case 1: %d  case 2: %d  case 3: %d\n", debug_info.godunov_cases[0], debug_info.godunov_cases[1], debug_info.godunov_cases[2] );

        // запись промежуточных результатов
        if ( paramsc.is_output_on_time  ) {
            // вывод по моментам времени
            if ( time_mom.curr_t >= ( files_counter + 1 ) * time_gap ) {
                write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, files_counter + 1,time_mom.curr_t, output_filename, beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );
		// запись показаний датчиков
		if ( paramsc.are_sensors ) {
		    write_sensors( &paramsc, u_prev, &time_mom, sensors_files, sensors_cells, M1D + params1d.number_of_scalars );
		}

                files_counter++;
            }
            sprintf_s( tmp_str, "%f", ( files_counter + 1 ) * time_gap );
            write_restart_info( argv[2], tmp_str, output_filename );
        }
        else {
            // вывод по количеству шагов по времени
            if ( time_mom.steps_num >= ( files_counter + 1 ) * steps_gap ) {
                write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, files_counter + 1 ,time_mom.curr_t, output_filename,  beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );
                files_counter++;
            }
            sprintf_s( tmp_str, "%d", ( files_counter + 1 ) * steps_gap );
            write_restart_info( argv[2], tmp_str, output_filename );
        }


    } // конец цикла по времени

    // освобождение памяти
    free( xc );
    free( x );
    for ( int i = 0; i < params1d.cells_number; i++ ) {
        free( u_prev[i] );
        free( u_next[i] );
    }
    for ( int i = 0; i < params1d.cells_number; i++ ) {
        free( slopes[i] );
    }
	free( B );
	free( beta );
	free( ignition_flag );
	free( initial_ignition);
        free( S );
    // закрытие файлов с показаниями датчиков
    if ( paramsc.are_sensors ) {
        for ( int i_sensor = 0; i_sensor < paramsc.sensors_num; i_sensor++ )
            fclose( sensors_files[i_sensor] );
    }
    // закрытие файлов для записи статистики
    fclose( debug_info.relaxation_out );
    fclose( debug_info.godunov_out );

    printf( "\nNormal finish: total time steps = %d, total time = %f\n\n", time_mom.steps_num, time_mom.curr_t );
    
    return 0;

}