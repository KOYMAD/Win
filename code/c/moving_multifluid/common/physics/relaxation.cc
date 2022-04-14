// relaxation.cc
// Релаксация скоростей и давлений фаз на межфазной границе

// Релаксация скоростей реализована по статье:
// Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467.

// Реализация давлений реализована по статье:
// Иванов И.Э. Численное моделирование многофазных течений с большим содержанием дисперсной фазы // Вестник МАИ. – 2009. – Т. 16, № 2. – С. 62 – 70.
// Алгоритм детально описан в \science\utkin\docs\Релаксация давления.docx
// Скорее всего, все восходит к работе:
// Saurel R., Lemetayer O. A multiphase model for compressible flows with interfaces, shocks,
// detonation waves and cavitation // Journal of Fluid Mechanics. – 2001. – V. 431. – P. 239 – 271.

// Возможные места изменений: - в работе (Saurel R., Lemetayer O.) используется другая формула для pi;
//                            - не рассматривать отдельно случае близкого к нулю дискриминанта;
//                            - возможна ли ситуация, когда оба корня положительные, и оба подходят?
//                            - переписать процедуру релаксации давлений, вычислять сразу два набора кандидатов по единым формулам, и отбирать один;
//                            - порисовать параболы для типичных параметров задач. 

// (c) Уткин Павел, 2017
// Создан: 29 марта 2017 г.

#include "relaxation.h"

// Релаксационный оператор для решения системы уравнений типа Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия релаксационного оператора
// n - реальный размер вектора
// curr_time - текущий момент времени
// configuration_pressure - конфигурационное давление
void Lr( const struct ParametersCommon* paramsс, const struct Parameters1d* params1d, struct DebugInfo *debug_info, const double left_ncons[M], const double center_ncons[M],
         const double right_ncons[M], const double dt, const double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure) {
    if ( paramsс->velocity_relaxation == true ) {
        double solution_ncons_Lrv[M]; // вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора и оператора релаксации скорости
        Lrv( paramsс, center_ncons, dt, h, solution_ncons_Lrv, n ); // действие оператора релаксации скорости
        if (!paramsс->pressure_relaxation_compaction)
            Lrp_Ivanov( paramsс, debug_info, solution_ncons_Lrv, solution_ncons, n ); // действие оператора релаксации давления
        else
            Lrp_compaction(paramsс, params1d, debug_info, solution_ncons_Lrv, solution_ncons, step_number, n, configuration_pressure);
    }
    else{
        if (!paramsс->pressure_relaxation_compaction)
            Lrp_Ivanov( paramsс, debug_info, center_ncons, solution_ncons, n ); // действие оператора релаксации давления
        else{
            if (step_number == 11200)
                write_time_dependent_information_in_a_cell(paramsс, params1d, curr_time, "relaxation_information", "before_pressure_relaxation_center.dat", solution_ncons, n);
            Lrp_compaction(paramsс, params1d, debug_info, center_ncons, solution_ncons, step_number, n, configuration_pressure);
            if (step_number == 11200)
                write_time_dependent_information_in_a_cell(paramsс, params1d, curr_time, "relaxation_information", "after_pressure_relaxation_center.dat", solution_ncons, n);
        }
    }
}

// Оператор релаксации скорости для решения системы уравнений типа Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия релаксационного оператора скорости
// n - реальный размер векторов
void Lrv( const struct ParametersCommon* paramsc, const double center_ncons[M], const double dt, const double h, double solution_ncons[M], int n ) {

    // инициализируем решение текущим вектором в рассчитываемой ячейке
    for ( int i = 0; i < n; i++ )
        solution_ncons[i] = center_ncons[i];

    // выравнивание скоростей фаз
    double t1 = center_ncons[B_DISP] * center_ncons[R_DISP];
    double t2 = ( 1.0 - center_ncons[B_DISP] ) * center_ncons[R_GAS];
    solution_ncons[V_DISP] = ( t1 * center_ncons[V_DISP] + t2 * center_ncons[V_GAS] ) / ( t1 + t2 );
    solution_ncons[V_GAS] = solution_ncons[V_DISP];

    // пересчет внутренней энергии фаз
    double eg_0 = e_gas( paramsc, center_ncons[P_GAS], center_ncons[R_GAS] ); // внутренняя энергия газа до коррекции
    double ed_0 = e_disp( paramsc, center_ncons[P_DISP], center_ncons[R_DISP] ); // внутренняя энергия дисперсной фазы до коррекции
    double u_i = calc_u_i( paramsc, center_ncons ); // усредненная межфазная скорость до коррекции
    // внутренняя энергия газа после коррекции
    double eg_1 = eg_0 + 0.5 * ( solution_ncons[V_GAS] - center_ncons[V_GAS] ) * ( u_i - center_ncons[V_GAS] );
    // внутренняя энергия дисперсной фазы после коррекции
    double ed_1 = ed_0 - 0.5 * ( solution_ncons[V_DISP] - center_ncons[V_DISP] ) * ( u_i - center_ncons[V_DISP] );

    // пересчет давления по измененной внутренней энергии
    solution_ncons[P_GAS] = p_gas( paramsc, eg_1, center_ncons[R_GAS] );
    solution_ncons[P_DISP] = p_disp( paramsc, ed_1, center_ncons[R_DISP] );

}

// Оператор релаксации давления для решения системы уравнений типа Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора и оператора релаксации скорости
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического и полного релаксационного операторов
// n - реальный размер векторов
void Lrp_Ivanov( const struct ParametersCommon* paramsc, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int n ) {

    // инициализируем решение текущим вектором в рассчитываемой ячейке
    for ( int i = 0; i < n; i++ )
        solution_ncons[i] = center_ncons[i];

    double pi = center_ncons[B_DISP] * center_ncons[P_DISP] + ( 1.0 - center_ncons[B_DISP] ) * center_ncons[P_GAS]; // давление на межфазной границе на текущем слое

    double g1 = paramsc->g1;
    double g2 = paramsc->g2;
    double p01 = paramsc->p01;
    double p02 = paramsc->p02;

    double g2m = g2 - 1.0; // показатель адиабаты газа без единицы
    double g1m = g1 - 1.0; // показатель адиабаты дисперсной фазы без единицы
    
    // расчет коэффициентов квадратного уравнения для плотности газовой фазы на следующем шаге по времени

    double a1 = 0.5 * center_ncons[R_GAS] * ( g2 + 1.0 );
    double b1 = - ( center_ncons[P_GAS] + p02 * paramsc->g2 ) - 0.5 * g2m * pi;
    double c1 = - 0.5 * g2m;
    double d1 = - p02 * g2 * center_ncons[R_GAS] - 0.5 * g2m * pi * center_ncons[R_GAS];

    double a2 = 0.5 * center_ncons[R_DISP] * ( g1 + 1.0 );
    double b2 = - ( center_ncons[P_DISP] + p01 * paramsc->g1 ) - 0.5 * g1m * pi;
    double c2 = - 0.5 * g1m;
    double d2 = - p01 * g1 * center_ncons[R_DISP] - 0.5 * g1m * pi * center_ncons[R_DISP];

    double a3 = center_ncons[B_DISP] * center_ncons[R_DISP];
    double b3 = ( 1.0 - center_ncons[B_DISP] ) * center_ncons[R_GAS];
    double c3 = - 1.0;

    // на всякий случай, можно вставить проверки на то, что коэффициенты ненулевые

    double G =  - a2 * b1 * c3 - a3 * b2 * c1 + a3 * b1 *c2 - c1 * c3 * d2;
    double H =  - a2 * b1 * b3 + a2 * c3 * d1 - a1 * a3 * b2 - a3 * c2 * d1 - b3 * c1 * d2 - a1 * c3 * d2;
    double I = a2 * b3 * d1 - a1 * b3 * d2;
    
    // решение уравнения G * y^2 + H * y + I = 0 и анализ корней

    if ( fabs( G ) < paramsc->eps_general ) {
        // уравнение на самом деле не является квадратным
        if ( fabs( H ) < paramsc->eps_general ) {
            printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
            vector_debug_print( debug_info );
            exit( EXIT_FAILURE );
        }
        else {
            solution_ncons[R_GAS] = - I / H;
            solution_ncons[P_DISP] = ( d1 - b1 * solution_ncons[R_GAS] ) / ( a1 + c1 * solution_ncons[R_GAS] );
            if ( solution_ncons[P_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[P_GAS] = solution_ncons[P_DISP];
            solution_ncons[R_DISP] = - a3 * solution_ncons[R_GAS] / ( b3 + c3 * solution_ncons[R_GAS] );
            if ( solution_ncons[R_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
            if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
        }
        debug_info->relaxation_cases[0]++;
    }
    else {
        double D = H * H - 4.0 * G * I;
        if ( D > 0 ) { // вариант - D > epsilon
            // корни квадратного уравнения
            double y1 = 0.5 * ( - H + sqrt( D ) ) / G;
            double y2 = 0.5 * ( - H - sqrt( D ) ) / G;
            bool first_root = false; // первый корень на данный момент решение задачи?
            bool second_root = false; // второй корень на данный момент решение задачи?
            if ( y1 < 0 ) {
                // первый корень заведомо не подходит
                if ( y2 < 0 ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
                else {
                    // рассчитываем все по второму корню
                    second_root = true;
                    solution_ncons[R_GAS] = y2;
                    solution_ncons[P_DISP] = ( d1 - b1 * y2 ) / ( a1 + c1 * y2 );
                    if ( solution_ncons[P_DISP] < 0.0 ) {
                        printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                        vector_debug_print( debug_info );
                        exit( EXIT_FAILURE );
                    }
                    solution_ncons[P_GAS] = solution_ncons[P_DISP];
                    solution_ncons[R_DISP] = - a3 * y2 / ( b3 + c3 * y2 );
                    if ( solution_ncons[R_DISP] < 0.0 ) {
                        printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                        vector_debug_print( debug_info );
                        exit( EXIT_FAILURE );
                    }
                    solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
                    if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                        printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                        vector_debug_print( debug_info );
                        exit( EXIT_FAILURE );
                    }
                    debug_info->relaxation_cases[1]++;
                }
            }
            else {
                // пытаемся рассчитать все по первому коорню
                first_root = true;
                solution_ncons[R_GAS] = y1;
                solution_ncons[P_DISP] = ( d1 - b1 * y1 ) / ( a1 + c1 * y1 );
                if ( solution_ncons[P_DISP] < 0.0 ) {
                    // первый корень нефизичен по найденному давлению на межфазной границе, нужно пробовать второй
                    first_root = false;
                }
                if ( first_root == true )
                    solution_ncons[P_GAS] = solution_ncons[P_DISP];
                if ( first_root == true ) {
                    solution_ncons[R_DISP] = - a3 * y1 / ( b3 + c3 * y1 );
                    if ( solution_ncons[R_DISP] < 0.0 ) {
                        // первый корень нефизичен по найденной плотности дисперсной фазы, нужно пробовать второй
                        first_root = false;
                    }
                }
                if ( first_root == true ) {
                    solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
                    if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                        // первый корень нефизичен по найденной объемной доле дисперсной фазы, нужно пробовать второй
                        first_root = false;
                    }
                    debug_info->relaxation_cases[2]++;
                }
            }
            if ( first_root == false && second_root == false ) {
                // первый корень не подошел в результате анализа, пробуем второй
                // нужна проверка y2 < 0 - если так, то ошибка
                solution_ncons[R_GAS] = y2;
                solution_ncons[P_DISP] = ( d1 - b1 * y2 ) / ( a1 + c1 * y2 );
                if ( solution_ncons[P_DISP] < 0.0 ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 4\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
                solution_ncons[P_GAS] = solution_ncons[P_DISP];
                solution_ncons[R_DISP] = - a3 * y2 / ( b3 + c3 * y2 );
                if ( solution_ncons[R_DISP] < 0.0 ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 4\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
                solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
                if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 4\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
            }
            debug_info->relaxation_cases[3]++;
        }
        else if ( fabs( D ) < paramsc->eps_general ) {
            // один корень
            solution_ncons[R_GAS] = - 0.5 * H / G;
            solution_ncons[P_DISP] = ( d1 - b1 * solution_ncons[R_GAS] ) / ( a1 + c1 * solution_ncons[R_GAS] );
            if ( solution_ncons[P_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 5\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[P_GAS] = solution_ncons[P_DISP];
            solution_ncons[R_DISP] = - a3 * solution_ncons[R_GAS] / ( b3 + c3 * solution_ncons[R_GAS] );
            if ( solution_ncons[R_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 5\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
            if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 6\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            debug_info->relaxation_cases[4]++;
        }
        else {
            printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, no roots\n", debug_info->current_cell );
            vector_debug_print( debug_info );
            exit( EXIT_FAILURE );
        }
    }

}

// Оператор релаксации давления для решения системы уравнений типа Baer-Nunziato с учетом компактирования Schwendeman и Saurel
// paramsc - структура с основными параметрами вычислительного эксперимента
// params1d - структура c параметрами вычислительного эксперимента, присущими 1d случаю
// debug_info - структура с отладочной информацией
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора и оператора релаксации скорости
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического и полного релаксационного операторов
// step_number - номер текущей ячейки
void Lrp_compaction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int step_number, int n, double *configuration_pressure){

    // инициализируем решение текущим вектором в рассчитываемой ячейке
    for ( int i = 0; i < n; i++ )
        solution_ncons[i] = center_ncons[i];    

    double g1 = paramsc->g1;
    double g2 = paramsc->g2;
    double p01 = paramsc->p01;
    double p02 = paramsc->p02;

    double z; // плотность дисперсной фазы - искомая переменная
    double x, y; // давление газовой фазы и плотность газовой фазы

    double dx, dy; // производные x и y по z

    double beta; // конфигурационное давление
    double f, df; // конфигурационное давление бета и проивзодная по z. Вычисляется на каждом шаге

    double A, a, b; // вспомогательные коэффициенты

    a = solution_ncons[B_DISP] * solution_ncons[R_DISP];
    b = (1.0 - solution_ncons[B_DISP]) * solution_ncons[R_GAS];

    double common_part; // общая часть (B / tau)^(1/n) в формуле SAUREL_COMPACTION

    double p0, p0_disp, b0, b0_disp, r0_disp;

    // определение блока, в котором сейчас находимся
    int number_of_block = 0;

    if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION){
        for (int i = 0; i < params1d->ic_blocks_number; i++){
	    if (step_number <= params1d->cell_end[i] && step_number >= params1d->cell_begin[i])
	        number_of_block = i;
        }

        p0 = params1d->block_values[number_of_block][P_GAS];
        p0_disp = params1d->block_values[number_of_block][P_DISP];
        b0 = 1.0 - params1d->block_values[number_of_block][B_DISP];
        b0_disp = params1d->block_values[number_of_block][B_DISP];
        r0_disp = params1d->block_values[number_of_block][R_DISP];
  
        if (solution_ncons[B_DISP] < paramsc->volume_fraction_compaction)
            A = 0.0;
        else
            A = -(p0 - p0_disp) / b0_disp / r0_disp * (2.0 - b0_disp) * (2.0 - b0_disp) / log(1.0 - b0_disp); // константа в выражении для конфигурационного давления

        beta = A * a * log(1.0 - solution_ncons[B_DISP]) / pow(2.0 - solution_ncons[B_DISP], 2.0);
    }
    else{
        
        common_part = (1.0 - solution_ncons[B_DISP]) * log( (1.0 - solution_ncons[B_DISP]) / (1.0 - paramsc->volume_fraction_compaction) ) + solution_ncons[B_DISP] - paramsc->volume_fraction_compaction;
        if (common_part < 0.0 && fabs(common_part) < paramsc->eps_relax_compaction){
            common_part = 0.0;
        }
        if (solution_ncons[B_DISP] > paramsc->volume_fraction_compaction)
            beta = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - solution_ncons[B_DISP]) / (1.0 - paramsc->volume_fraction_compaction) ) * pow (common_part , paramsc->n_parameter - 1.0);
        else
            beta = 0.0;
    }

    double fi, dfi;

    z = solution_ncons[R_DISP];
 //   printf("Before relaxation procedure\n");
//    for (int i = 0; i < n; i++)
//        printf("solution_ncons[%d] = %lf \n", i, solution_ncons[i]);
//    printf("\n");
 //   if (solution_ncons[P_DISP] - beta > solution_ncons[P_GAS]){

    do{
		p0 = params1d->block_values[number_of_block][P_GAS];
		p0_disp = params1d->block_values[number_of_block][P_DISP];
		b0 = 1.0 - params1d->block_values[number_of_block][B_DISP];
		b0_disp = params1d->block_values[number_of_block][B_DISP];
		r0_disp = params1d->block_values[number_of_block][R_DISP];
        if (a / z < paramsc->volume_fraction_compaction)
            A = 0.0;
        else
            A = -(p0 - p0_disp) / b0_disp / r0_disp * (2.0 - b0_disp) * (2.0 - b0_disp) / log(1.0 - b0_disp); // константа в выражении для конфигурационного давления

        y = b * z / (z - a);
        x = ( 1.0 / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) + (solution_ncons[P_GAS] + g2 * p02) / (g2 - 1.0) / solution_ncons[R_GAS] + solution_ncons[P_GAS] / 2.0 / solution_ncons[R_GAS] ) / 
            ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 0.5 / solution_ncons[R_GAS] );
        
        if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION){
            f = A * a * log(1.0 - a / z) / pow ( 2.0 - a / z , 2.0 ); 
            df = A * a * a / z / z / (2 - a / z) / (2 - a / z) * ( 1.0 / (1.0 - a / z) - 2.0 * log (1.0 - a / z) / (2.0 - a / z) );// (2 - a/z)
        }
        else{
            if (a / z > paramsc->volume_fraction_compaction){
                common_part = ((1.0 - a / z) * log ( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ) + a / z - paramsc->volume_fraction_compaction);
                if (common_part < 0.0 && fabs(common_part) < paramsc->eps_relax_compaction){
                    common_part = 0.0;
                }
                f = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ) * pow ( common_part , paramsc->n_parameter - 1.0);
                df = - a * a * paramsc->tau_parameter * paramsc->n_parameter / z / z * ( z / (z - a) * pow(common_part, paramsc->n_parameter - 1.0) + 
                    (paramsc->n_parameter - 1.0) * pow(log( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ), 2.0) * pow(common_part, paramsc->n_parameter - 2.0) );

                if (f == 0.0)
                    df = 0.0;

            }
            else{
                f = 0;
                df = 0;
            }
        }
        fi = (x + f + g1 * p01) / (g1 - 1.0) / z - (solution_ncons[P_DISP] + g1 * p01) / (g1 - 1.0) / solution_ncons[R_DISP] + 0.5 * (x + f + solution_ncons[P_GAS] + beta) * (1.0 / z - 1.0 / solution_ncons[R_DISP]);

        dy = - b * a / (z - a) / (z - a);
        dx = (- dy / y / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) * ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 1.0 / 2.0 / solution_ncons[R_GAS] ) + 
            (1.0 / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) + (solution_ncons[P_GAS] + g2 * p02) / (g2 - 1.0) / solution_ncons[R_GAS] + solution_ncons[P_GAS] / 2.0 / solution_ncons[R_GAS] ) * (- dy / y / y * (g2 + 1.0) / 2.0 / (g2 - 1.0) ) ) /
            pow ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 1.0 / 2.0 / solution_ncons[R_GAS] , 2.0 );
        
        dfi = (dx + df) / (g1 - 1.0) / z - 1.0 / z / z * (x + f + g1 * p01) / (g1 - 1.0) - 1.0 / 2.0 / z / z * (x + f + solution_ncons[P_GAS] + beta) + 0.5 * (dx + df) * (1.0 / z - 1.0 / solution_ncons[R_DISP]);

        z = z - fi / dfi;

    }while(fabs(fi/dfi / solution_ncons[R_DISP]) > 10 * paramsc->eps_general);

    y = b * z / (z - a);
    x = ( 1.0 / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) + (solution_ncons[P_GAS] + g2 * p02) / (g2 - 1.0) / solution_ncons[R_GAS] + solution_ncons[P_GAS] / 2.0 / solution_ncons[R_GAS] ) / 
            ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 0.5 / solution_ncons[R_GAS] );

    double v_ncons_new[M];
    double beta_new;
    v_ncons_new[R_DISP] = z;
    v_ncons_new[R_GAS] = y;
    v_ncons_new[P_GAS] = x;
    v_ncons_new[B_DISP] = solution_ncons[B_DISP] * solution_ncons[R_DISP] / v_ncons_new[R_DISP];
    if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION)
        beta_new = A * a * log(1.0 - v_ncons_new[B_DISP]) / pow(2.0 - v_ncons_new[B_DISP], 2.0);
    else{
        if (v_ncons_new[B_DISP] > paramsc->volume_fraction_compaction)
            beta_new = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - a / v_ncons_new[R_DISP]) / (1.0 - paramsc->volume_fraction_compaction) ) * pow ( common_part , paramsc->n_parameter - 1.0);
        else
            beta_new = 0.0;
    }
    v_ncons_new[P_DISP] = beta_new + v_ncons_new[P_GAS];
    v_ncons_new[V_DISP] = solution_ncons[V_DISP];
    v_ncons_new[V_GAS] = solution_ncons[V_GAS];
    
   /* double check_f;
    if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION)
        check_f = A * a * log(1.0 - a / z) / pow ( 2.0 - a / z , 2.0 );
    else
        if (solution_ncons[B_DISP] > paramsc->volume_fraction_compaction )
            check_f = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ) * pow ( common_part , paramsc->n_parameter - 1.0);
        else
            check_f = 0.0;
    double check_fi = (x + check_f + g1 * p01) / (g1 - 1.0) / z - (solution_ncons[P_DISP] + g1 * p01) / (g1 - 1.0) / solution_ncons[R_DISP] + 0.5 * (x + check_f + solution_ncons[P_GAS] + beta) * (1.0 / z - 1.0 / solution_ncons[R_DISP]);
*/
   // if (v_ncons_new[P_DISP] <= 0.0 )
        //printf("current step is %d, check_fi = %lf, v_ncons_new[P_DISP] = %lf\n", step_number, check_fi, v_ncons_new[P_DISP]);

   // if (v_ncons_new[P_GAS] <= 0.0 )
        //printf("current step is %d, check_fi = %lf, v_ncons_new[P_DISP] = %lf\n", step_number, check_fi, v_ncons_new[P_GAS]);

    //if (check_fi > 1.e-5)
        //printf("current step is %d, check_fi = %lf\n", step_number, check_fi);

    solution_ncons[B_DISP] = v_ncons_new[B_DISP];          
    solution_ncons[R_DISP] = v_ncons_new[R_DISP];  
    solution_ncons[V_DISP] = v_ncons_new[V_DISP];  
    solution_ncons[P_DISP] = v_ncons_new[P_DISP];  
    solution_ncons[R_GAS] = v_ncons_new[R_GAS];  
    solution_ncons[V_GAS] = v_ncons_new[V_GAS];  
    solution_ncons[P_GAS] = v_ncons_new[P_GAS];  
    *configuration_pressure = beta_new;

  //  printf("After relaxation procedure\n");
  //  for (int i = 0; i < n; i++)
 //       printf("solution_ncons[%d] = %lf \n", i, solution_ncons[i]);
 //   printf("\n");
}