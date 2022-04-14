// godunov.cc
// Метод Годунова расчета потоков для одномерной системы уравнений Баера-Нунциато.
// Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
// two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526.
// (c) Уткин Павел, 2013
// Создан: 4 июля 2013 г.

#include "godunov_bn.h"

// Метод Годунова интегрирования одномерной системы уравнений Баера-Нунциато,
// расчет одного шага по времени в одной ячейке
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_params[M] - вектор примитивных переменных в ячейке слева от рассчитываемой (in)
// center_params[M] - вектор примитивных переменных в рассчитываемой ячейке (in)
// right_params[M] - вектор примитивных переменных в ячейке справа от рассчитываемой (in)
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой (in)
// slopes_center - вектор наклонов в рассчитываемой ячейке (in)
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге (out)
// n - реальный размер векторов
// configuration_pressure - конфигурационное давление
void godunov( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_params[M], double center_params[M],
              double right_params[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
              double dt, double h, double solution[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, double curr_time, double *configuration_pressure ) {

    double left_godunov_flux[M];            /* полный "поток" через левую грань ячейки */
    double right_godunov_flux[M];           /* полный "поток" через правую грань ячейки */
    double center_cons_params[M];           /* вектор консервативных переменных в рассчитываемой ячейке на текущем шаге */
    double solution_cons[M];                /* вектор консервативных переменных в рассчитываемой ячейке на новом шаге */
    ReturnCodes left_flux_code, right_flux_code, pressure_code; /* код возвратов функций расчета потоков и перевода
                                                                   "консервативных переменных" в примитивные */
    Disp_phase_cases solver_part_left;      /* константа, которая определяет, какую из четырех частей солвера использовать
                                               для расчета "потока" через левую грань ячейки */
    Disp_phase_cases solver_part_right;     /* константа, которая определяет, какую из четырех частей солвера использовать
                                               для расчета "потока" через правую грань ячейки */

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
        convert_noncons_to_cons( paramsc, left_params, left_minus_cons, 0 );
        left_minus_cons[j] += 0.5 * h * slopes_left[j];
        convert_cons_to_noncons( paramsc, left_minus_cons, left_minus_ncons, 0 );
        // реконструкция значений на левой границе справа
        convert_noncons_to_cons( paramsc, center_params, left_plus_cons, 0 );
        left_plus_cons[j] -= 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, left_plus_cons, left_plus_ncons, 0 );
        // реконструкция значений на правой границе слева
        convert_noncons_to_cons( paramsc, center_params, right_minus_cons, 0 );
        right_minus_cons[j] += 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, right_minus_cons, right_minus_ncons, 0 );
        // реконструкция значений на правой границе справа
        convert_noncons_to_cons( paramsc, right_params, right_plus_cons, 0 );
        right_plus_cons[j] -= 0.5 * h * slopes_right[j];
        convert_cons_to_noncons( paramsc, right_plus_cons, right_plus_ncons, 0 );
    }

    // анализ соотношения между величинами объемной доли дисперсной фазы слева от текущей ячейки и в текущей
    solver_part_left = what_case( paramsc, left_minus_ncons[B_DISP], left_plus_ncons[B_DISP] );
    
    // анализ соотношения между величинами объемной доли дисперсной фазы в текущей ячейке и справа от нее
    solver_part_right = what_case( paramsc, right_minus_ncons[B_DISP], right_plus_ncons[B_DISP] );

    if ( solver_part_left == NO_GRAD && solver_part_right == NO_GRAD ) {
        // полное расщепление фаз - находим решение с использованием сокращенного вектора
        full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution, n );
        debug_info->godunov_cases[0]++;
        debug_info->godunov_cases[0]++;
    }
    else { 
        // используем полный вектор, даже если содержание объемной фазы в ячейке мало
        
        // расчет полного "потока" через левую грань
        
        // заполнение части полей структуры для отладки кода
        debug_info->neighbour_cell = debug_info->current_cell - 1;
        for ( int i = 0; i < n; i++ )
            debug_info->current_cell_vncons[i] = left_plus_ncons[i];
        for ( int i = 0; i < n; i++ )
            debug_info->neighbour_cell_vncons[i] = left_minus_ncons[i];
        
        left_flux_code = godunov_flux( paramsc, debug_info, left_minus_ncons, left_plus_ncons,
            solver_part_left, LEFT, center_params[B_DISP], left_godunov_flux, n );
        if ( left_flux_code == GODUNOV_FAILS ) {
            // метод Годунова не сработал - пробуем полностью расщепить и найти решение
            full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution, n );
            return;
        }

        // расчет полного "потока" через правую грань
        
        // заполнение части полей структуры для отладки кода, касающихся данных о соседней ячейке
        debug_info->neighbour_cell = debug_info->current_cell + 1;
        for ( int i = 0; i < n; i++ )
            debug_info->current_cell_vncons[i] = right_minus_ncons[i];
        for ( int i = 0; i < n; i++ )
            debug_info->neighbour_cell_vncons[i] = right_plus_ncons[i];
        
        right_flux_code = godunov_flux( paramsc, debug_info, right_minus_ncons, right_plus_ncons,
            solver_part_right, RIGHT, center_params[B_DISP], right_godunov_flux, n );
        if ( right_flux_code == GODUNOV_FAILS ) {
            // метод Годунова не сработал - пробуем полностью расщепить и найти решение
            full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution , n);
            return;
        }

        // обновление параметров в ячейке
        convert_noncons_to_cons( paramsc, center_params, center_cons_params, 0 );
        for ( int i = 0; i < n; i++ ) {
            solution_cons[i] = center_cons_params[i] - dt * ( right_godunov_flux[i] - left_godunov_flux[i] ) / h;
        }
        pressure_code = convert_cons_to_noncons( paramsc, solution_cons, solution, 0 );
        if ( pressure_code == NEGATIVE_PRESSURE ) {
            /* в результате расчета получилось отрицательное давление в дисперсной фазе - пробуем полностью расщепить */
            full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution, n );
        }
    }

    // при необходимости делаем релаксацию скоростей и давлений
    double solution_ncons_Lh[M]; // вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
    for ( int i = 0; i < n; i++ )
        solution_ncons_Lh[i] = solution[i];
    if ( paramsc->pressure_relaxation == true && is_pressure_relaxation_after_this_step == true )
        Lr( paramsc, params1d, debug_info, left_params, solution_ncons_Lh, right_params, dt, h, solution, step_number, n, curr_time, configuration_pressure ); // действие релаксационного оператора

}

// Анализ перепада значений объемной доли дисперсной фазы слева и справа от разрыва
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// left_beta - величина объемной доли дисперсной фазы слева от разрыва (in)
// right_beta - величина объемной доли дисперсной фазы справа от разрыва (in)
// Возвращает одну из констант перечисления Disp_phase_cases
Disp_phase_cases what_case( struct ParametersCommon *paramsc, double left_beta, double right_beta ) {
    
    if ( fabs( left_beta - right_beta ) <= paramsc->eps_decouple ) {
        // перепад значений объемной доли слева и справа от разрыва незначительный - система уравнений расщепляется,
        // в уравнениях пропадает объемная доля дисперсной фазы; предполагается, что данный случай включает ситуацию
        // отсутствия дисперсной фазы по обе стороны от разрыва, то есть константа params->eps_decouple значительно
        // больше константы params->eps_disp_abs
        return NO_GRAD;
    }
    else { 
        // перепад значений объемной доли слева и справа от разрыва значительный - рассматриваем уравнения с неконсервативными
        // членами в правой части
        if ( left_beta > paramsc->eps_disp_abs ) {
            // дисперсная фаза слева от разрыва присутствует
            if ( right_beta > paramsc->eps_disp_abs ) {
                // дисперсная фаза справа от разрыва присутствует
                return BOTH_GRAD;
            }
            else {
                // дисперсная фаза слева от разрыва присутствует, а справа отсутствует
                return LEFT_ONLY_GRAD;
            }
        }
        else {
            // дисперсная фаза слева от разрыва отсутствует
            if ( right_beta > paramsc->eps_disp_abs ) {
                // дисперсная фаза справа от разрыва присутствует
                return RIGHT_ONLY_GRAD;
            }
            else {
                // дисперсная фаза отсутствует по обе стороны от разрыва - попадать сюда не должны, это случай
                // отсутствия перепада значений объемной доли слева и справа от разрыва, но подстраховываемся
                return NO_GRAD;
            }
        }
    }

}

// Получение решения на следующем шаге по времени для случая отсутствия неконсервативного члена в правой части - 
// полное расщепление фаз
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// center_params[M] - вектор примитивных переменных в рассчитываемой ячейке (in)
// left_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грании ячейки слева от рассчитываемой (in)
// left_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани рассчитываемой ячейки (in)
// right_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грани рассчитываемой ячейки (in)
// right_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани ячейки справа от рассчитываемой (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге (out)
// n (in)
void full_decouple_case_sol( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double center_params[M], double left_minus_ncons[M], 
                             double left_plus_ncons[M], double right_minus_ncons[M], double right_plus_ncons[M], double dt, double h, double solution[M], int n ) {

    double v_ncons_res_solid_left[M];
    double v_ncons_res_solid_right[M];

    // обновляем переменные вектора solution, соответствующие газовой фазе
    godunov_classical_one_phase( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, GAS_PHASE, dt, h, 
        v_ncons_res_solid_left, v_ncons_res_solid_right, solution, n );
    // обновляем переменные вектора solution, соответствующие дисперсной фазе
    godunov_classical_one_phase( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, DISPERSED_PHASE, dt, h, 
        v_ncons_res_solid_left, v_ncons_res_solid_right, solution, n );
    // объемную долю оставляем без изменения
    solution[B_DISP] = center_params[B_DISP];

}

// Если в ячейке фактически отсутствует дисперсная фаза, то установить в ней фоновые параметры
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке (in/out)
void set_background_state( struct ParametersCommon *paramsc, double solution[M] ) {

    if ( solution[B_DISP] < paramsc->eps_disp_abs ) {
        solution[R_DISP] = paramsc->background_density;
        solution[V_DISP] = paramsc->background_velocity;
        solution[P_DISP] = paramsc->background_pressure;
    }

}

// Расчет полного "потока" через грань ячейки методом Годунова
// params - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_params[M] - вектор примитивных переменных слева от разрыва (in)
// right_params[M] - вектор примитивных переменных справа от разрыва (in)
// solver_part - константа, которая определяет, какую из четырех частей солвера использовать для расчета "потока" (in)
// edge - идентификатор грани, через которую считается "поток" - LEFT или RIGHT (in)
// curr_cell_beta - значение объемной доли в рассчитываемой ячейке (in)
// godunov_flux[M] - полный "поток" Годунова (out)
// n - реальный размер вектора left_params, right_params и godunov_dlux (in)
// Возвращает: SUCCESS          в случае успеха
//             GODUNOV_FAILS    в случае невозможности построить решение задачи Римана
ReturnCodes godunov_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                          Disp_phase_cases solver_part, Direction edge, double curr_cell_beta, double godunov_flux[M], int n ) {

    double cons_flux[M]; // консервативная составляющая вектора "потока"
    double v_ncons_res[M]; // вектор-решение задачи о распаде разрыва, не используется, но нужен в качестве
                           // аргумента для godunov_cons_flux
    double cont_ncons[M]; //нужен для аргумента
    double s = 0.0; // значение автомодельной переменной x/t
    double ncons_flux[M]; // неконсервативная составляющая вектора "потока"
    double v_solid_cont; // скорость контактного разрыва в дисперсной фазе
    double p_solid_cont_l, p_solid_cont_r; // давления слева и справа от контактного разрыва в дисперсной фазе
    int i; // индекс компонента вектора
    ReturnCodes return_code; // код возврата функции нахождения параметров на распаде
    
    if ( solver_part == NO_GRAD ) {
        /* требуется расщепить систему уравнений и рассчитать сокращенные вектора потоков, но затем
           сформировать полный "поток", используя в качестве объемной доли значение в текущей ячейке */
        full_decouple_case_flux( paramsc, debug_info, left_params, right_params, curr_cell_beta, godunov_flux );
        debug_info->godunov_cases[0]++;
    }
    else {
        /* система уравнение не расщепляется, общий случай */

        /* расчет консервативной составляющей вектора "потока" */
        return_code = godunov_cons_flux( paramsc, debug_info, left_params, right_params, s, solver_part,
            cons_flux, v_ncons_res, &v_solid_cont, &p_solid_cont_l, &p_solid_cont_r, n, cont_ncons );
        if ( return_code == GODUNOV_FAILS ) {
            /* задачу Римана решить не удалось */
            return return_code;
        }

        /* расчет неконсервативной составляющей вектора "потока" */
        for ( i = 0; i < n; i++ )
            ncons_flux[i] = 0.0;   /* инициализация */
        if ( v_solid_cont > 0.0 && edge == LEFT )
            calc_ncons_term( v_solid_cont, left_params[B_DISP], right_params[B_DISP], p_solid_cont_l, p_solid_cont_r, ncons_flux );
        if ( v_solid_cont < 0.0 && edge == RIGHT )
            calc_ncons_term( v_solid_cont, left_params[B_DISP], right_params[B_DISP], p_solid_cont_l, p_solid_cont_r, ncons_flux );


        /* расчет полного "потока" */
        if ( edge == LEFT ) {
            for ( i = 0; i < n; i++ )
                godunov_flux[i] = cons_flux[i] + ncons_flux[i];
        }
        else {
            for ( i = 0; i < n; i++ )
                godunov_flux[i] = cons_flux[i] - ncons_flux[i];
        }
    }

    return SUCCESS;

}

/* Тебуется расщепить систему уравнений и рассчитать сокращенные вектора потоков, но затем
   сформировать полный "поток", используя в качестве объемной доли значение в текущей ячейке.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   left_params_full[M] - полный вектор примитивных переменных слева от разрыва (in)
   right_params_full[M] - полный вектор примитивных переменных справа от разрыва (in)
   beta - значение объемной доли в рассчитываемой ячейке (in)
   
   flux_full[M] - полный "поток" (out) */
void full_decouple_case_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params_full[M], double right_params_full[M],
                              double beta, double flux_full[M] ) {

    double flux_reduced[M_REDUCTION];           /* сокращенный однофазный вектор потока */
    double left_params_reduced[M_REDUCTION];    /* однофазный вектор примитивных переменных в ячейке слева от рассматриваемого разрыва */
    double right_params_reduced[M_REDUCTION];   /* однофазный вектор примитивных переменных в ячейке справа от рассматриваемого разрыва */
    double v_ncons_res[M_REDUCTION];            /* однофазный вектор примитивных переменных, полученный в результате решения задачи о распаде разрыва */
    double cont_red[M_REDUCTION];
    /* расчет вектора потока в газовой фазе */
    convert_full_to_reduced( left_params_full, GAS_PHASE, left_params_reduced );
    convert_full_to_reduced( right_params_full, GAS_PHASE, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, GAS_PHASE, v_ncons_res, flux_reduced,cont_red );
    convert_reduced_to_full( flux_reduced, GAS_PHASE, flux_full );

    /* расчет вектора потока в дисперсной фазе */
    convert_full_to_reduced( left_params_full, DISPERSED_PHASE, left_params_reduced );
    convert_full_to_reduced( right_params_full, DISPERSED_PHASE, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, DISPERSED_PHASE, v_ncons_res, flux_reduced,cont_red );
    convert_reduced_to_full( flux_reduced, DISPERSED_PHASE, flux_full );

    double beta_new = full_decouple_case_flux_volume_fraction( 0, left_params_full[B_DISP], beta, right_params_full[B_DISP], v_ncons_res[V]);

    if (beta_new != beta){
        printf("here\n");
    }

    /* учет объемной доли */
    flux_full[B_DISP] = 0.0;
    flux_full[R_DISP] *= beta_new;
    flux_full[V_DISP] *= beta_new;
    flux_full[P_DISP] *= beta_new;
    flux_full[R_GAS] *= 1.0 - beta_new;
    flux_full[V_GAS] *= 1.0 - beta_new;
    flux_full[P_GAS] *= 1.0 - beta_new;
}

/*  Функция расчета консервативной составляющей "потока" методом Годунова
    для одномерной системы уравнений Баера-Нунциато

    paramsc - структура с основными параметрами вычислительного эксперимента (in)
    debug_info - структура с отладочной информацией (in)
    left_params[M] - вектор примитивных переменных слева от разрыва (in)
    right_params[M] - вектор примитивных переменных справа от разрыва (in)
    s - значение автомодельной переменной x/t (in)
    solver_part - константа, которая определяет, какую из четырех частей солвера использовать для расчета "потока" (in)

    flux[M] - консервативная составляющая вектора "потока" (out)
    v_ncons_res[M] - вектор-решение задачи о распаде разрыва, нужен в выходных параметрах для построения точного решения (out)
    v_solid_cont - скорость контактного разрыва в дисперсной фазе (out)
    p_solid_cont_l - давление в дисперсной фазе слева от контактного разрыва в дисперсной фазе (out)
    p_solid_cont_r - давление в дисперсной фазе справа от контактного разрыва в дисперсной фазе (out)
    n - реальный размер вектора left_params, right_params и godunov_dlux (in)
    Возвращает: SUCCEESS        успешная сходимость итераций
                GODUNOV_FAILS   итерации не сошлись, решение задачи Римана построить не удалось */
ReturnCodes godunov_cons_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                               double s, Disp_phase_cases solver_part, double flux[M], double v_ncons_res[M], double *v_solid_cont,
                               double *p_solid_cont_l, double *p_solid_cont_r, int n, double cont_ncons[M] ) {

    double c1l, c2l;                    /* скорости звука в дисперсной и газовой фазах слева от разрыва */
    double c1r, c2r;                    /* скорости звука в дисперсной и газовой фазах справа от разрыва */
    double p_cont_disp, p_cont_gas;     /* давления на контактном разрыве в газовой и дисперсной фазах без учета разрыва
                                           объемной доли дисперсной фазы */
    double v_cont_disp, v_cont_gas;     /* скорости на контактном разрыве в газовой и дисперсной фазах без
                                           учета разрыва объемной доли дисперсной фазы */
    double solid_discontinuity_pressures[K_GENERAL_CASE];   /* [GAS_LEFT] - давление слева в газе, [GAS_RIGHT] - давление справа в газе,
                                                               [DISP_LEFT] - давление слева в дисперсной фазе,
                                                               [DISP_RIGHT] - давление справа в дисперсной фазе */
    double sound_velocities[K_GENERAL_CASE];    /* вектор со скоростями звука фаз по разные стороны от разрыва для
                                                   передачи в функции */
    double v1, v2;                      /* скорости контактных разрывов в дисперсной и газовой фазах */
    double v21, v22;                    /* скорости газа слева и справа от контактного разрыва дисперсной фазы */
    int i_comp;                         /* счетчик компонент вектора-решения */
    int return_code;                    /* код возврата функции нахождения параметров на распаде */
    
    /* 0. Рассчитываем скорости звука в обеих фазах слева и справа от разрыва для последующего использования */
    calc_sound_velocity( paramsc, left_params, &c1l, &c2l );
    calc_sound_velocity( paramsc, right_params, &c1r, &c2r );
    fill_sound_velocities( c2l, c2r, c1l, c1r, sound_velocities );

    /* 1. Построение начального приближения - давления и скорости на контактном разрыве в дисперсной фазе для обеих
          фаз для итерационного процесса */
    return_code = calc_p_v_initial_guess( paramsc, debug_info, left_params, right_params, sound_velocities,
        solver_part, &p_cont_gas, &v_cont_gas, &p_cont_disp, &v_cont_disp );
    if ( return_code == GODUNOV_FAILS )
        /* начальное приближение для итерационного процесса построить не удалось */
        return GODUNOV_FAILS;
    
    /* 2. Расчет давлений и скоростей обеих фаз на контактном разрыве в дисперсной фазе */
    /* инициализация давлений и скоростей по разные стороны от разрыва в дисперсной фазе */
    init_solid_discontinuity_pressures( p_cont_gas, p_cont_disp, solid_discontinuity_pressures );
    v1 = v_cont_disp;
    v2 = v_cont_gas;
    v21 = v_cont_gas;
    v22 = v_cont_gas;
    return_code = calc_p_v( paramsc, debug_info, left_params, right_params, sound_velocities,
        solver_part, solid_discontinuity_pressures, &v1, &v2, &v21, &v22 );
    if ( return_code == GODUNOV_FAILS )
        /* давления и скорости в задаче Римана найти не удалось */
        return GODUNOV_FAILS;
    cont_ncons[P_GAS] = solid_discontinuity_pressures[GAS_RIGHT];
    cont_ncons[P_DISP] = solid_discontinuity_pressures[DISP_RIGHT];
    cont_ncons[V_GAS] = v_cont_gas;
    cont_ncons[V_DISP] = v_cont_disp;
    /* 3. Отбор решения в дисперсной фазе */
    /* инициализация вектора-решения нулями */
    for ( i_comp = 0; i_comp < n; i_comp++ )
        v_ncons_res[i_comp] = 0.0;
    sample_solid_solution( paramsc, left_params, right_params, c1l, c1r, solid_discontinuity_pressures[DISP_LEFT],
        solid_discontinuity_pressures[DISP_RIGHT], v1, s, v_ncons_res, cont_ncons );

    /* 4. Отбор решения в газовой фазе */
    sample_gas_solution( paramsc, left_params, right_params, c2l, c2r, solid_discontinuity_pressures[GAS_LEFT],
        solid_discontinuity_pressures[GAS_RIGHT], v1, v2, v21, v22, s, v_ncons_res, cont_ncons );

    // 5. Расчет консервативной составляющей "потока"
	array1D tmp_flux( M );
    diff_flux_ncons( paramsc, v_ncons_res, &tmp_flux, 0 );
	for ( int i = 0; i < n; i++ )
		flux[i] = tmp_flux[i];

    /* 6. Сохранение параметров для аппроксимации неконсервативной составляющей вектора "потока" */
    *v_solid_cont = v1;
    *p_solid_cont_l = solid_discontinuity_pressures[DISP_LEFT];
    *p_solid_cont_r = solid_discontinuity_pressures[DISP_RIGHT];

    return SUCCESS;

}

// Построение начального приближения - давления и скорости на контактном разрыве в дисперсной фазе для обеих фаз
// для итерационного процесса
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_params[M] - вектор примитивных переменных слева от разрыва (in)
// right_params[M] - вектор примитивных переменных справа от разрыва (in)
// c - вектор со скоростями звука фаз по разные стороны от разрыва (in)
// p_cont_gas - начальное приближение для давлений в газе на контактном разрыве в дисперсной фазе (out)
// v_cont_gas - начальное приближение для скорости газа на контактном разрыве в дисперсной фазе (out)
// p_cont_disp - начальное приближение для давлений в дисперсной фазе на контактном разрыве в дисперсной фазе (out)
// v_cont_disp - начальное приближение для скорости дисперсной фазы на контактном разрыве в дисперсной фазе (out)
// solver_part - константа, которая определяет, какая из конфигураций объемной доли дисперсной фазы реализуется (in)
// Возвращает: SUCCEESS         начальное приближение построено успешно
//             GODUNOV_FAILS    начальное приближение построить не удалось
ReturnCodes calc_p_v_initial_guess( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M],
                                    double right_params[M], double c[K_GENERAL_CASE], Disp_phase_cases solver_part,
                                    double *p_cont_gas, double *v_cont_gas, double *p_cont_disp, double *v_cont_disp ) {

    double F, DF;                                   // функции, определяющей скорость среды на контактном разрыве, и ее производная
    double c1l = c[DISP_LEFT], c1r = c[DISP_RIGHT]; // скорости звука в дисперсной фазе слева и справа от первоначального разрыва
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   // скорости звука в газовой фазе слева и справа от первоначального разрыва
    int return_code;                                // код возврата

    // всегда находим давление и скорость на контактном разрыве в газовой фазе по стандартному алгоритму метода Годунова
    return_code = calc_contact_pressure_velocity( paramsc, debug_info, left_params, right_params, M, c2l, c2r, GAS_PHASE,
        p_cont_gas, v_cont_gas );
    if ( return_code == GODUNOV_FAILS )
        // не удалось построить начальное приближение для решения задачи Римана
        return GODUNOV_FAILS;
    
    // для дисперсной фазы - в зависимости от конфигурации объемной доли дисперной фазы
    switch ( solver_part ) {
        case BOTH_GRAD:
            // дисперсная фаза присутствует и справа, и слева от разрыва. Находим давление и скорость на контактном
            // разрыве в дисперсной фазе по стандартному алгоритму метода Годунова для случая двучленного уравнения состояния
            // Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
            // two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 499, formula (24).
            return_code = calc_contact_pressure_velocity( paramsc, debug_info, left_params, right_params, M, c1l, c1r,
                DISPERSED_PHASE, p_cont_disp, v_cont_disp );
            if ( return_code == GODUNOV_FAILS )
                // не удалось построить начальное приближение для решения задачи Римана
                return GODUNOV_FAILS;
            break;
        case LEFT_ONLY_GRAD:
            // дисперная фаза присутствует слева от разрыва и отсутствует справа от разрыва, причем перепад значительный, особый случай
            // Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
            // two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 501 + диссертация
            *p_cont_disp = *p_cont_gas;
            calc_F_and_DF( paramsc, *p_cont_disp, left_params, M, c1l, DISPERSED_PHASE, &F, &DF ); 
            *v_cont_disp = left_params[V_DISP] - F;
            break;
        case RIGHT_ONLY_GRAD:
            // дисперная фаза присутствует справа от разрыва и отсутствует слева от разрыва, причем перепад значительный, особый случай
            // Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
            // two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 501 + диссертация
            *p_cont_disp = *p_cont_gas;
            calc_F_and_DF( paramsc, *p_cont_disp, right_params, M, c1r, DISPERSED_PHASE, &F, &DF );
            *v_cont_disp = right_params[V_DISP] + F;
            break;
        default:
            printf( "\ncalc_p_v_initial_guess -> wrong solver_part value\n\n" );
            exit( EXIT_FAILURE );
    }

    return SUCCESS;

}

/* Расчет давлений и скоростей обеих фаз на контактном разрыве в дисперсной фазе
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   sound_velocities[K_GENERAL_CASE] - вектор со скоростями звука фаз по разные стороны от разрыва - 
                                      sound_velocities[GAS_LEFT] - в газе слева, sound_velocities[GAS_RIGHT] - в газе справа,
                                      sound_velocities[DISP_LEFT] - в дисперсной фазе слева, sound_velocities[DISP_RIGHT] - в дисперсной фазе справа (in)
   solver_part - константа, которая определяет, какая из конфигураций объемной доли дисперсной фазы реализуется (in)
   
   solid_discontinuity_pressures[K_GENERAL_CASE] - давление газа и дисперсной фазы по разные стороны от контактного разрыва
                                                   в дисперсной фазе, на входе - начальное приближение, на выходе - результат;
                                                   [GAS_LEFT] - давление слева в газе, [GAS_RIGHT] - давление справа в газе,
                                                   [DISP_LEFT] - давление слева в дисперсной фазе, [DISP_RIGHT] - давление справа в дисперсной фазе (in/out)
   v1 - скорость контактного разрыва в дисперсной фазе (out)
   v2 - скорость контактного разрыва в газовой фазе (out)
   v21 - скорость газа слева от контактного разрыва дисперсной фазы (out)
   v22 - скорость газа справа от контактного разрыва дисперсной фазы (out)
   
   Возвращает: SUCCEESS         успешная сходимость итераций
               GODUNOV_FAILS    итерации не сошлись, решение задачи Римана построить не удалось */ 
ReturnCodes calc_p_v( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                      double sound_velocities[K_GENERAL_CASE], Disp_phase_cases solver_part,
                      double solid_discontinuity_pressures[K_GENERAL_CASE], double *v1, double *v2, double *v21, double *v22 ) {
    
    int return_code;    /* код возврата функции */

    /* для дисперсной фазы - в зависимости от конфигурации объемной доли дисперной фазы */
    switch ( solver_part ) {
        case BOTH_GRAD:
            return_code = calc_solid_discontinuity_pressures( paramsc, debug_info, left_params, right_params, BOTH_GRAD,
                sound_velocities, solid_discontinuity_pressures, v1, v2, v21, v22 );
            debug_info->godunov_cases[1]++;
            break;
        case LEFT_ONLY_GRAD:
            /* дисперная фаза присутствует слева от разрыва и отсутствует справа от разрыва, особый случай
               Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of
               compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 502 + диссертация */
            return_code = calc_solid_discontinuity_pressures( paramsc, debug_info, left_params, right_params, LEFT_ONLY_GRAD,
                sound_velocities, solid_discontinuity_pressures, v1, v2, v21, v22 );
            /* давление в дисперсной фазе справа от разрыва не имеет смысла, оно кладется равным фоновому значению */
            solid_discontinuity_pressures[DISP_RIGHT] = paramsc->background_pressure;
            debug_info->godunov_cases[2]++;
            break;
        case RIGHT_ONLY_GRAD:
            /* дисперная фаза присутствует справа от разрыва и отсутствует слева от разрыва, особый случай
               Аналогично Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a
               model of compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 502 + диссертация */
            return_code = calc_solid_discontinuity_pressures( paramsc, debug_info, left_params, right_params, RIGHT_ONLY_GRAD,
                sound_velocities, solid_discontinuity_pressures, v1, v2, v21, v22 );
             /* давление в дисперсной фазе слева от разрыва не имеет смысла, оно кладется равным фоновому значению */
             solid_discontinuity_pressures[DISP_LEFT] = paramsc->background_pressure;
             debug_info->godunov_cases[2]++;
            break;
        default:
            printf( "\ncalc_p_v -> wrong solver_part value\n\n" );
            exit( EXIT_FAILURE );

    }

    // в случае стационарного контактного разрыва скорость контактного разрыва кладется равным малому числу,
    // чтобы избежать неопределенности в процедуре выбора решения
    if ( fabs( *v1 ) < paramsc->eps_contact ) *v1 = paramsc->eps_contact;
    if ( fabs( *v2 ) < paramsc->eps_contact ) *v2 = paramsc->eps_contact;
    if ( fabs( *v21 ) < paramsc->eps_contact ) *v21 = paramsc->eps_contact;
    if ( fabs( *v22 ) < paramsc->eps_contact ) *v22 = paramsc->eps_contact;

    return SUCCESS;

}

/* Рабочая функция заполнения вспомогательного массива со скоростями звука фаз слева и справа от разрыва

   c2l - скорость звука в газовой фазе слева от разрыва (in)
   c2r - скорость звука в газовой фазе справа от разрыва (in)
   c1l - скорость звука в дисперсной фазе слева от разрыва (in)
   c1r - скорость звука в дисперсной фазе справа от разрыва (in)

   sound_velocities[K_GENERAL_CASE] - вектор скоростей звука (out) */
void fill_sound_velocities( double c2l, double c2r, double c1l, double c1r, double sound_velocities[K_GENERAL_CASE] ) {

    sound_velocities[GAS_LEFT] = c2l;
    sound_velocities[GAS_RIGHT] = c2r;
    sound_velocities[DISP_LEFT] = c1l;
    sound_velocities[DISP_RIGHT] = c1r;

}

/* Рабочая функция инициализации вектора давлений слева и справа от разрыва в дисперсной фазе

   p_cont_gas - давление газа на контактном разрыве без учета разрыва объемной доли дисперсной фазы (in)
   p_cont_solid - давление дисперсной фазы на контактном разрыве без учета разрыва объемной доли дисперсной фазы (in)

   solid_discontinuity_pressures[K_GENERAL_CASE] - вектор начальных давлений слева и справа от разрыва объемной доли
   дисперсной фазы (out) */
void init_solid_discontinuity_pressures( double p_cont_gas, double p_cont_solid,
                                         double solid_discontinuity_pressures[K_GENERAL_CASE] ) {

    solid_discontinuity_pressures[GAS_LEFT] = p_cont_gas;
    solid_discontinuity_pressures[GAS_RIGHT] = p_cont_gas;
    solid_discontinuity_pressures[DISP_LEFT] = p_cont_solid;
    solid_discontinuity_pressures[DISP_RIGHT] = p_cont_solid;

}

/* Итерационная процедура расчета давления и скорости на контактном разрыве в среде без разрыва объемной доли дисперсной фазы

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 155. - Subroutine STARPU.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   v_ncons_l - вектор примитивных переменных слева от разрыва (in)
   v_ncons_r - вектор примитивных переменных справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   phase - идентификатор фазы, для которой ищется начальное приближение - газовая или дисперсная (in)

   p_cont - давление на контактном разрыве (out)
   v_cont - скорость на контактном разрыве (out)
   
   Возвращает: SUCCEESS         успешная сходимость итераций
               GODUNOV_FAILS    итерации не сошлись, решение задачи Римана построить не удалось */
int calc_contact_pressure_velocity( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double *v_ncons_l, double *v_ncons_r,
                                    int vector_size, double cl, double cr, Phase phase, double *p_cont, double *v_cont ) {

    double vl, vr;      /* скорости слева и справа от разрыва */
    double p_old;       /* значение давления на предыдущей итерации */
    double fl, fr;      /* значения функций */
    double fld, frd;    /* значения производных */
    int iter_num = 0;   /* количество проведенных итераций */
    double criteria;    /* переменная для определения сходимости */
    double g;           /* показатель адиабаты */
    double p_prev = 0;// давление на  предыдущий итерации через одну;
    if ( vector_size == M_REDUCTION ) {
        /* сокращенный однофазный вектор */
        vl = v_ncons_l[V];
        vr = v_ncons_r[V];
        printf("\n cl = %lf cr = %lf", cl, cr);
    }
    switch ( phase ) {
        case GAS_PHASE:
            if ( vector_size == M ) {
                /* полный двухфазный вектор */
                vl = v_ncons_l[V_GAS];
                vr = v_ncons_r[V_GAS];
            }
            g = paramsc->g2;
            break;
        case DISPERSED_PHASE:
            if ( vector_size == M ) {
                vl = v_ncons_l[V_DISP];
                vr = v_ncons_r[V_DISP];
            }
            g =paramsc->g1;
            break;
        default:
            printf( "\ncalc_contact_pressure_velocity -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    if ( 2.0 * ( cl + cr ) / ( g - 1.0 ) <= vr - vl ) {
        /* случай возникновения вакуума */
        printf( "\ncalc_contact_pressure_velocity -> vacuum is generated in Godunov flux calculation " );
        debug_print( debug_info, phase );
        exit( EXIT_FAILURE );
    }

    /* расчет начального приближения для давления */
    p_old = pressure_initial_guess( paramsc, v_ncons_l, v_ncons_r, vector_size, cl, cr, phase );
    if ( p_old < 0.0 ) {
        printf( "\ncalc_contact_pressure_velocity -> initial pressure guess is negative " );
        debug_print( debug_info, phase );
        printf(" \n pressure_old = %lf", p_old);
        exit( EXIT_FAILURE );
    }
    
    /* решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона */
    do {
        calc_F_and_DF( paramsc, p_old, v_ncons_l, vector_size, cl, phase, &fl, &fld );
        calc_F_and_DF( paramsc, p_old, v_ncons_r, vector_size, cr, phase, &fr, &frd );
        *p_cont = p_old - ( fl + fr + vr - vl ) / ( fld + frd );
        printf("\n p_cont = %lf, iter_num = %d fl = %lf, fr = %lf, vl = %lf, vr = %lf, fld = %lf, frd = %lf", *p_cont, iter_num, fl, fr, vl, vr, fld, frd);
        criteria = 2.0 * fabs( ( *p_cont - p_old ) / ( *p_cont + p_old ) );
        if(iter_num == 0)
            p_prev = *p_cont;
        printf( "\n %lf %lf %lf", p_prev, p_old, *p_cont);
        if( iter_num>1){
            if(*p_cont - p_prev < paramsc->eps_general){
                *p_cont = max(*p_cont, p_old);
                criteria = paramsc->eps_general;
                printf("check");
            }
            p_prev = p_old;
        }
        
        iter_num++;
        if ( iter_num > paramsc->max_iter_num ) {
            printf( "\ncalc_contact_pressure_velocity -> number of iterations exceeds the maximum value " );
            debug_print( debug_info, phase );
            return GODUNOV_FAILS;
        }
        if ( *p_cont < 0.0 && p_old < 0) {
            printf( "\ncalc_contact_pressure_velocity -> pressure is negative \n" );
            debug_print( debug_info, phase );
            printf(" \n pressure_cont = %lf", *p_cont);
            /* построение адиабат для проверки */
            /* draw_adiabatic_curve( params, v_ncons_l, vector_size, cl, phase, LEFT );
            draw_adiabatic_curve( params, v_ncons_r, vector_size, cr, phase, RIGHT );
             */exit( EXIT_FAILURE );
            return GODUNOV_FAILS;
        }
        p_old = *p_cont;
    } while ( criteria > paramsc->eps_general );

    /* скорость контактного разрыва */
    *v_cont = 0.5 * ( vl + vr + fr - fl );
    return SUCCESS;

}

/* Определение начального приближения для расчета давления на контактном разрыве в среде без разрыва объемной доли дисперсной фазы

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l - вектор примитивных переменных слева от разрыва (in)
   v_ncons_r - вектор примитивных переменных справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   phase - идентификатор фазы, для которой ищется начальное приближение - газовая или дисперсная (in)

   Возвращает искомое начальное приближения */
double pressure_initial_guess( struct ParametersCommon *paramsc, double *v_ncons_l, double *v_ncons_r, int vector_size,
                               double cl, double cr, int phase ) {

    double rl, vl, pl;                  /* примитивные переменные слева от разрыва */
    double rr, vr, pr;                  /* примитивные переменные справа от разрыва */
    double g;                           /* показатель адиабаты */
    double p01;                          /* поправка для двучленного уравнения состояния */
    /* начальное приближение, рассчитанное на освановании рассмотрения линеаризованной системы
       в примитивных переменных */
    double p_lin;
    double p_min, p_max;                /* минимальное и максимальное давления слева и справа от разрыва */
    double p_ratio;                     /* перепад по давлению слева и справа от разрыва */
    double p1, p2, g1, g2;              /* вспомогательные переменные для промежуточных расчетов */
    double p_cand;                      /* кандидат на роль начального приближения */
    
    if ( vector_size == M_REDUCTION ) {
        /* сокращенный однофазный вектор */
        rl = v_ncons_l[R];
        vl = v_ncons_l[V];
        pl = v_ncons_l[P];
        rr = v_ncons_r[R];
        vr = v_ncons_r[V];
        pr = v_ncons_r[P];
    }
    switch ( phase ) {
        case GAS_PHASE:
            if ( vector_size == M ) {
                rl = v_ncons_l[R_GAS];
                vl = v_ncons_l[V_GAS];
                pl = v_ncons_l[P_GAS];
                rr = v_ncons_r[R_GAS];
                vr = v_ncons_r[V_GAS];
                pr = v_ncons_r[P_GAS];
            }
            g = paramsc->g2;
            p01 = 0.0;
            break;
        case DISPERSED_PHASE:
            if ( vector_size == M ) {
                rl = v_ncons_l[R_DISP];
                vl = v_ncons_l[V_DISP];
                pl = v_ncons_l[P_DISP];
                rr = v_ncons_r[R_DISP];
                vr = v_ncons_r[V_DISP];
                pr = v_ncons_r[P_DISP];
            }
            g = paramsc->g1;
            p01 =paramsc->p01;
            break;
        default:
            printf( "\npressure_initial_guess -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    /* Начальное приближение из линейной задачи
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = min( pl, pr );
    p_max = max( pl, pr );
    p_ratio = p_max / p_min;
    printf("\n p_ratio = %lf", p_ratio);
    if ( ( p_ratio <= paramsc->p_max_ratio ) &&
        ( ( p_min < p_lin && p_lin < p_max ) || ( fabs( p_min - p_lin ) < paramsc->eps_general || fabs( p_max - p_lin ) < paramsc->eps_general ) ) ) {
        /* Начальное приближение из линеаризованной задачи */
        p_cand = p_lin;
        printf("\n pressure initial guess linear p_lin = %lf", p_lin);
    } else {
        if ( p_lin < p_min ) {
            /* Начальное приближение по двум волнам разрежения
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 301. - Formula (9.32) + поправка на двучленное уравнение состояния */
            g1 = 0.5 * ( g - 1.0 ) / g;
            p_cand = pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / pow( pl + p01, g1 ) + cr / pow( pr + p01, g1 ) ) ),
                1 / g1 ) - p01;
            printf("\n pressure initial guess two fan p_cand = %lf  pr = %lf  pl = %lf cr = %lf cl = %lf vr = %lf, vl = %lf", p_cand, pr, pl, cr, cl, vr, vl);
        } else {
            /* Начальное приближение по двум ударным волнам
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48) + поправка на двучленное уравнение состояния */
            g1 = 2.0 / ( g + 1.0 );
            g2 = ( g - 1.0 ) / ( g + 1.0 );
            p1 = sqrt( g1 / rl / ( g2 * ( pl + p01 ) + p_lin + p01 ) );
            p2 = sqrt( g1 / rr / ( g2 * ( pr + p01 ) + p_lin + p01 ) );
            p_cand = ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
            printf("\n pressure initial guess two shock p_cand = %lf  pr = %lf  pl = %lf", p_cand, pr, pl);
        }
    }
    
    if ( p_cand < 0.0 ) {
        /* Если не сработали методики Toro, то пробуем давление из "звукового распада разрыва"
           Годунов С.К. и др. Численное решение многомерных задач газовой динамики. - 
           М.: Наука, 1976. - С. 113. - Формула (13.26). */
        p_cand = ( pl * rr * cr + pr * rl * cl + ( vl - vr ) * rl * cl * rr * cr ) / ( rl * cl + rr * cr );
        printf("\n pressure initial guess sonic discontinue p_cand = %lf  pr = %lf  pl = %lf ", p_cand, pr, pl);
    }

    if ( p_cand < 0.0 ) {
        /* Если ничего не помогло, пытаемся построить хоть какое-то начальное приближение */
        p_cand = 0.5 * ( pl + pr );
        printf("\n pressure initial guess average p_cand = %lf  pr = %lf  pl = %lf ", p_cand, pr, pl);
    }
    printf("\n final guess p_cand = %lf", p_cand);
    return p_cand;

}

/* Расчет функции F, определяющей скорость среды на контактном разрыве в среде без разрыва объемной доли дисперсной фазы, и ее производной по давлению среды DF

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + поправка на двучленное уравнение состояния:
   Годунов С.К. и др. Численное решение многомерных задач газовой динамики. - 
   М.: Наука, 1976. - С. 110 - 111. - Формулы (13.16), (13.17).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   curr_press - давление с предыдущей итерации (in)
   v_ncons - вектор примитивных переменных (in)
   vector_size - размер вектора переменных (in)
   c - скорость звука в невозмущенной среде (in)
   phase - идентификатор фазы, для которой ищется начальное приближение - газовая или дисперсная (in)

   F - значение функции (out)
   DF - значение производной (out) */
void calc_F_and_DF( struct ParametersCommon *paramsc, double curr_press, double *v_ncons, int vector_size, double c, int phase, double *F, double *DF ) {

    double r, v, p; // примитивные переменные
    double g; // показатель адиабаты
    double p01; // поправка для двучленного уравнения состояния
    double p_ratio, fg, q; // вспомогательные переменные

    if ( vector_size == M_REDUCTION ) {
        // сокращенный однофазный вектор
        r = v_ncons[R];
        v = v_ncons[V];
        p = v_ncons[P];
     }
    switch ( phase ) {
        case GAS_PHASE:
            if ( vector_size == M ) {
                r = v_ncons[R_GAS];
                v = v_ncons[V_GAS];
                p = v_ncons[P_GAS];
            }
            g = paramsc->g2;
            p01 = 0.0;
            break;
        case DISPERSED_PHASE:
            if ( vector_size == M ) {
                r = v_ncons[R_DISP];
                v = v_ncons[V_DISP];
                p = v_ncons[P_DISP];
            }
            g = paramsc->g1;
            p01 =paramsc->p01;
            break;
        default:
            printf( "\ncalc_F_and_DF -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    p_ratio = ( curr_press + p01 ) / ( p + p01 );
    if ( curr_press <= p ) {
        // волна разрежения
        fg = 2.0 / ( g - 1.0 );
        *F = fg * c * ( pow( p_ratio, 1.0 / fg / g ) - 1.0 );
        *DF = ( 1.0 / r / c ) * pow( p_ratio, - 0.5 * ( g + 1.0 ) / g );
    }
    else {
        // ударная волна
        double a, b;
        a = 2.0 / r / ( g + 1.0 );
        b = ( g - 1.0 ) * p / ( g + 1.0 );
        q = sqrt( a / ( b + curr_press ) );
        *F = ( curr_press - p ) * q;
        *DF = ( 1.0 - 0.5 * ( curr_press - p ) / ( b + curr_press ) ) * q;
        /*q = sqrt( 0.5 * ( g + 1.0 ) / g * p_ratio + 0.5 * ( g - 1.0 ) / g );
        *F = ( curr_press - p ) / c / r / q;
        *DF = 0.25 * ( ( g + 1.0 ) * p_ratio + 3 * g - 1.0 ) / g / r / c / pow( q, 3.0 );*/
    }

}

/* Расчет функции G, определяющей плотность газа по разные стороны от контактного разрыва в среде без разрыва объемной доли дисперсной фазы,
   и ее производной по давлению DG

   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - Функция G определяется формулами (7).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   curr_press - давление с предыдущей итерации (in)
   v_ncons[M] - вектор примитивных переменных (in)
   
   G - значение функции (out)
   DG - значение производной (out) */
void calc_G_and_DG( struct ParametersCommon *paramsc, double curr_press, double v_ncons[M], double *G, double *DG ) {

    double r = v_ncons[R_GAS], p = v_ncons[P_GAS]; // плотность и давление газа в невозмущенной среде
    double p_ratio; // отношение текущего давления на контактном разрыве к давлению в невозмущенной среде
    double g1 = paramsc->g2 - 1.0, g2 = paramsc->g2 + 1.0; // переменные для экономии вычислений
    double A, B;
        
    p_ratio = curr_press / p;
    if ( curr_press < p ) {
        // волна разрежения
        *G = r * pow( p_ratio, 1.0 / paramsc->g2 );
        *DG = *G / paramsc->g2 / curr_press;
    }
    else {
        // ударная волна
        A = g1 + g2 * p_ratio;
        B = g1 * p_ratio + g2;
        *G = r * A / B;
        *DG = r * ( g2 * B - g1 * A ) / p / pow( B, 2.0 );
    }

}

/* Итерационная процедура для поиска давлений в газовой и дисперсной фазах на контактном разрыве в дисперсной фазе
   методом Ньютона-Рафсона
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   dsolver_part - константа, которая определяет, какая из конфигураций объемной доли дисперсной фазы реализуется (in)
   sound_velocities[K_GENERAL_CASE] - вектор со скоростями звука фаз по разные стороны от разрыва - 
                                      sound_velocities[GAS_LEFT] - в газе слева, sound_velocities[GAS_RIGHT] - в газе справа,
                                      sound_velocities[DISP_LEFT] - в дисперсной фазе слева,
                                      sound_velocities[DISP_RIGHT] - в дисперсной фазе справа (in)

   solid_discontinuity_pressures[K_GENERAL_CASE] - давление газа и дисперсной фазы по разные стороны от контактного разрыва
                                                   в дисперсной фазе, на входе - начальное приближение, на выходе - результат;
                                                   [GAS_LEFT] - давление слева в газе, [GAS_RIGHT] - давление справа в газе,
                                                   [DISP_LEFT] - давление слева в дисперсной фазе, [DISP_RIGHT] - давление справа в дисперсной фазе (in/out)
   v1 - скорость контактного разрыва в дисперсной фазе (out)
   v2 - скорость контактного разрыва в газовой фазе (out)
   v21 - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v22 - скорость газа справа от контактного разрыва в дисперсной фазе (out)
   
   Возвращает: SUCCEESS         успешная сходимость итераций
               GODUNOV_FAILS    итерации не сошлись, решение задачи Римана построить не удалось */
int calc_solid_discontinuity_pressures( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M],
                                        double right_params[M], Disp_phase_cases solver_part, double sound_velocities[K_GENERAL_CASE],
                                        double solid_discontinuity_pressures[K_GENERAL_CASE], double *v1, double *v2,
                                        double *v21, double *v22 ) {

    int actual_sol_cont_pres_size; /* количество рассчитываемых давлений - K_GENERAL_CASE в общем случае и K_SPECIAL_CASE в особом */
    // all arrays are of the maximum possible constant size
    double dp[M]; /* приращение вектора-решения в итерационном процессе */
    double P[M]; /* вектор-функция, определяющая соотношения на разрыве объемной доли в дисперсной фазе */
    double DP[M][M]; 

    double criteria; /* переменная для определения сходимости */
    double average_pressures[M]; /* вектор со средними давлениями между итерациями - для расчета критерия останова */
    double check_system_sol[M]; /* вектор для проверки решения СЛАУ */

    // solution of system of non-linear algebraic equations to find both phases pressures on the solid discontinuity
    // by Newton-Raphson method

    /* определение размера решаемой системы */
    switch ( solver_part ) {
        case BOTH_GRAD:
            actual_sol_cont_pres_size = K_GENERAL_CASE; /* общий случай наличия дисперсной фазы по обе стороны от разрыва */
            break;
        case LEFT_ONLY_GRAD:
        case RIGHT_ONLY_GRAD:
            actual_sol_cont_pres_size = K_SPECIAL_CASE; /* особый случай отсутствия дисперсной фазы по одну из сторон от разрыва */
            break;
        default:
            printf( "\ncalc_solid_discontinuity_pressures -> wrong value of solver_part variable.\n" );
            exit( EXIT_FAILURE );
    }

    do {
        
        /* расчет определяющей функции и матрицы Якобы */
        switch ( solver_part ) {
            case BOTH_GRAD:
                /* общий случай наличия дисперсной фазы по обе стороны от разрыва */
                calc_P_and_DP_general_case( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case LEFT_ONLY_GRAD:
                /* особый случай - дисперсная фаза только слева от разрыва */
                calc_P_and_DP_left_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case RIGHT_ONLY_GRAD:
                /* особый случай - дисперсная фаза только справа от разрыва */
                calc_P_and_DP_right_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            default:
                printf( "\ncalc_solid_discontinuity_pressures -> wrong value of solver_part variable.\n" );
                exit( EXIT_FAILURE );
        }
                
        // right-hand side vector calculation for Newton iterations
        for ( int i = 0; i < actual_sol_cont_pres_size; i++ )
            P[i] = - P[i];
        
        // solution of the system for residuals
        solve_linear_system( DP, P, actual_sol_cont_pres_size, paramsc->eps_general, dp );
        if ( paramsc->is_debug ) {
            // check the correctness of system solution
            mult_matrix_vector( actual_sol_cont_pres_size, DP, dp, check_system_sol );
            if ( !compare_vectors( actual_sol_cont_pres_size, check_system_sol, P, paramsc->eps_general ) ) {
                printf( "\ncalc_solid_discontinuity_pressures -> linear system of equations is solved incorrectly.\n" );
                exit( EXIT_FAILURE );
            }
        }
        
        // solution vector update
        if ( solver_part == BOTH_GRAD || solver_part == LEFT_ONLY_GRAD ) {
            for ( int i = 0; i < actual_sol_cont_pres_size; i++ )
                /* обновляются 4 компоненты в общем случае или 3 подряд в случае наличия дисперной фазы только слева от разрыва */
                solid_discontinuity_pressures[i] += dp[i];
        }
        else {
            solid_discontinuity_pressures[GAS_LEFT] += dp[0];   /* давление газа слева от контактного разрыва в дисперсной фазе */
            solid_discontinuity_pressures[GAS_RIGHT] += dp[1];  /* давление газа справа от контактного разрыва в дисперсной фазе */
            solid_discontinuity_pressures[DISP_RIGHT] += dp[2]; /* давление дисперсной фазы справа от контактного разрыва в дисперсной фазе */
        }

        // stop criteria calculation
        if ( solver_part == BOTH_GRAD || solver_part == LEFT_ONLY_GRAD ) {
            for ( int i = 0; i < actual_sol_cont_pres_size; i++ )
                average_pressures[i] = solid_discontinuity_pressures[i] - 0.5 * dp[i];
        }
        else {
            average_pressures[0] = solid_discontinuity_pressures[GAS_LEFT] - 0.5 * dp[0];   /* соответствует давлению газа слева от
                                                                                               контактного разрыва в дисперсной фазе */
            average_pressures[1] = solid_discontinuity_pressures[GAS_RIGHT] - 0.5 * dp[1];  /* соответствует давлению газа справа от
                                                                                               контактного разрыва в дисперсной фазе */
            average_pressures[2] = solid_discontinuity_pressures[DISP_RIGHT] - 0.5 * dp[2]; /* соответствует давлению дисперсной фазы справа от
                                                                                               контактного разрыва в дисперсной фазе */
        }
        criteria = norm( dp, actual_sol_cont_pres_size ) / norm( average_pressures, actual_sol_cont_pres_size );

        // the necessary check of the correctness of the resulats
        for ( int i = 0; i < actual_sol_cont_pres_size; i++ ) {
            if ( solid_discontinuity_pressures[i] < 0.0 ) {
                
                /* возможность детального отслеживания ситуации отсутствия сходимости и сбора полной трассы */
                /* printf( "\ncalc_solid_discontinuity_pressures -> pressure is negative " );
                debug_print( debug_info, i ); */
                
                return GODUNOV_FAILS;
            }
        }

    } while ( criteria > paramsc->eps_general );

    if ( paramsc->is_debug ) {
        /* проверка, что найденное решение соответствует решаемой нелинейной системе */
        switch ( solver_part ) {
            case BOTH_GRAD:
                /* общий случай наличия дисперсной фазы по обе стороны от разрыва */
                calc_P_and_DP_general_case( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case LEFT_ONLY_GRAD:
                /* особый случай - дисперсная фаза только слева от разрыва */
                calc_P_and_DP_left_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case RIGHT_ONLY_GRAD:
                /* особый случай - дисперсная фаза только справа от разрыва */
                calc_P_and_DP_right_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            default:
                printf( "\ncalc_solid_discontinuity_pressures -> wrong value of disp_phase variable.\n" );
                exit( EXIT_FAILURE );
        }
        for ( int i = 0; i < actual_sol_cont_pres_size; i++ ) {
            if ( fabs( P[i] ) > paramsc->eps_thin_layer ) {
                printf( "\ncalc_solid_discontinuity_pressures -> the 'thin-layer' algebraic system of non-linear equations is solved wrong " );
                debug_print( debug_info, i );
                printf( "calc_solid_discontinuity_pressures -> P[%d] = %e instead of acceptable tolerance %e.\n", i, P[i], paramsc->eps_thin_layer );
                exit( EXIT_FAILURE );
            }
        }
    }
    return SUCCESS;

}

/* Расчет вектор-функции для определения соотношений на контактном разрыве в дисперсной фазе, ее матрицы Якоби, а также скоростей
   контактных разрывов и скоростей газа слева и справа от контактного разрыва в дисперсной фазе.
   Общий случай наличия дисперсной фазы по обе стороны от разрыва.
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - Искомая функция определяется формулами
   (23), (25).

   Общее правило обозначения переменных с двумя индексами:
   - первый индекс - номер фазы: 1 - дисперсная фаза, 2 - газовая фаза
   - второй индекс - положение относительно контактного разрыва: 1 - слева, 2 - справа
   Индексы l и r соответствуют параметрам в невозмущенной среде слева и справа от разрыва, соответственно.
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   c[K_GENERAL_CASE] - вектор скоростей звука в дисперсной и газовой фазах по разные стороны от разрыва (in)
   curr_p[K_GENERAL_CASE] - текущий вектор давлений дисперсной и газовой фаз слева и справа от разрыва пористости (in)

   P[K_GENERAL_CASE] - искомая вектор-функция (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - искомая матрица Якоби (out)
   v_cont_disp - скорость контактного разрыва в дисперсной фазе (out)
   v_cont_gas - скорость контактного разрыва в газовой фазе (out)
   v_gas_left - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v_gas_right - скорость газа справа от контактного разрыва в дисперсной фазе (out) */
void calc_P_and_DP_general_case( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                 double c[M], double curr_p[M], double P[M],
                                 double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                 double *v_gas_left, double *v_gas_right ) {

    double p11 = curr_p[DISP_LEFT], p12 = curr_p[DISP_RIGHT]; /* текущие давления в дисперсной фазе слева и справа от контактного разрыва в дисперсной фазе */
    double p21 = curr_p[GAS_LEFT], p22 = curr_p[GAS_RIGHT];   /* текущие давления в газовой фазе слева и справа от контактного разрыва в дисперсной фазе */
    
    double c1l = c[DISP_LEFT], c1r = c[DISP_RIGHT]; /* скорости звука в дисперсной фазе слева и справа от первоначального разрыва */
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   /* скорости звука в газовой фазе слева и справа от первоначального разрыва */

    double v1l = left_params[V_DISP], v1r = right_params[V_DISP];  /* скорости дисперсной фазы слева и справа от первоначального разрыва */
    double v2l = left_params[V_GAS], v2r = right_params[V_GAS];  /* скорости газовой фазы слева и справа от первоначального разрыва */

    double b1l = left_params[B_DISP], b1r = right_params[B_DISP];  /* объемные доли дисперсной фазы слева и справа от первоначального разрыва */
    double b2l = 1.0 - b1l, b2r = 1.0 - b1r;               /* объемные доли газовой фазы слева и справа от первоначального разрыва */

    double f11, df11;   /* функция и производная для определения скорости дисперсной фазы слева от контактного разрыва в дисперсной фазе */
    double f12, df12;   /* функция и производная для определения скорости дисперсной фазы справа от контактного разрыва в дисперсной фазе */

    double f21, df21;   /* функция и производная для определения скорости газовой фазы слева от контактного разрыва в дисперсной фазе */
    double f22, df22;   /* функция и производная для определения скорости газовой фазы справа от контактного разрыва в дисперсной фазе */

    double g21, dg21;   /* функция и производная для определения плотности газовой фазы слева от контактного разрыва в газе */
    double g22, dg22;   /* функция и производная для определения плотности газовой фазы слева от контактного разрыва в газе */

    double v11, v12;    /* текущие скорости дисперсной фазы слева и справа от контактного разрыва в дисперсной фазе */
    double v21, v22;    /* текущие скорости газовой фазы слева и справа от контактного разрыва в дисперсной фазе */

    double r1, r2;  /* плотности газа слева и справа от контактного разрыва в газе */

    double g, A, B, C, D, E, dv2, dv1, dv2c;    /* вспомогательные переменные для сокращения объема вычислений */

    /* по разности текущих скоростей контактных разрывов в дисперсной и газовой фазах определяем, какая из двух возможных
       конфигураций реализуется */
    /* дисперсная фаза */
    calc_F_and_DF( paramsc, p11, left_params, M, c1l, DISPERSED_PHASE, &f11, &df11 );
    v11 = v1l - f11;
    calc_F_and_DF( paramsc, p12, right_params, M, c1r, DISPERSED_PHASE, &f12, &df12 );
    v12 = v1r + f12;
    /* газовая фаза */
    calc_F_and_DF( paramsc, p21, left_params, M, c2l, GAS_PHASE, &f21, &df21 );
    v21 = v2l - f21;
    calc_F_and_DF( paramsc, p22, right_params, M, c2r, GAS_PHASE, &f22, &df22 );
    v22 = v2r + f22;

    /* расчет плотностей газа слева и справа от контактного разрыва в газе по текущему давлению на контактном разрыве */
    calc_G_and_DG( paramsc, p21, left_params, &g21, &dg21 );
    r1 = g21;
    calc_G_and_DG( paramsc, p22, right_params, &g22, &dg22 );
    r2 = g22;

    /* расчет вспомогательных переменных */
    g = paramsc->g2 / ( paramsc->g2 - 1.0 );
    A = pow( p22 / p21, 1.0 / paramsc->g2 );
    B = b1r * p12 + b2r * p22 - b1l * p11 - b2l * p21;
    dv2 = v22 - v12;
    dv1 = v21 - v11;
    dv2c = v22 - v21;
    C = 0.5 * ( pow( dv2, 2.0 ) - pow( dv1, 2.0 ) );
    D = b2r * A * dv2;
    E = g / r1 / A;

    /* расчет компонентов вектор-функции P */
    P[0] = v12 - v11;
    if ( v21 > v11 ) {
        P[1] = D - b2l * dv1;
        P[2] = B + b2l * r1 * dv1 * dv2c;
        P[3] = E * p22 - g * p21 / r1 + C;
        *v_cont_gas = v22;  /* текущая скорость контактного разрыва в газовой фазе */
    }
    else {
        P[1] = - b2l / A * dv1 + b2r * dv2;
        P[2] = B + b2r * r2 * dv2 * dv2c;
        P[3] = - g * p21 * A / r2 + g * p22 / r2 + C;
        *v_cont_gas = v21;  /* текущая скорость контактного разрыва в газовой фазе */
    }

    /* расчет компонентов матрицы Якоби */
    /* первая строка */
    DP[0][0] = 0.0;
    DP[0][1] = 0.0;
    DP[0][2] = df11;
    DP[0][3] = df12;
    if ( v21 > v11 ) {
        /* вторая строка */
        DP[1][0] = - D / p21 / paramsc->g2 + b2l * df21;
        DP[1][1] = D / p22 / paramsc->g2 + b2r * A * df22;
        DP[1][2] = - b2l * df11;
        DP[1][3] = - b2r * A * df12;
        /* третья строка */
        DP[2][0] = b2l * ( - 1.0 + dv1 * dv2c * dg21 + r1 * df21 * ( dv1 - dv2c ) );
        DP[2][1] = b2r + b2l * r1 * dv1 * df22;
        DP[2][2] = - b1l + b2l * r1 * dv2c * df11;
        DP[2][3] = b1r;
        /* четвертая строка */
        DP[3][0] = g / r1 * ( p22 / A / paramsc->g2 / p21 - 1.0 ) + dv1 * df21 + dg21 * g / pow( r1, 2.0 ) * ( p21 - p22 / A );
        DP[3][1] = 1.0 / r1 / A + dv2 * df22;
        DP[3][2] = - dv1 * df11;
        DP[3][3] = - dv2 * df12;
    }
    else {
        /* вторая строка */
        DP[1][0] = b2l / A * ( df21 - dv1 / paramsc->g2 / p21 );
        DP[1][1] = b2r * df22 + b2l * dv1 / A / paramsc->g2 / p22;
        DP[1][2] = - b2l / A * df11; DP[1][3] = - b2r * df12;
        /* третья строка */
        DP[2][0] = - b2l + b2r * r2 * dv2 * df21;
        DP[2][1] = b2r * ( 1.0 + dg22 * dv2 * dv2c + r2 * df22 * ( dv2c + dv2 ) );
        DP[2][2] = - b1l;
        DP[2][3] = b1r - b2r * r2 * dv2c * df12;
        /* четвертая строка */
        DP[3][0] = - A / r2 + dv1 * df21;
        DP[3][1] = g / r2 + dv2 * df22 - A * p21 / ( paramsc->g2 - 1.0 ) / r2 / p22 + dg22 * g / pow( r2, 2.0 ) * ( p21 * A - p22 );
        DP[3][2] = - dv1 * df11;
        DP[3][3] = - dv2 * df12;
    }

    *v_cont_disp = v11; /* текущая скорость контактного разрыва в дисперсной фазе */
    *v_gas_left = v21;  /* текущая скорость газа слева от контактного разрыва в дисперсной фазе */
    *v_gas_right = v22; /* текущая скорость газа справа от контактного разрыва в дисперсной фазе */

}

/* Функция отбора решения в дисперсной фазе

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
   + поправка на двучленное уравнение состояния: Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann
   problem and a high-resolution Godunov method for a model of compressible two-phase flow // Journal of Computational Physics.
   - 2006. - V. 212. - P. 490 - 526. - Формулы (9) - (11).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор неконсервативных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор неконсервативных переменных справа от разрыва (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   p1 - давление слева от контактного разрыва в дисперсной фазе (in)
   p2 - давление справа от контактного разрыва в дисперсной фазе (in)
   v_cont - скорость на контактном разрыве (in)
   s - значение x/t, для которого отбирается решение (in)

   v_ncons_res[M] - вектор неконсервативных переменных с отобранными компонентами для дисперсной фазы (out) */
void sample_solid_solution( struct ParametersCommon *paramsc, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr,
                            double p1, double p2, double v_cont, double s, double v_ncons_res[M], double cont_ncons[M] ) {

    double betal, rl, vl, pl;           /* примитивные переменные слева от разрыва */
    double betar, rr, vr, pr;           /* примитивные переменные справа от разрыва */
    double g1, g2, g3, g4, g5, g6, g7;  /* вспомогательные переменные, производные от показателя адиабаты,
                                           в соответствии с Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics.
                                           - 2nd Edition. - Springer, 1999. - P. 153. */

    /* скорости левых волн */
    double shl, stl;    /* скорости "головы" и "хвоста" левой волны разрежения */
    double sl;          /* скорость левой ударной волны */

    /* скорости правых волн */
    double shr, str;    /* скорости "головы" и "хвоста" правой волны разрежения */
    double sr;          /* скорость правой ударной волны */

    double cml, cmr;        /* скорости звука слева и справа от контактного разрыва */
    double c;               /* локальная скорость звука внутри волны разрежения */
    double p_ratio;
    double beta, r, v, p;   /* отобранные значения объемной доли, плотности, скорости и давления */

    /* вспомогательные переменные */
    /* параметры слева от разрыва */
    betal = v_ncons_l[B_DISP];
    rl = v_ncons_l[R_DISP];
    vl = v_ncons_l[V_DISP];
    pl = v_ncons_l[P_DISP];
    /* параметры справа от разрыва */
    betar = v_ncons_r[B_DISP];
    rr = v_ncons_r[R_DISP];
    vr = v_ncons_r[V_DISP];
    pr = v_ncons_r[P_DISP];
    /* производные от показателя адиабаты */
    g1 = 0.5 * ( paramsc->g1 - 1.0 ) / paramsc->g1;
    g2 = 0.5 * ( paramsc->g1 + 1.0 ) / paramsc->g1;
    g3 = 2.0 * paramsc->g1 / ( paramsc->g1 - 1.0 );
    g4 = 2.0 / ( paramsc->g1 - 1.0 );
    g5 = 2.0 / ( paramsc->g1 + 1.0 );
    g6 = ( paramsc->g1 - 1.0 ) / ( paramsc->g1 + 1.0 );
    g7 = 0.5 * ( paramsc->g1 - 1.0 );
    cont_ncons[B_DISP] = betal;
    if ( s <= v_cont ) {
        /* рассматриваемая точка - слева от контактного разрыва */
        beta = betal;
        if ( p1 <= pl ) {
            /* левая волна разрежения */
            cont_ncons[R_DISP] = rl * pow( p1 / pl, 1.0 / paramsc->g1 );
            shl = vl - cl;
            if ( s <= shl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( ( p1 + paramsc->p01 ) / ( pl + paramsc->p01 ), g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* параметры слева от контактного разрыва */
                    r = rl * pow( ( p1 + paramsc->p01 ) / ( pl + paramsc->p01 ), 1.0 / paramsc->g1 );
                    v = v_cont;
                    p = p1;
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = ( pl + paramsc->p01 ) * pow( c / cl, g3 ) - paramsc->p01;
                }
            }
        }
        else {
            /* левая ударная волна */
            p_ratio = ( p1 + paramsc->p01 ) / ( pl + paramsc->p01 );
            cont_ncons[R_DISP] = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p1;
            }
        }
    }
    else {
        /* рассматриваемая точка - справа от контактного разрыва */
        beta = betar;
        if ( p2 > pr ) {
            /* правая ударная волна */
            p_ratio = ( p2 + paramsc->p01 ) / ( pr + paramsc->p01 );
            cont_ncons[R_DISP] = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p2;
            }
        }
        else {
            /* правая волна разрежения */
            cont_ncons[R_DISP] = rr * pow( p1 / pr, 1.0 / paramsc->g1 );
            shr = vr + cr;
            if ( s >= shr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( ( p2 + paramsc->p01 ) / ( pr + paramsc->p01 ), g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* параметры справа от контактного разрыва */
                   r = rr * pow( ( p2 + paramsc->p01 ) / ( pr + paramsc->p01 ), 1.0 / paramsc->g1 );
                   v = v_cont;
                   p = p2;
               }
               else {
                    /* параметры внутри правой волны разрежения */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = ( pr + paramsc->p01 ) * pow( c / cr, g3 ) - paramsc->p01;
               }
            }
        }
    }
    
    /* формирование выходного вектора с результатом */
    v_ncons_res[B_DISP] = beta;
    v_ncons_res[R_DISP] = r;
    v_ncons_res[V_DISP] = v;
    v_ncons_res[P_DISP] = p;
    
}

/* Функция отбора решения в газовой фазе

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор неконсервативных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор неконсервативных переменных справа от разрыва (in)
   cl - скорость звука в газе слева от разрыва (in)
   cr - скорость звука в газе справа от разрыва (in)
   p1 - давление газа слева от контактного разрыва в дисперсной фазе (in)
   p2 - давление газа справа от контактного разрыва в дисперсной фазе (in)
   v_cont_solid - скорость контактного разрыва в дисперсной фазе (in)
   v_cont_gas - скорость контактного разрыва в газовой фазе (in)
   v1 - скорость газа слева от контактного разрыва в дисперсной фазе (in)
   v2 - скорость газа справа от контактного разрыва в дисперсной фазе (in)
   s - значение x/t, для которого отбирается решение (in)

   v_ncons_res[M] - отобранный вектор неконсервативных переменных (out) */
void sample_gas_solution( struct ParametersCommon *paramsc, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr, double p1,
                          double p2, double v_cont_solid, double v_cont_gas, double v1, double v2, double s, double v_ncons_res[M], double cont_ncons[M] ) {

    double rl, vl, pl;                  /* примитивные переменные слева от разрыва */
    double rr, vr, pr;                  /* примитивные переменные справа от разрыва */
    double g1, g2, g3, g4, g5, g6, g7;  /* вспомогательные переменные, производные от показателя адиабаты,
                                           в соответствии с Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics.
                                           - 2nd Edition. - Springer, 1999. - P. 153. */

    /* скорости левых волн */
    double shl, stl;    /* скорости "головы" и "хвоста" левой волны разрежения */
    double sl;          /* скорость левой ударной волны */

    /* скорости правых волн */
    double shr, str;    /* скорости "головы" и "хвоста" правой волны разрежения */
    double sr;          /* скорость правой ударной волны */

    double cml, cmr;    /* скорости звука слева и справа от контактного разрыва */
    double c;           /* локальная скорость звука внутри волны разрежения */
    double p_ratio;
    double r, v, p;     /* отобранные значения плотности, скорости и давления */

    /* вспомогательные переменные */
    /* параметры слева от разрыва */
    rl = v_ncons_l[R_GAS];
    vl = v_ncons_l[V_GAS];
    pl = v_ncons_l[P_GAS];
    /* параметры справа от разрыва */
    rr = v_ncons_r[R_GAS];
    vr = v_ncons_r[V_GAS];
    pr = v_ncons_r[P_GAS];
    /* производные от показателя адиабаты */    
    g1 = 0.5 * ( paramsc->g2 - 1.0 ) / paramsc->g2;
    g2 = 0.5 * ( paramsc->g2 + 1.0 ) / paramsc->g2;
    g3 = 2.0 * paramsc->g2 / ( paramsc->g2 - 1.0 );
    g4 = 2.0 / ( paramsc->g2 - 1.0 );
    g5 = 2.0 / ( paramsc->g2 + 1.0 );
    g6 = ( paramsc->g2 - 1.0 ) / ( paramsc->g2 + 1.0 );
    g7 = 0.5 * ( paramsc->g2 - 1.0 );

    if ( s <= v_cont_gas ) {
        /* рассматриваемая точка - слева от контактного разрыва в газовой фазе */
        if ( p1 <= pl ) {
            /* левая волна разрежения */
            cont_ncons[R_GAS] = rl * pow( p1 / pl, 1.0 / paramsc->g2 );
            shl = vl - cl;
            if ( s <= shl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( p1 / pl, g1 );
                stl = v1 - cml;
                if ( s > stl ) {
                    /* параметры слева от контактного разрыва */
                    if ( s <= v_cont_solid ) {
                        r = rl * pow( p1 / pl, 1.0 / paramsc->g2 );
                        v = v1;
                        p = p1;                        
                    }
                    else {
                        r = rl * pow( p2 / pl, 1.0 / paramsc->g2 );
                        v = v2;
                        p = p2;
                    }
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = pl * pow( c / cl, g3 );
                }
            }
        }
        else {
            /* левая ударная волна */
            if ( s <= v_cont_solid )
                p_ratio = p1 / pl;
            else
                p_ratio = p2 / pl;
            cont_ncons[R_GAS] = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                if ( s <= v_cont_solid ) {
                    v = v1;
                    p = p1;
                }
                else {
                    v = v2;
                    p = p2;
                }
            }
        }
    }
    else {
        /* рассматриваемая точка - справа от контактного разрыва */
        if ( p2 > pr ) {
            /* правая ударная волна */
            if ( s <= v_cont_solid )
                p_ratio = p1 / pr;
            else
                p_ratio = p2 / pr;
            cont_ncons[R_GAS] = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                if ( s <= v_cont_solid ) {
                    v = v1;
                    p = p1;
                }
                else {
                    v = v2;
                    p = p2;
                }
            }
        }
        else {
            /* правая волна разрежения */
            cont_ncons[R_GAS] = rr * pow( p1 / pr, 1.0 / paramsc->g2 );
            shr = vr + cr;
            if ( s >= shr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( p2 / pr, g1 );
               str = v2 + cmr;
               if ( s <= str ) {
                   /* параметры справа от контактного разрыва */
                   if ( s <= v_cont_solid ) {
                       r = rr * pow( p1 / pr, 1.0 / paramsc->g2 );
                       v = v1;
                       p = p1;
                   }
                   else {
                       r = rr * pow( p2 / pr, 1.0 / paramsc->g2 );
                       v = v2;
                       p = p2;
                   }
               }
               else {
                    /* параметры внутри правой волны разрежения */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = pr * pow( c / cr, g3 );
               }
            }
        }
    }
    
    /* формирование выходного вектора с результатом */
    v_ncons_res[R_GAS] = r;
    v_ncons_res[V_GAS] = v;
    v_ncons_res[P_GAS] = p;
    
}

/* Расчет неконсервативной составляющей вектора "потока"
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model
   of compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - Формула (30).
   
   v_solid_cont - скорость контактного разрыва в дисперсной фазе (in)
   beta_l - объемная доля дисперсной фазы слева от разрыва (in)
   beta_r - объемная доля дисперсной фазы справа от разрыва (in)
   p_solid_l - давление в дисперсной фазе слева от контактного разрыва в дисперсной фазе (in)
   p_solid_r - давление в дисперсной фазе справа от контактного разрыва в дисперсной фазе (in)
   
   v_ncons_term[M] - неконсервативная составляющая вектора "потока" (out) */
void calc_ncons_term( double v_solid_cont, double beta_l, double beta_r, double p_solid_l, double p_solid_r,
                      double v_ncons_term[M] ) {

    v_ncons_term[B_DISP] = - v_solid_cont * ( beta_r - beta_l );    /* поправка в уравнение переноса объемной доли дисперсной фазы */
    v_ncons_term[R_DISP] = 0.0;                                     /* поправка в уравнение неразрывности для дисперсной фазы */
    v_ncons_term[V_DISP] = p_solid_r * beta_r - p_solid_l * beta_l; /* поправка в ЗСИ для дисперсной фазы */
    v_ncons_term[P_DISP] = v_solid_cont * v_ncons_term[V_DISP];         /* поправка в ЗСЭ для дисперсной фазы */
    v_ncons_term[R_GAS] = 0.0;                                     /* поправка в уравнение неразрывности для газовой фазы */
    v_ncons_term[V_GAS] = - v_ncons_term[V_DISP];                      /* поправка в ЗСИ для газовой фазы */
    v_ncons_term[P_GAS] = - v_ncons_term[P_DISP];                      /* поправка в ЗСЭ для газовой фазы */

}

/* Расчет вектор-функции для определения соотношений на контактном разрыве в дисперсной фазе, ее матрицы Якоби, а также скоростей
   контактных разрывов и скоростей газа слева и справа от контактного разрыва в дисперсной фазе.
   ОСОБЫЙ СЛУЧАЙ ОТСУТСТВИЯ ДИСПЕРСНОЙ ФАЗЫ СПРАВА ОТ РАЗРЫВА.
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 501 - 502.

   Общее правило обозначения переменных с двумя индексами:
   - первый индекс - номер фазы: 1 - дисперсная фаза, 2 - газовая фаза
   - второй индекс - положение относительно контактного разрыва: 1 - слева, 2 - справа
   Индексы l и r соответствуют параметрам в невозмущенной среде слева и справа от разрыва, соответственно.
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   c[K_GENERAL_CASE] - вектор скоростей звука в дисперсной и газовой фазах по разные стороны от разрыва (in)
   curr_p[K_GENERAL_CASE] - текущий вектор давлений дисперсной и газовой фаз слева и справа от разрыва пористости (in)

   curr_p[GAS_LEFT] - давление газа слева от контактного разрыва в дисперсной фазе
   curr_p[GAS_RIGHT] - давление газа справа от контактного разрыва в дисперсной фазе
   curr_p[DISP_LEFT] - давление дисперсной фазы слева от контактного разрыва в дисперсной фазе

   P[K_GENERAL_CASE] - искомая вектор-функция (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - искомая матрица Якоби (out)
   v_cont_disp - скорость контактного разрыва в дисперсной фазе (out)
   v_cont_gas - скорость контактного разрыва в газовой фазе (out)
   v_gas_left - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v_gas_right - скорость газа справа от контактного разрыва в дисперсной фазе (out) */
void calc_P_and_DP_left_disp_phase( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                    double c[K_GENERAL_CASE], double curr_p[K_GENERAL_CASE], double P[M],
                                    double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                    double *v_gas_left, double *v_gas_right ) {

    /* давления  */
                                        
    double p11 = curr_p[DISP_LEFT];                         /* текущее давление в дисперсной фазе слева от контактного разрыва
                                                               в дисперсной фазе */
    double p21 = curr_p[GAS_LEFT], p22 = curr_p[GAS_RIGHT]; /* текущие давления в газовой фазе слева и справа от контактного разрыва
                                                               в дисперсной фазе */
    
    /* скорости звука */

    double c1l = c[DISP_LEFT];                      /* скорость звука в дисперсной фазе слева от первоначального разрыва */
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   /* скорости звука в газовой фазе слева и справа от первоначального разрыва */

    /* скорости фаз */

    double v1l = left_params[V_DISP];                           /* скорость дисперсной фазы слева от первоначального разрыва */
    double v2l = left_params[V_GAS], v2r = right_params[V_GAS];   /* скорости газовой фазы слева и справа от первоначального разрыва */

    double b1l = left_params[B_DISP];   /* объемная доля дисперсной фазы слева от первоначального разрыва */
    double b2l = 1.0 - b1l;         /* объемная доля газовой фазы слева от первоначального разрыва */

    double f11, df11;  /* функция и производная для определения скорости дисперсной фазы слева
                           от контактного разрыва в дисперсной фазе */
    
    double f21, df21;  /* функция и производная для определения скорости газовой фазы слева от контактного разрыва в дисперсной фазе */
    double f22, df22;  /* функция и производная для определения скорости газовой фазы справа от контактного разрыва в дисперсной фазе */

    double g21, dg21;  /* функция и производная для определения плотности газовой фазы слева от контактного разрыва в газе */
    double g22, dg22;  /* функция и производная для определения плотности газовой фазы слева от контактного разрыва в газе */

    double v11;       /* текущая скорость дисперсной фазы слева от контактного разрыва в дисперсной фазе */
    double v21, v22;  /* текущие скорости газовой фазы слева и справа от контактного разрыва в дисперсной фазе */

    double r1, r2;  /* плотности газа слева и справа от контактного разрыва в газе */

    double g, A, B, C, D, E, dvL, dv1, dv2c;    /* вспомогательные переменные для сокращения объема вычислений */

    /* по разности текущих скоростей контактных разрывов в дисперсной и газовой фазах будет определено, какая из двух возможных
       конфигураций будет реализовываться. Поэтому - ищем скорости контактных разрывов. */

    /* дисперсная фаза */
    calc_F_and_DF( paramsc, p11, left_params, M, c1l, DISPERSED_PHASE, &f11, &df11 );
    v11 = v1l - f11;
    /* газовая фаза */
    calc_F_and_DF( paramsc, p21, left_params, M, c2l, GAS_PHASE, &f21, &df21 );
    v21 = v2l - f21;
    calc_F_and_DF( paramsc, p22, right_params, M, c2r, GAS_PHASE, &f22, &df22 );
    v22 = v2r + f22;

    /* расчет плотностей газа слева и справа от контактного разрыва в газе по текущему давлению на контактном разрыве */
    calc_G_and_DG( paramsc, p21, left_params, &g21, &dg21 );
    r1 = g21;
    calc_G_and_DG( paramsc, p22, right_params, &g22, &dg22 );
    r2 = g22;

    /* расчет вспомогательных переменных */
    g = paramsc->g2 / ( paramsc->g2 - 1.0 );
    A = pow( p22 / p21, 1.0 / paramsc->g2 );
    B = p22 - b1l * p11 - b2l * p21;
    dvL = v22 - v11;
    dv1 = v21 - v11;
    dv2c = v22 - v21;
    C = 0.5 * ( pow( dvL, 2.0 ) - pow( dv1, 2.0 ) );
    D = A * dvL;
    E = g / r1 / A;

    /* расчет компонентов вектор-функции P */
    if ( v21 > v11 ) {
        P[0] = D - b2l * dv1;
        P[1] = B + b2l * r1 * dv1 * dv2c;
        P[2] = E * p22 - g * p21 / r1 + C;
        *v_cont_gas = v22;  /* текущая скорость контактного разрыва в газовой фазе */
    }
    else {
        /* ветка отлажена, тест 16 */
        P[0] = dvL - b2l * dv1 / A;
        P[1] = B + r2 * dvL * dv2c;
        P[2] = - g * p21 * A / r2 + g * p22 / r2 + C;
        *v_cont_gas = v21;  /* текущая скорость контактного разрыва в газовой фазе */
    }

    /* расчет компонентов матрицы Якоби */
    if ( v21 > v11 ) {
        /* первая строка */
        DP[0][0] = - D / p21 / paramsc->g2 + b2l * df21;
        DP[0][1] = D / p22 / paramsc->g2 + A * df22;
        DP[0][2] = ( A - b2l ) * df11;
        /* вторая строка */
        DP[1][0] = b2l * ( - 1.0 + dv1 * dv2c * dg21 + r1 * df21 * ( dv1 - dv2c ) );
        DP[1][1] = 1.0 + b2l * r1 * dv1 * df22;
        DP[1][2] = - b1l + b2l * r1 * dv2c * df11;
        /* третья строка */
        DP[2][0] = g / r1 * ( p22 / A / paramsc->g2 / p21 - 1.0 ) + dv1 * df21 + dg21 * g / pow( r1, 2.0 ) * ( p21 - p22 / A );
        DP[2][1] = 1.0 / r1 / A + dvL * df22;
        DP[2][2] = dv2c * df11;
    }
    else {
        /* ветка отлажена, тест 16 */
        /* первая строка */
        DP[0][0] = b2l / A * ( df21 - dv1 / paramsc->g2 / p21 );
        DP[0][1] = df22 + b2l * dv1 / A / paramsc->g2 / p22;
        DP[0][2] = ( 1.0 - b2l / A ) * df11;
        /* вторая строка */
        DP[1][0] = - b2l + r2 * dvL * df21;
        DP[1][1] = 1.0 + dg22 * dvL * dv2c + r2 * df22 * ( dv2c + dvL );
        DP[1][2] = - b1l + r2 * dv2c * df11;
        /* третья строка */
        DP[2][0] = - 1.0 / r2 * A + dv1 * df21;
        DP[2][1] = g / r2 + dvL * df22 - A * p21 / ( paramsc->g2 - 1.0 ) / r2 / p22 + dg22 * g / pow( r2, 2.0 ) * ( p21 * A - p22 );
        DP[2][2] = dv2c * df11;
    }

    *v_cont_disp = v11;     /* текущая скорость контактного разрыва в дисперсной фазе */
    *v_gas_left = v21;      /* текущая скорость газа слева от контактного разрыва в дисперсной фазе */
    *v_gas_right = v22;     /* текущая скорость газа справа от контактного разрыва в дисперсной фазе */

}

/* Расчет вектор-функции для определения соотношений на контактном разрыве в дисперсной фазе, ее матрицы Якоби, а также скоростей
   контактных разрывов и скоростей газа слева и справа от контактного разрыва в дисперсной фазе.
   ОСОБЫЙ СЛУЧАЙ ОТСУТСТВИЯ ДИСПЕРСНОЙ ФАЗЫ СЛЕВА ОТ РАЗРЫВА.
   
   Формулы приведены в диссертации.

   Общее правило обозначения переменных с двумя индексами:
   - первый индекс - номер фазы: 1 - дисперсная фаза, 2 - газовая фаза
   - второй индекс - положение относительно контактного разрыва: 1 - слева, 2 - справа
   Индексы l и r соответствуют параметрам в невозмущенной среде слева и справа от разрыва, соответственно.
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   c[K_GENERAL_CASE] - вектор скоростей звука в дисперсной и газовой фазах по разные стороны от разрыва (in)
   curr_p[K_GENERAL_CASE] - текущий вектор давлений дисперсной и газовой фаз слева и справа от разрыва пористости (in)
   
   curr_p[GAS_LEFT] - давление газа слева от контактного разрыва в дисперсной фазе
   curr_p[GAS_RIGHT] - давление газа справа от контактного разрыва в дисперсной фазе
   curr_p[DISP_RIGHT] - давление дисперсной фазы справа от контактного разрыва в дисперсной фазе

   P[K_GENERAL_CASE] - искомая вектор-функция (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - искомая матрица Якоби - последние строка и столбец не имеют смысла (out)
   v_cont_disp - скорость контактного разрыва в дисперсной фазе (out)
   v_cont_gas - скорость контактного разрыва в газовой фазе (out)
   v_gas_left - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v_gas_right - скорость газа справа от контактного разрыва в дисперсной фазе (out) */
void calc_P_and_DP_right_disp_phase( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                     double c[M], double curr_p[M], double P[M],
                                     double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                     double *v_gas_left, double *v_gas_right ) {

    double p12 = curr_p[DISP_RIGHT];                        /* текущее давление в дисперсной фазе справа от контактного разрыва
                                                               в дисперсной фазе */
    double p21 = curr_p[GAS_LEFT], p22 = curr_p[GAS_RIGHT]; /* текущие давления в газовой фазе слева и справа от контактного разрыва
                                                               в дисперсной фазе */
    
    double c1r = c[DISP_RIGHT];                     /* скорость звука в дисперсной фазе справа от первоначального разрыва */
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   /* скорости звука в газовой фазе слева и справа от первоначального разрыва */

    double v1r = right_params[V_DISP];                          /* скорость дисперсной фазы справа от первоначального разрыва */
    double v2l = left_params[V_GAS], v2r = right_params[V_GAS];   /* скорости газовой фазы слева и справа от первоначального разрыва */

    double b1r = right_params[B_DISP];  /* объемная доля дисперсной фазы справа от первоначального разрыва */
    double b2r = 1.0 - b1r;         /* объемная доля газовой фазы слева от первоначального разрыва */

    double f12, df12;   /* функция и производная для определения скорости дисперсной фазы справа от контактного разрыва в дисперсной фазе */
    
    double f21, df21;   /* функция и производная для определения скорости газовой фазы слева от контактного разрыва в дисперсной фазе */
    double f22, df22;   /* функция и производная для определения скорости газовой фазы справа от контактного разрыва в дисперсной фазе */

    double g21, dg21;   /* функция и производная для определения плотности газовой фазы слева от контактного разрыва в газе */
    double g22, dg22;   /* функция и производная для определения плотности газовой фазы слева от контактного разрыва в газе */

    double v12;         /* текущая скорость дисперсной фазы справа от контактного разрыва в дисперсной фазе */
    double v21, v22;    /* текущие скорости газовой фазы слева и справа от контактного разрыва в дисперсной фазе */

    double r1, r2;  /* плотности газа слева и справа от контактного разрыва в газе */

    double g, A, B, C, D, E, dv2, dvR, dv2c;  /* вспомогательные переменные для сокращения объема вычислений */

    // по разности текущих скоростей контактных разрывов в дисперсной и газовой фазах будет определено, какая из двух возможных
    // конфигураций будет реализовываться. Поэтому - ищем скорости контактных разрывов.

    // дисперсная фаза
    calc_F_and_DF( paramsc, p12, right_params, M, c1r, DISPERSED_PHASE, &f12, &df12 );
    v12 = v1r + f12;
    // газовая фаза
    calc_F_and_DF( paramsc, p21, left_params, M, c2l, GAS_PHASE, &f21, &df21 );
    v21 = v2l - f21;
    calc_F_and_DF( paramsc, p22, right_params, M, c2r, GAS_PHASE, &f22, &df22 );
    v22 = v2r + f22;

    // расчет плотностей газа слева и справа от контактного разрыва в газе по текущему давлению на контактном разрыве
    calc_G_and_DG( paramsc, p21, left_params, &g21, &dg21 );
    r1 = g21;
    calc_G_and_DG( paramsc, p22, right_params, &g22, &dg22 );
    r2 = g22;

    // расчет вспомогательных переменных
    g = paramsc->g2 / ( paramsc->g2 - 1.0 );
    A = pow( p22 / p21, 1.0 / paramsc->g2 );
    B = b1r * p12 + b2r * p22 - p21;
    dv2 = v22 - v12;
    dvR = v21 - v12;
    dv2c = v22 - v21;
    C = 0.5 * ( pow( dv2, 2.0 ) - pow( dvR, 2.0 ) );
    D = A * dv2;
    E = g / r1 / A;

    // расчет компонентов вектор-функции P
    if ( v21 > v12 ) {
        // ветка отлажена, тест 17
        P[0] = b2r * D - dvR;
        P[1] = B + r1 * dvR * dv2c;
        P[2] = E * p22 - g * p21 / r1 + C;
        *v_cont_gas = v22; // текущая скорость контактного разрыва в газовой фазе
    }
    else {
        P[0] = b2r * dv2 - dvR / A;
        P[1] = B + b2r * r2 * dv2 * dv2c;
        P[2] = - g * p21 * A / r2 + g * p22 / r2 + C;
        *v_cont_gas = v21; // текущая скорость контактного разрыва в газовой фазе
    }

    // расчет компонентов матрицы Якоби
    if ( v21 > v12 ) {
        // ветка отлажена, тест 17
        // первая строка
        DP[0][0] = - b2r * D / p21 / paramsc->g2 + df21;
        DP[0][1] = b2r * ( D / p22 / paramsc->g2 + A * df22 );
        DP[0][2] = ( 1.0 - b2r * A ) * df12;
        // вторая строка
        DP[1][0] = - 1.0 - r1 * dv2c * df21 + r1 * dvR * df21 + dv2c * dvR * df21;
        DP[1][1] = b2r + r1 * dvR * df22;
        DP[1][2] = b1r - r1 * dv2c * df12;
        // третья строка
        DP[2][0] = g / r1 * ( p22 / A / paramsc->g2 / p21 - 1.0 ) + dg21 * g / pow( r1, 2.0 ) * ( p21 - p22 / A ) + dvR * df21;
        DP[2][1] = 1.0 / r1 / A + dv2 * df22;
        DP[2][2] = - dv2c * df12;
    }
    else {
        // первая строка
        DP[0][0] = ( df21 - dvR / paramsc->g2 / p21 ) / A;
        DP[0][1] = b2r * df22 + dvR / A / paramsc->g2 / p22;
        DP[0][2] = ( 1.0 / A - b2r ) * df12;
        // вторая строка
        DP[1][0] = - 1.0 + b2r * r2 * dv2 * df21;
        DP[1][1] = b2r * ( 1.0 + dv2 * dv2c * dg22 + r2 * df22 * ( dv2c + dv2 ) );
        DP[1][2] = b1r - b2r * r2 * dv2c * df12;
        // третья строка
        DP[2][0] = - ( paramsc->g2 + 1.0 ) / ( paramsc->g2 - 1.0 ) / r2 / A + dvR * df21;
        DP[2][1] = g / r2 * ( 1.0 + p21 / paramsc->g2 / p22 / A ) - g * dg22 / pow( r2, 2.0 ) * ( p22 - p21 / A ) + dv2 * df22;
        DP[2][2] = - dv2c * df12;
    }

    *v_cont_disp = v12; // текущая скорость контактного разрыва в дисперсной фазе
    *v_gas_left = v21; // текущая скорость газа слева от контактного разрыва в дисперсной фазе
    *v_gas_right = v22; // текущая скорость газа справа от контактного разрыва в дисперсной фазе

}

// Метод Годунова интегрирования одномерной системы уравнений Баера-Нунциато при отсутствия градиента объемной доли
// дисперсной фазы, расчет одного шага по времени в одной ячейке для одной фазы
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// center_params_full[M] - вектор примитивных переменных в рассчитываемой ячейке для полной системы (in)
// left_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грании ячейки слева от рассчитываемой для полной системы (in)
// left_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани рассчитываемой ячейки для полной системы (in)
// right_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грани рассчитываемой ячейки для полной системы (in)
// right_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани ячейки справа от рассчитываемой для полной системы (in)
// phase - идентификатор фазы, для которой интегрируется система уравнений - газовая или дисперсная (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// v_ncons_res (out)
// solution_full[M] - вектор примитивных переменных в рассчитываемой ячейке для полной системы на следующем шаге (out)
// n - реальный размер векторов center_params_full, left_minus_ncons_full, left_plus_ncons_full, right_minus_ncons_full, right_plus_ncons_full, solution_full
void godunov_classical_one_phase( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double center_params_full[M],
                                  double left_minus_ncons_full[M], double left_plus_ncons_full[M], double right_minus_ncons_full[M],
                                  double right_plus_ncons_full[M], Phase phase, double dt, double h, double v_ncons_res_left[M], 
                                  double v_ncons_res_right[M], double solution_full[M], int n ) {

    double left_flux_reduced[M_REDUCTION]; // однофазный поток через левую грань ячейки
    double right_flux_reduced[M_REDUCTION]; // однофазный поток через правую грань ячейки
    double left_params_reduced[M_REDUCTION]; // однофазный вектор примитивных переменных в ячейке слева от рассматриваемого разрыва
    double right_params_reduced[M_REDUCTION]; // однофазный вектор примитивных переменных в ячейке справа от рассматриваемого разрыва
    double center_ncons_params_reduced[M_REDUCTION]; // однофазный вектор примитивных переменных в рассчитываемой ячейке
    double center_cons_params_reduced[M_REDUCTION]; // однофазный вектор консервативных переменных в рассчитываемой ячейке
    double solution_ncons_reduced[M_REDUCTION]; // обновленный однофазный вектор примитивных переменных в рассчитываемой ячейке
    double solution_cons_reduced[M_REDUCTION]; // обновленный однофазный вектор консервативных переменных в рассчитываемой ячейке
    int i_comp; // индекс компонента вектора
    double cont_red[M_REDUCTION];
    // расчет однофазного потока через левую грань
    // заполнение полей структуры с отладочной информацией
    debug_info->neighbour_cell = debug_info->current_cell - 1;
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->current_cell_vncons[i_comp] = left_plus_ncons_full[i_comp];
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->neighbour_cell_vncons[i_comp] = left_minus_ncons_full[i_comp];
    // формирование сокращенных однофазных векторов из полных двухфазных
    convert_full_to_reduced( left_minus_ncons_full, phase, left_params_reduced );
    convert_full_to_reduced( left_plus_ncons_full, phase, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, phase, v_ncons_res_left, left_flux_reduced,cont_red );

    // расчет однофазного потока через правую грань
    // заполнение полей структуры с отладочной информацией
    debug_info->neighbour_cell = debug_info->current_cell + 1;
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->current_cell_vncons[i_comp] = right_minus_ncons_full[i_comp];
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->neighbour_cell_vncons[i_comp] = right_plus_ncons_full[i_comp];
    // формирование сокращенных однофазных векторов из полных двухфазных
    convert_full_to_reduced( right_minus_ncons_full, phase, left_params_reduced );
    convert_full_to_reduced( right_plus_ncons_full, phase, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, phase, v_ncons_res_right, right_flux_reduced,cont_red );
    
    // обновление параметров в ячейке
    convert_full_to_reduced( center_params_full, phase, center_ncons_params_reduced );
    convert_noncons_to_cons_reduction( paramsc, center_ncons_params_reduced, phase, center_cons_params_reduced );
    for ( i_comp = 0; i_comp < M_REDUCTION; i_comp++ ) {
        solution_cons_reduced[i_comp] = center_cons_params_reduced[i_comp] - dt * ( right_flux_reduced[i_comp] - left_flux_reduced[i_comp] ) / h;
    }
    convert_cons_to_noncons_reduction( paramsc, solution_cons_reduced, phase, solution_ncons_reduced );

    // запись решения в полный вектор
    convert_reduced_to_full( solution_ncons_reduced, phase, solution_full );
    
}

// Расчет классического однофазного потока С.К. Годунова
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных слева от разрыва (in)
// right_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных справа от разрыва (in)
// phase - идентификатор фазы, для которой рассчитывается поток - газовая или дисперсная (in)
// v_ncons_res[M_REDUCTION] - рассчитанный вектор неконсервативных переменных (in)
// flux[M_REDUCTION] - рассчитанный вектор потока (out)
void godunov_flux_classical( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                             double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons_res[M_REDUCTION], double flux[M_REDUCTION], double cont_ncons_red[M_REDUCTION] ) {
    
    double cl, cr;                      // скорости звука слева и справа от разрыва
    double p_cont, v_cont;              // давление и скорость на контактном разрыве
   
    cl = calc_sound_velocity_reduced( paramsc, left_ncons_params, phase );
    cr = calc_sound_velocity_reduced( paramsc, right_ncons_params, phase );

    // итерационная процедура расчета давления и скорости газа на контактном разрыве
    calc_contact_pressure_velocity( paramsc, debug_info, left_ncons_params, right_ncons_params, M_REDUCTION,
        cl, cr, phase, &p_cont, &v_cont );

    // отбор решения
    sample_reduced( paramsc, left_ncons_params, right_ncons_params, M_REDUCTION, cl, cr, phase,
        p_cont, v_cont, 0.0, v_ncons_res, cont_ncons_red );
    cont_ncons_red[P] = p_cont;
    cont_ncons_red[V] = v_cont;
    // расчет потока Годунова по вектору неконсервативных переменных
    diff_flux_ncons_reduced( paramsc, v_ncons_res, phase, flux );
        
}

// Построение решение классической однофазной задачи Римана
// params - структура с параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных слева от разрыва (in)
// right_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных справа от разрыва (in)
// phase - идентификатор фазы, для которой рассчитывается поток - газовая или дисперсная (in)
// v_ncons[M_REDUCTION] - рассчитанный вектор-решение (out)
void get_classical_Riemann_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                                     double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] ) {

    double cl = calc_sound_velocity_reduced( paramsc, left_ncons_params, phase ); // скорость звука слева от разрыва
    double cr = calc_sound_velocity_reduced( paramsc, right_ncons_params, phase ); // скорость звука справа от разрыва
    double cont_red[M_REDUCTION];
    double p_cont, v_cont; // давление и скорость на контактном разрыве
    // итерационная процедура расчета давления и скорости газа на контактном разрыве
    calc_contact_pressure_velocity( paramsc, debug_info, left_ncons_params, right_ncons_params, M_REDUCTION,
        cl, cr, phase, &p_cont, &v_cont );

    // отбор решения
    sample_reduced( paramsc, left_ncons_params, right_ncons_params, M_REDUCTION, cl, cr, phase,
        p_cont, v_cont, 0.0, v_ncons,cont_red );

}

/* Функция отбора решения в одной из фаз по сокращенному однофазному вектору

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
   + поправка на двучленное уравнение состояния: Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann
   problem and a high-resolution Godunov method for a model of compressible two-phase flow // Journal of Computational Physics.
   - 2006. - V. 212. - P. 490 - 526. - Формулы (9) - (11).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l - вектор неконсервативных переменных слева от разрыва (in)
   v_ncons_r - вектор неконсервативных переменных справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   phase - фаза, для которой отбирается решение - газовая или дисперсная (in)
   p_cont - давление на контактном разрыве (in)
   v_cont - скорость на контактном разрыве (in)
   s - значение x/t, для которого отбирается решение (in)

   v_ncons_res - вектор неконсервативных переменных с отобранными компонентами (out) */
void sample_reduced( struct ParametersCommon *paramsc, double *v_ncons_l, double *v_ncons_r, int vector_size, double cl, double cr,
                     Phase phase, double p_cont, double v_cont, double s, double *v_ncons_res, double cont_ncons_red[M_REDUCTION]  ) {

    double rl, vl, pl;                  /* примитивные переменные слева от разрыва */
    double rr, vr, pr;                  /* примитивные переменные справа от разрыва */
    double g;                           /* показатель адиабаты */
    double p01;                          /* параметр в уравнении состояния */
    double g1, g2, g3, g4, g5, g6, g7;  /* вспомогательные переменные, производные от показателя адиабаты,
                                           в соответствии с Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics.
                                           - 2nd Edition. - Springer, 1999. - P. 153. */

    /* скорости левых волн */
    double shl, stl;    /* скорости "головы" и "хвоста" левой волны разрежения */
    double sl;          /* скорость левой ударной волны */

    /* скорости правых волн */
    double shr, str;    /* скорости "головы" и "хвоста" правой волны разрежения */
    double sr;          /* скорость правой ударной волны */

    double cml, cmr;    /* скорости звука слева и справа от контактного разрыва */
    double c;           /* локальная скорость звука внутри волны разрежения */
    double p_ratio;
    double r, v, p;     /* отобранные значения объемной доли, плотности, скорости и давления */

    /* вспомогательные переменные */
    /* параметры слева от разрыва */
    rl = v_ncons_l[R];
    vl = v_ncons_l[V];
    pl = v_ncons_l[P];
    /* параметры справа от разрыва */
    rr = v_ncons_r[R];
    vr = v_ncons_r[V];
    pr = v_ncons_r[P];
    switch ( phase ) {
        case GAS_PHASE:
            g = paramsc->g2;
            p01 = 0.0;
            break;
        case DISPERSED_PHASE:
            g = paramsc->g1;
            p01 = paramsc->p01;
            break;
        default:
            printf( "\nsample_reduced -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    /* производные от показателя адиабаты */
    g1 = 0.5 * ( g - 1.0 ) / g;
    g2 = 0.5 * ( g + 1.0 ) / g;
    g3 = 2.0 * g / ( g - 1.0 );
    g4 = 2.0 / ( g - 1.0 );
    g5 = 2.0 / ( g + 1.0 );
    g6 = ( g - 1.0 ) / ( g + 1.0 );
    g7 = 0.5 * ( g - 1.0 );

    if ( s <= v_cont ) {
        /* рассматриваемая точка - слева от контактного разрыва */
        if ( p_cont <= pl ) {
            /* левая волна разрежения */
            cont_ncons_red[R] = rl * pow(( p_cont + p01) / (pl + p01), 1.0 / g );
            shl = vl - cl;
            if ( s <= shl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( ( p_cont + p01 ) / ( pl + p01 ), g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* параметры слева от контактного разрыва */
                    r = rl * pow( ( p_cont + p01 ) / ( pl + p01 ), 1.0 / g );
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = ( pl + p01 ) * pow( c / cl, g3 ) - p01;
                }
            }
        }
        else {
            /* левая ударная волна */
            p_ratio = ( p_cont + p01 ) / ( pl + p01 );
            cont_ncons_red[R] = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        /* рассматриваемая точка - справа от контактного разрыва */
        if ( p_cont > pr ) {
            /* правая ударная волна */
            p_ratio = ( p_cont + p01 ) / ( pr + p01 );
            cont_ncons_red[R] = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* правая волна разрежения */
            cont_ncons_red[R] = rr * pow(( p_cont + p01) / (pr + p01), 1.0 / g );
            shr = vr + cr;
            if ( s >= shr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( ( p_cont + p01 ) / ( pr + p01 ), g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* параметры справа от контактного разрыва */
                   r = rr * pow( ( p_cont + p01 ) / ( pr + p01 ), 1.0 / g );
                   v = v_cont;
                   p = p_cont;
               }
               else {
                    /* параметры внутри правой волны разрежения */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = ( pr + p01 ) * pow( c / cr, g3 ) - p01;
               }
            }
        }
    }
    
    /* формирование выходного вектора с результатом */
    v_ncons_res[R] = r;
    v_ncons_res[V] = v;
    v_ncons_res[P] = p;
    
}

/* Функция для построения массива точек, соответствующих множеству состояний, куда можно перевести однофазную среду
   слева или справа от разрыва, пустив по ней ударную волну или волну разрежения

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons - вектор неконсервативных переменных слева или справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   c - скорость звука слева или справа от разрыва (in)
   phase - фаза, для которой строится адиабата - газовая или дисперсная (in)
   dir - направление, для которого строится адиабата - слева или справа от разрыва (in) */
void draw_adiabatic_curve( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *v_ncons, int vector_size, double c, Phase phase, Direction dir ) {

    int i_pt;                       /* счетчик количества точек */
    /* шаг вывода по оси давлений (ординат) */
    double h = ( ADAIABATIC_CURVE_P_T - ADAIABATIC_CURVE_P_B ) / ( ADIABATIC_CURVE_PNUM - 1 );
    double curr_p;                  /* текущее давление */
    double curr_v;                  /* текущая скорость */
    int v_index;                    /* индекс компоненты, соответствующей скорости рассматриваемой фазы */
    double f, fd;                   /* значения функции и производной */
    char fname[MAX_STRING_SIZE];    /* имя файла для записи данных для построения адиабат */
    FILE *adiab_curve;              /* дескриптор файла */    

    /* формирование имени файла и открытие файла для записи данных для построения адиабат */
    strcpy_s( fname, paramsc->output_file_directory );
    if ( dir == LEFT ) {
        strcat_s( fname, "\\adiabatic_curve_left.dat" );
    }
    else {
        strcat_s( fname, "\\adiabatic_curve_right.dat" );
    }
    if ( ( fopen_s( &adiab_curve, fname, "wt" ) ) != 0 ) {
        printf( "\ndraw_adiabatic_curve -> can't open file %s for writing\n\n", fname );
        exit( EXIT_FAILURE );
    }
    
    /* определение индекса компоненты, соответствующей скорости */
    if ( vector_size == M_REDUCTION ) {
        v_index = V;
    }
    else {
        if ( phase == GAS_PHASE )
            v_index = V_GAS;
        else
            v_index = V_DISP;
    }

    for ( i_pt = 0; i_pt < ADIABATIC_CURVE_PNUM; i_pt++ ) {
        curr_p = ADAIABATIC_CURVE_P_B + i_pt * h;   /* текущее значение давления */
        calc_F_and_DF( paramsc, curr_p, v_ncons, vector_size, c, phase, &f, &fd );
        
        if ( dir == LEFT ) curr_v = v_ncons[v_index] - f;   
        else curr_v = v_ncons[v_index] + f;
        
        fprintf_s( adiab_curve, "%e %e\n", curr_v, curr_p );
    }

    fclose( adiab_curve );

}


// Возвращает значение объемной доли твердой фазы, используемое в full_decouple_case_flux (случай NO_GRAD) при расчете потока
// case_beta - номер рассматриваемого варианта ( 0 - берется значение в центре ячейки, 1 - берется значение слева 
//             или справа от ребра, через которое считается поток, в зависимости от знака скорости дисперсной фазы) (in)
// left_edge_beta - значение объемной доли дисперсной фазы слева от ребра, через которое считается поток (in)
// center_beta - значение объемной доли дисперсной фазы в центре рассчитываемой ячейки (in)
// right_edge_beta - значение объемной доли дисперсной фазы справа от ребра, через которое считается поток (in)
// solid_velocity - значение скорости дисперсной фазы, получаемое в результате решения задачи Римана о распаде разрыва (in)
double full_decouple_case_flux_volume_fraction ( int case_beta, double left_edge_beta, double center_beta, double right_edge_beta, double solid_velocity ){

    switch ( case_beta ){
    
    case 0:
        return center_beta;
        break;
    case 1:
        if ( solid_velocity > 0.0 ){
            return left_edge_beta;
            break;
        }
        else{
            return right_edge_beta;
            break;
        }
    default:
        printf( "\nfull_decouple_case_flux_volume_fraction -> wrong case_beta.\n\n" );
        system ( "Pause" );
        exit( EXIT_FAILURE );
    }

}