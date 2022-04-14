// rusanov_sa_1d2phc.cc
// Метод Русанова численного интегрирования уравнений Saurel-Abgrall, 1D случай
// Реализовано по: Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) Уткин Павел, 2018
// Создан: 4 марта 2018 г.

#include "rusanov_sa_1d2phc.h"

// Метод Русанова численного интегрирования уравнений Saurel-Abgrall, 1D случай
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
// configuration_pressure - конфигурационное давление
void rusanov_1d( struct ParametersCommon* paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
                 double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                 double dt, double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure ) {

    double solution_ncons_Lh[M]; // вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора

    Lh_rusanov_1d( paramsc, left_ncons, center_ncons, right_ncons, slopes_left, slopes_center, slopes_right,
        dt, h, solution_ncons_Lh ); // действие гиперболического оператора
    if ( paramsc->pressure_relaxation == true )
        Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_pressure  ); // действие релаксационного оператора
    else {
        for ( int i = 0; i < n; i++ )
            solution_ncons[i] = solution_ncons_Lh[i];
    }

}

// Гиперболический оператор метода Русанова, 1D случай
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
void Lh_rusanov_1d( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
                    const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                    const double dt, const double h, double solution_ncons[M] ) {

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
        convert_noncons_to_cons( paramsc, left_ncons, left_minus_cons, 0 );
        left_minus_cons[j] += 0.5 * h * slopes_left[j];
        convert_cons_to_noncons( paramsc, left_minus_cons, left_minus_ncons, 0 );
        // реконструкция значений на левой границе справа
        convert_noncons_to_cons( paramsc, center_ncons, left_plus_cons, 0 );
        left_plus_cons[j] -= 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, left_plus_cons, left_plus_ncons, 0 );
        // реконструкция значений на правой границе слева
        convert_noncons_to_cons( paramsc, center_ncons, right_minus_cons, 0 );
        right_minus_cons[j] += 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, right_minus_cons, right_minus_ncons, 0 );
        // реконструкция значений на правой границе справа
        convert_noncons_to_cons( paramsc, right_ncons, right_plus_cons, 0 );
        right_plus_cons[j] -= 0.5 * h * slopes_right[j];
        convert_cons_to_noncons( paramsc, right_plus_cons, right_plus_ncons, 0 );
    }

    // расчет потока через левую границу ячейки
    array1D left_flux( M );
    double splus_left;
    rusanov_flux_1d( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_left );

    // расчет потока через правую границу ячейки
    array1D right_flux( M );
    double splus_right;
    rusanov_flux_1d( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_right );

    // аппроксимация неконсервативных составляющих в правых частях
    array1D rhst( M ); // неконсервативный вектор правых частей
    rhst_ncons( paramsc, center_ncons, &rhst, 0 );
    double dB1dx = 0.5 * ( right_plus_ncons[B_DISP] - left_minus_ncons[B_DISP] ) / h; // аппроксимация градиента объемной доли дисперсной фазы
 
    // рассчитываем компоненты вектора-решения на следующем шаге по времени без первой компоненты (объемная доля дисперсной фазы)
    double center_cons[M]; // вектор консервативных переменных в текущей рассчитываемой ячейке
    convert_noncons_to_cons( paramsc, center_ncons, center_cons, 0 );
    double solution_cons[M]; // вектор-решение консервативных переменных на следующем шаге по времени
    for ( int i = 1; i < M; i++ )
        solution_cons[i] = center_cons[i] - dt * ( right_flux[i] - left_flux[i] ) / h + dt * dB1dx * rhst[i];

    // расчет переноса объемной доли дисперсной фазы
    double u_i = calc_u_i( paramsc, center_ncons );
    double t1 = - 0.5 * u_i * ( right_plus_ncons[B_DISP] - left_minus_ncons[B_DISP] );
    double t2 = 0.5 * ( splus_right * ( right_plus_ncons[B_DISP] - right_minus_ncons[B_DISP] ) - splus_left * ( left_plus_ncons[B_DISP] - left_minus_ncons[B_DISP] ) );
    solution_cons[B_DISP] = center_ncons[B_DISP] + dt * ( t1 + t2 ) / h;
    
    convert_cons_to_noncons( paramsc, solution_cons, solution_ncons, 0 );

}

// Метод Русанова расчета потоков в двухфазной среде, 1D случай
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных слева от разрыва
// right_ncons - вектор примитивных переменных справа от разрыва
// flux - рассчитанный вектор потока
// smax - оценка для скорости S+
void rusanov_flux_1d( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M], array1D* flux, double *splus ) {

    // расчет векторов дифференциального потока по параметрам слева и справа от разрыва
    array1D left_diff_flux( M );
    diff_flux_ncons( paramsc, left_ncons, &left_diff_flux, 0 );
    array1D right_diff_flux( M );
    diff_flux_ncons( paramsc, right_ncons, &right_diff_flux, 0 );

    // расчет оценки для скорости волны S+ с учетом обеих фаз
    // терминология соответствует Торо (издание 2009 г., стр. 329)
    calc_splus_rusanov( paramsc, left_ncons, right_ncons, splus );

    // расчет векторов консервативных переменных слева и справа от разрыва
    double left_cons[M];
    convert_noncons_to_cons( paramsc, left_ncons, left_cons, 0 );
    double right_cons[M];
    convert_noncons_to_cons( paramsc, right_ncons, right_cons, 0 );
	
    for ( int i = 1; i < M; i++ )
        (*flux)[i] = 0.5 * ( left_diff_flux[i] + right_diff_flux[i] - (*splus) * ( right_cons[i] - left_cons[i] ) );

}

// Расчет оценки для скорости волны S+ с учетом обеих фаз
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons[M] - вектор примитивных переменных слева от разрыва
// right_ncons[M] - вектор примитивных переменных справа от разрыва
// splus - оценка для скорости волны S+
void calc_splus_rusanov( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M], double* splus ) {

    double c_gas_left, c_disp_left; 
    calc_sound_velocity( paramsc, left_ncons, &c_gas_left, &c_disp_left ); // скорость звука слева от разрыва для обеих фаз
    double c_gas_right, c_disp_right; 
    calc_sound_velocity( paramsc, right_ncons, &c_gas_right, &c_disp_right ); // скорость звука справа от разрыва для обеих фаз

    double t1 = max( max( max( fabs( left_ncons[V_GAS] + c_gas_left ), fabs( right_ncons[V_GAS] + c_gas_right ) ), fabs( left_ncons[V_DISP] + c_disp_left ) ),
        fabs( right_ncons[V_DISP] + c_disp_right ) );
    double t2 = max( max( max( fabs( left_ncons[V_GAS] - c_gas_left ), fabs( right_ncons[V_GAS] - c_gas_right ) ), fabs( left_ncons[V_DISP] - c_disp_left ) ),
        fabs( right_ncons[V_DISP] - c_disp_right ) );

    *splus = max( t1, t2 );

}