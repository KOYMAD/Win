// hll.cc
// Метод HLL из Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) Уткин Павел, 2015
// Создан: 8 июля 2015 г.

#include "hll.h"

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
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure, double S_center, double S_left, double S_right, int l ) {
    double solution_ncons_Lh[M]; // вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
if (l == 1049)
    printf("\n %lf", center_ncons[P_GAS]);
    Lh_HLL( paramsc, left_ncons, center_ncons, right_ncons, slopes_left, slopes_center, slopes_right,
        dt, h, solution_ncons_Lh, n, number_of_scalars, S_center, S_left, S_right );
    if (l == 1049)
    printf("\n %lf", solution_ncons_Lh[P_GAS]);// действие гиперболического оператора
    if ( paramsc->pressure_relaxation == true && is_pressure_relaxation_after_this_step == true ){
        Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // действие релаксационного оператора
        if (l == 1049)
    printf("\n %lf", solution_ncons[P_GAS]);
    }
    else {
        for ( int i = 0; i < n; i++ )
            solution_ncons[i] = solution_ncons_Lh[i];
    }

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
void Lh_HLL( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars, double S_center, double S_left, double S_right  ) {

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

    // расчет потока через левую границу ячейки

    array1D left_flux( M );
    double splus_l, sminus_l; // расчет оценок для скоростей волн S+ и S- для полной системы на левом ребре
    hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars);

    // расчет потока через правую границу ячейки
    array1D right_flux( M );
    double splus_r, sminus_r; // расчет оценок для скоростей волн S+ и S- для полной системы на правом ребре
    hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
    //printf( "\n %lf %lf", left_flux[P_GAS],right_flux[P_GAS]);
    // аппроксимация неконсервативных составляющих в правых частях
    array1D rhst( M ); // неконсервативный вектор правых частей
    array1D additional( M ); // добавка к правой части для переменного сечения
    rhst_ncons(h, splus_l, sminus_l, splus_r, sminus_r, paramsc, center_ncons, &rhst, &additional, number_of_scalars, S_center, S_left, S_right );
    double dB2dx; // аппроксимация градиента объемной доли 
    calc_dB2( left_minus_cons, left_plus_cons, right_minus_cons, right_plus_cons, splus_l, sminus_l, splus_r, sminus_r, &dB2dx, S_center, S_left, S_right);
    // рассчитываем компоненты вектора-решения на следующем шаге по времени без первой компоненты (объемная доля дисперсной фазы)
    double center_cons[M]; // вектор консервативных переменных в текущей рассчитываемой ячейке
    convert_noncons_to_cons( paramsc, center_ncons, center_cons, number_of_scalars );
    double solution_cons[M]; // вектор-решение консервативных переменных на следующем шаге по времени
    for ( int i = 1; i < n; i++ ){
        solution_cons[i] = center_cons[i] - dt * ( right_flux[i] - left_flux[i] ) / h + dt *( dB2dx / h * rhst[i] + additional[i] );
    }

 
    // расчет переноса объемной доли дисперсной фазы
    double u_i = calc_u_i( paramsc, center_ncons );
    double t1 = ( u_i * ( splus_r * right_minus_cons[B_DISP] - sminus_r * right_plus_cons[B_DISP] ) + splus_r * sminus_r * ( right_plus_cons[B_DISP] - right_minus_cons[B_DISP] ) )
        / ( splus_r - sminus_r );
    double t2 = ( u_i * ( splus_l * left_minus_cons[B_DISP] - sminus_l * left_plus_cons[B_DISP] ) + splus_l * sminus_l * ( left_plus_cons[B_DISP] - left_minus_cons[B_DISP] ) )
        / ( splus_l - sminus_l );        
        solution_cons[B_DISP] = center_ncons[B_DISP] - dt *( t1 - t2 ) / h ;
        if (solution_cons[B_DISP] > 0.94)
            solution_cons[B_DISP] = 0.94;
  
    convert_cons_to_noncons( paramsc, solution_cons, solution_ncons, number_of_scalars);

   
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
    diff_flux_ncons( paramsc, left_ncons, &left_diff_flux, number_of_scalars);
    array1D right_diff_flux( M );
    diff_flux_ncons( paramsc, right_ncons, &right_diff_flux, number_of_scalars );
    //printf( "\n %lf %lf", left_diff_flux[P_GAS],right_diff_flux[P_GAS]);
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
    convert_noncons_to_cons( paramsc, left_ncons, left_cons, number_of_scalars);
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

// убрать производную по S -- done
void calc_dB2( const double left_minus[M], const double left_plus[M], const double right_minus[M], const double right_plus[M],
               double splus_l, double sminus_l, double splus_r, double sminus_r, double* dB2dx, double S_center, double S_left, double S_right) {

    *dB2dx = ( ( splus_r * right_minus[B_DISP] - sminus_r * right_plus[B_DISP] ) / ( splus_r - sminus_r )
        - ( splus_l * left_minus[B_DISP] - sminus_l * left_plus[B_DISP] ) / ( splus_l - sminus_l ) );



}

