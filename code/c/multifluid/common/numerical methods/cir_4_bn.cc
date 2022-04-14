// cir_4.cc
// Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунциато. Четвертая модификация.
// Куликовский А.Г., Погорелов Н.В., Семенов А.Ю. Математические вопросы численного решения гиперболических систем
// уравнений. - М.: Физматлит, 2001. - С. 67. - Формулы (2.3.44), (2.3.45).
// (c) Уткин Павел, 2014
// Создан: 08 августа 2014 г.

#include "cir_4_bn.h"

#include "utils.h"
#include "math_utils.h"

// Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунциато, четвертая модификация
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_params - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_params - вектор примитивных переменных в рассчитываемой ячейке
// right_params - вектор примитивных переменных в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге
void cir_4( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M],
            double dt, double h, double solution[M] ) {

    // "потоки" через грани ячейки
    array1D left_flux( M );
    array1D right_flux( M );

    // вклады неконсервативных составляющих на левой и правой гранях ячейки
    array1D left_rhst( M );
    array1D right_rhst( M );

    // переводим все в консервативные переменные
    double left_cons_params[M];
    convert_noncons_to_cons( paramsc, left_params, left_cons_params, 0 );
    double center_cons_params[M];
    convert_noncons_to_cons( paramsc, center_params, center_cons_params, 0 );
    double right_cons_params[M];
    convert_noncons_to_cons( paramsc, right_params, right_cons_params, 0 );

    // расчет "потоков" через левую и правую грани ячейки
    calc_flux_cir4( paramsc, left_cons_params, center_cons_params, &left_flux );
    calc_flux_cir4( paramsc, center_cons_params, right_cons_params, &right_flux );

    // расчет вклада неконсервативных составляющих на левой и правой гранях
    calc_rhst_term_cir4( paramsc, left_cons_params, center_cons_params, &left_rhst );
    calc_rhst_term_cir4( paramsc, center_cons_params, right_cons_params, &right_rhst );

    double solution_cons[M];
    // обновление вектора "консервативных" переменных
    for ( int iComp = 0; iComp < M; iComp++ )
        solution_cons[iComp] = center_cons_params[iComp] - dt * ( right_flux[iComp] - left_flux[iComp] ) / h +
            0.5 * ( left_rhst[iComp] + right_rhst[iComp] ) / h;

    convert_cons_to_noncons( paramsc, solution_cons, solution, 0 );
}

// Расчет "потока" через грань ячейки
// params - структура с основными параметрами вычислительного эксперимента
// left_params - вектор "консервативных" переменных слева от разрыва
// right_params - вектор "консервативных" переменных справа от разрыва
// flux - "поток" через грань ячейки
void calc_flux_cir4( const struct ParametersCommon* paramsc, const double left_params[M], const double right_params[M], array1D* flux ) {

    vector<double> aver_flux( M ); // усредненный дифференциальный "поток"
    vector<double> left_diff_flux( M ); // вектор дифференциального "потока" по параметрам слева от разрыва
    vector<double> right_diff_flux( M ); // вектор дифференциального "потока" по параметрам справа от разрыва
    vector<double> charact_part( M ); // характеристическая составляющая "потока"
    double aver_params[M]; // усредненные параметры на ребре
    double aver_ncons_params[M]; // примитивный вектор усредненных параметров на ребре

    // расчет усредненных параметров на ребре
    for ( int iComp = 0; iComp < M; iComp++ )
        aver_params[iComp] = 0.5 * ( left_params[iComp] + right_params[iComp] );
    
    // расчет первой, дифференциальной составляющей "потока"
    diff_flux_cons( paramsc, left_params, &left_diff_flux, 0 );
    diff_flux_cons( paramsc, right_params, &right_diff_flux, 0 );
    for ( int iComp = 0; iComp < M; iComp++ )   
        aver_flux[iComp] = 0.5 * ( left_diff_flux[iComp] + right_diff_flux[iComp] );

    // расчет второй, характеристической составляющей
    // расчет максимального по модулю собственного числа матрицы Якоби
    convert_cons_to_noncons( paramsc, aver_params, aver_ncons_params, 0 );
    double c1, c2;
    calc_sound_velocity(paramsc, aver_ncons_params, &c1, &c2 );
    double eigenvalues[M-1] = { fabs( aver_ncons_params[V_DISP] ), fabs( aver_ncons_params[V_GAS] ),
                                fabs( aver_ncons_params[V_DISP] + c1 ), fabs( aver_ncons_params[V_DISP] - c1 ),
                                fabs( aver_ncons_params[V_GAS] + c2 ), fabs( aver_ncons_params[V_GAS] - c2 ) };
    double max_eigenvalue = 0.0;
    for ( int iComp = 0; iComp < M - 1; iComp++ ) {
        if ( eigenvalues[iComp] > max_eigenvalue )
            max_eigenvalue = eigenvalues[iComp];
    }
    for ( int iComp = 0; iComp < M; iComp++ )
        charact_part[iComp] = 0.5 * max_eigenvalue * ( left_params[iComp] - right_params[iComp] );

    for ( int iComp = 0; iComp < M; iComp++ )
        (*flux)[iComp] = aver_flux[iComp] + charact_part[iComp];
}

// Расчет вклада неконсервативных членов в правой части для одной из граней расчетной ячейки
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_params - вектор "консервативных" переменных слева от разрыва
// right_params - вектор "консервативных" переменных справа от разрыва
// rhst - искомый вклад неконсервативных членов в правой части
void calc_rhst_term_cir4( const struct ParametersCommon* paramsc, const double left_params[M], const double right_params[M], array1D* rhst ) {

    double aver_params[M]; // усредненные параметры на ребре

    // расчет усредненных параметров на ребре
    for ( int iComp = 0; iComp < M; iComp++ )
        aver_params[iComp] = 0.5 * ( left_params[iComp] + right_params[iComp] );

    // расчет неконсервативного вектора правых частей
    rhst_cons( paramsc, aver_params, rhst, 0 );

    // домножили на разность объемной доли дисперсной фазы
    for ( int iComp = 0; iComp < M; iComp++ )
        (*rhst)[iComp] *= right_params[B_DISP] - left_params[B_DISP];

}