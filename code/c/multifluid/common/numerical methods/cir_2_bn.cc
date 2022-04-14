/*
 * cir_2.cc
 *
 * Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунзиато. Вторая модификация.
 *
 * Куликовский А.Г., Погорелов Н.В., Семенов А.Ю. Математические вопросы численного решения гиперболических систем
 * уравнений. - М.: Физматлит, 2001. - С. 62. - Формулы (2.3.22), (2.3.23) + аппроксимация матрицы S на грани как 0.5 * ( S_left + S_right ).
 *
 * (c) Уткин Павел, 2013
 *
 * Создан: 24 мая 2013 г.
 *
 */

#include "cir_2_bn.h"

#include "utils.h"
#include "math_utils.h"

// Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунзиато, вторая модификация
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// left_params[M] - вектор примитивных переменных в ячейке слева от рассчитываемой (in)
// center_params[M] - вектор примитивных переменных в рассчитываемой ячейке (in)
// right_params[M] - вектор примитивных переменных в ячейке справа от рассчитываемой (in)
// dt - временной шаг (in)
//  h - пространственный шаг (in)
//  solution[M] - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге (out)
void cir_2( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M], double dt, double h, double solution[M] ) {

    /* матрица системы */
    double A[M][M];

    /* аппроксимации вектора-решения на левой и правой гранях ячейки */
    double left_edge_sol[M];
    double right_edge_sol[M];

    /* "потоки" через грани ячейки */
    double left_flux[M];
    double right_flux[M];
    
    if( !check_matrix_decomposition( paramsc, center_params ) ) {
        printf( "Error: the matrix decomposition is incorrect.\n" );
        exit( EXIT_FAILURE );
    }

    /* расчет аппроксимации вектора-решения на гранях ячейки */
    calc_edge_solution_cir2( paramsc, left_params, center_params, left_edge_sol );
    calc_edge_solution_cir2( paramsc, center_params, right_params, right_edge_sol );

    /* расчет матрицы системы по параметрам в рассчитываемой ячейке */
    calc_A_ncons( paramsc, center_params, A );

    /* определение "потоков" через левую и правую грани ячейки */
    calc_flux_cir2( left_edge_sol, A, left_flux );
    calc_flux_cir2( right_edge_sol, A, right_flux );

    for ( int i = 0; i < M; i++ )
        solution[i] = center_params[i] - dt * ( right_flux[i] - left_flux[i] ) / h;

}

// Расчет аппроксимации вектора-решения на грани ячейки
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// left_params[M] - вектор примитивных переменных в ячейке слева от грани (in)
// right_params[M] - вектор примитивных переменных в ячейке справа от грани (in)
// edge_sol[M] - аппроксимация вектора-решения в примитивных переменных на грани (out)
void calc_edge_solution_cir2( struct ParametersCommon *paramsc, double left_params[M], double right_params[M], double edge_sol[M] ) {

    // результат перемножения всех матриц для параметров слева и справа от грани
    double m_left[M][M];
    double m_right[M][M];

    // расчет матриц, входящих в аппроксимацию вектора-решения на грани
    cir_util_cir2( paramsc, left_params, m_left );
    cir_util_cir2( paramsc, right_params, m_right );

    // расчет аппроксимации вектора-решения на грани
    for ( int i = 0; i < M; i++ ) {
        edge_sol[i] = 0.5 * ( left_params[i] + right_params[i] );
        for ( int j = 0; j < M; j++ )
            edge_sol[i] += 0.5 * ( 0.5 * ( m_left[i][j] + m_right[i][j] ) * ( left_params[j] - right_params[j] ) );
    }

}

// Расчет матрицы, входящей в аппроксимацию вектора-решения на грани
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// ncons_params[M] - вектор примитивных переменных (in)
// m[M][M] - искомая матрица (out)
void cir_util_cir2( struct ParametersCommon *paramsc, double ncons_params[M], double m[M][M] ) {

    double omega_ncons[M][M]; // матрица из собственных векторов одномерной системы уравнений Баера-Нунзиато
    double omega_ncons_inverse[M][M]; // обратная к матрице из собственных векторов одномерной системы уравнений Баера-Нунзиато
    double sign_lambda[M][M]; // диагональная матрица со знаками собственных чисел одномерной системы уравнений Баера-Нунзиато
    double m_tmp[M][M]; // результат перемножения матриц omega_ncons и sign_lambda

    calc_omega_ncons( paramsc, ncons_params, omega_ncons );
    calc_sign_lambda( paramsc, ncons_params, sign_lambda );
    calc_omega_ncons_inverse( paramsc, ncons_params, omega_ncons_inverse );
    if ( !check_inverse_matrix( omega_ncons, omega_ncons_inverse, paramsc->eps_general ) ) {
        printf( "Error: the inverse matrix is incorrect.\n" );
        exit( EXIT_FAILURE );
    }
    mult_matrixes( omega_ncons, sign_lambda, m_tmp, M );
    mult_matrixes( m_tmp, omega_ncons_inverse, m, M );

}

/*  Расчет "потока" через грань ячейки
    edge_sol[M] - аппроксимация вектора-решения в примитивных переменных на грани (in)
    A_ncons[M][M] - матрица системы в примитивных переменных в рассчитываемой ячейке (in)
    flux[M] - "поток" через грань ячейки (out) */
void calc_flux_cir2( double edge_sol[M], double A[M][M], double flux[M] ) {

    int i, j;

    for ( i = 0; i < M; i++ ) {
        flux[i] = 0.0;
        for ( j = 0; j < M; j++ )
            flux[i] += A[i][j] * edge_sol[j];
    }

}

	