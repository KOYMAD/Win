/*
 * cir_3.cc
 *
 * Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунзиато. Третья модификация.
 *
 * Куликовский А.Г., Погорелов Н.В., Семенов А.Ю. Математические вопросы численного решения гиперболических систем
 * уравнений. - М.: Физматлит, 2001. - С. 62. - Формулы (2.3.22), (2.3.23) + аппроксимация матрицы S на грани как S ( 0.5 * ( u_left + u_right ) ).
 *
 * (c) Уткин Павел, 2013
 *
 * Создан: 31 мая 2013 г.
 *
 */

#include "cir_3_bn.h"

#include "utils.h"
#include "math_utils.h"

// Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунзиато, третья модификация
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// left_params[M] - вектор примитивных переменных в ячейке слева от рассчитываемой (in)
// center_params[M] - вектор примитивных переменных в рассчитываемой ячейке (in)
// right_params[M] - вектор примитивных переменных в ячейке справа от рассчитываемой (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге (out)
void cir_3( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M], double dt, double h, double solution[M] ) {

    double A[M][M]; // матрица системы

    // аппроксимации вектора-решения на левой и правой гранях ячейки
    double left_edge_sol[M];
    double right_edge_sol[M];

    // "потоки" через грани ячейки
    double left_flux[M];
    double right_flux[M];
    
    if( !check_matrix_decomposition( paramsc, center_params ) ) {
        printf( "\ncir_3 -> matrix decomposition is incorrect\n\n" );
        exit( EXIT_FAILURE );
    }

    // расчет аппроксимации вектора-решения на гранях ячейки
    calc_edge_solution_cir3( paramsc, left_params, center_params, left_edge_sol );
    calc_edge_solution_cir3( paramsc, center_params, right_params, right_edge_sol );

    calc_A_ncons( paramsc, center_params, A ); // расчет матрицы системы по параметрам в рассчитываемой ячейке

    // определение "потоков" через левую и правую грани ячейки
    calc_flux_cir3( left_edge_sol, A, left_flux );
    calc_flux_cir3( right_edge_sol, A, right_flux );

    for ( int i = 0; i < M; i++ )
        solution[i] = center_params[i] - dt * ( right_flux[i] - left_flux[i] ) / h;

}

// Расчет аппроксимации вектора-решения на грани ячейки
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// left_params[M] - вектор примитивных переменных в ячейке слева от грани (in)
// right_params[M] - вектор примитивных переменных в ячейке справа от грани (in)
// edge_sol[M] - аппроксимация вектора-решения в примитивных переменных на грани (out)
void calc_edge_solution_cir3( struct ParametersCommon *paramsc, double left_params[M], double right_params[M], double edge_sol[M] ) {

    double S[M][M]; // матрица, входящая в аппроксимацию вектора-решения на грани
    double edge_params[M]; // интерполяция примитивных переменных на грань

    for ( int i = 0; i < M; i++ )
        edge_params[i] = 0.5 * ( left_params[i] + right_params[i] );

    cir_util_cir3( paramsc, edge_params, S ); // расчет матрицы, входящей в аппроксимацию вектора-решения на грани

    // расчет аппроксимации вектора-решения на грани
    for ( int i = 0; i < M; i++ ) {
        edge_sol[i] = 0.5 * ( left_params[i] + right_params[i] );
        for ( int j = 0; j < M; j++ )
            edge_sol[i] += 0.5 * S[i][j] * ( left_params[j] - right_params[j] );
    }

}

// Расчет матрицы, входящей в аппроксимацию вектора-решения на грани
// params - структура с параметрами вычислительного эксперимента (in)
// ncons_params[M] - вектор примитивных переменных (in)
// m[M][M] - искомая матрица (out)
void cir_util_cir3( struct ParametersCommon *paramsc, double ncons_params[M], double m[M][M] ) {

    double omega_ncons[M][M]; // матрица из собственных векторов одномерной системы уравнений БН
    double omega_ncons_inverse[M][M]; // обратная к матрице из собственных векторов одномерной системы уравнений БН
    double sign_lambda[M][M]; // диагональная матрица со знаками собственных чисел одномерной системы уравнений БН
    double m_tmp[M][M]; // результат перемножения матриц omega_ncons и sign_lambda

    calc_omega_ncons( paramsc, ncons_params, omega_ncons );
    calc_sign_lambda( paramsc, ncons_params, sign_lambda );
    calc_omega_ncons_inverse( paramsc, ncons_params, omega_ncons_inverse );
    if ( !check_inverse_matrix( omega_ncons, omega_ncons_inverse, paramsc->eps_general ) ) {
        printf( "Error: the inverse matrix is incorrect.\n" );
        exit( EXIT_FAILURE );
    }
    mult_matrixes( omega_ncons, sign_lambda, m_tmp , M);
    mult_matrixes( m_tmp, omega_ncons_inverse, m, M );

}

/*  Расчет "потока" через грань ячейки
    edge_sol[M] - аппроксимация вектора-решения в примитивных переменных на грани (in)
    A_ncons[M][M] - матрица системы в примитивных переменных в рассчитываемой ячейке (in)
    flux[M] - "поток" через грань ячейки (out) */
void calc_flux_cir3( double edge_sol[M], double A[M][M], double flux[M] ) {

    int i, j;

    for ( i = 0; i < M; i++ ) {
        flux[i] = 0.0;
        for ( j = 0; j < M; j++ )
            flux[i] += A[i][j] * edge_sol[j];
    }

}

	