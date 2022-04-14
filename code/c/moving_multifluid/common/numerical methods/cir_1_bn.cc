// cir_1.cc
// Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунзиато. Первая модификация.
// Куликовский А.Г., Погорелов Н.В., Семенов А.Ю. Математические вопросы численного решения гиперболических систем
// уравнений. - М.: Физматлит, 2001. - С. 62. - Формулы (2.3.18), (2.3.19).
// (c) Уткин Павел, 2013
// Создан: 24 мая 2013 г.

#include "cir_1_bn.h"

#include "utils.h"
#include "math_utils.h"

// Метод Куранта-Изаксона-Рис интегрирования одномерной системы уравнений Баера - Нунзиато, первая модификация
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// left_params[M] - вектор примитивных переменных в ячейке слева от рассчитываемой (in)
// center_params[M] - вектор примитивных переменных в рассчитываемой ячейке (in)
// right_params[M] - вектор примитивных переменных в ячейке справа от рассчитываемой (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге (out)
void cir_1( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M],
            double dt, double h, double solution[M] ) {

    double A[M][M]; // матрица системы
    double S[M][M]; // матрица для аппроксимации вектора-решения на грани ячейки

    // аппроксимации вектора-решения на левой и правой гранях ячейки
    double left_edge_sol[M];
    double right_edge_sol[M];

    // "потоки" через грани ячейки
    double left_flux[M];
    double right_flux[M];
    
    if( !check_matrix_decomposition( paramsc, center_params ) ) {
        printf( "Error: the matrix decomposition is incorrect.\n" );
        exit( EXIT_FAILURE );
    }

    // расчет матрицы для аппроксимации вектора-решения на грани ячейки по параметрам в рассчитываемой ячейке
    cir_util_cir1( paramsc, center_params, S );

    // расчет аппроксимации вектора-решения на гранях ячейки
    calc_edge_solution_cir1( left_params, center_params, S, left_edge_sol );
    calc_edge_solution_cir1( center_params, right_params, S, right_edge_sol );

    // расчет матрицы системы по параметрам в рассчитываемой ячейке
    calc_A_ncons( paramsc, center_params, A );

    // определение "потоков" через левую и правую грани ячейки
    calc_flux_cir1( left_edge_sol, A, left_flux );
    calc_flux_cir1( right_edge_sol, A, right_flux );

    for ( int i = 0; i < M; i++ )
        solution[i] = center_params[i] - dt * ( right_flux[i] - left_flux[i] ) / h;

}

// Расчет матрицы, входящей в аппроксимацию вектора-решения на грани
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// ncons_params[M] - вектор примитивных переменных (in)
// S[M][M] - искомая матрица (out)
void cir_util_cir1( struct ParametersCommon *paramsc, double ncons_params[M], double S[M][M] ) {

    /* матрица из собственных векторов одномерной системы уравнений Баера-Нунзиато */
    double omega_ncons[M][M];

    /* обратная к матрице из собственных векторов одномерной системы уравнений Баера-Нунзиато */
    double omega_ncons_inverse[M][M];
	
    /* диагональная матрица со знаками собственных чисел одномерной системы уравнений Баера-Нунзиато */
    double sign_lambda[M][M];

    /* результат перемножения матриц omega_ncons и sign_lambda */
    double m_tmp[M][M];

    calc_omega_ncons( paramsc, ncons_params, omega_ncons );
    calc_sign_lambda( paramsc, ncons_params, sign_lambda );
    calc_omega_ncons_inverse( paramsc, ncons_params, omega_ncons_inverse );
    if ( !check_inverse_matrix( omega_ncons, omega_ncons_inverse, paramsc->eps_general ) ) {
        printf( "Error: the inverse matrix is incorrect.\n" );
        exit( EXIT_FAILURE );
    }
    mult_matrixes( omega_ncons, sign_lambda, m_tmp, M );
    mult_matrixes( m_tmp, omega_ncons_inverse, S, M );

}

/*  Расчет аппроксимации вектора-решения на грани ячейки
    left_params[M] - вектор примитивных переменных в ячейке слева от грани (in)
    right_params[M] - вектор примитивных переменных в ячейке справа от грани (in)
    S[M][M] - матрица, определяемая по параметрам в рассчитываемой ячейке (in)
    edge_sol[M] - аппроксимация вектора-решения в примитивных переменных на грани (out) */
void calc_edge_solution_cir1( double left_params[M], double right_params[M], double S[M][M], double edge_sol[M] ) {

    int i, j;

    /* расчет аппроксимации вектора-решения на грани */
    for ( i = 0; i < M; i++ ) {
        edge_sol[i] = 0.5 * ( left_params[i] + right_params[i] );
        for ( j = 0; j < M; j++ )
            edge_sol[i] += 0.5 * ( S[i][j] * ( left_params[j] - right_params[j] ) );
    }

}

/*  Расчет "потока" через грань ячейки
    edge_sol[M] - аппроксимация вектора-решения в примитивных переменных на грани (in)
    A_ncons[M][M] - матрица системы в примитивных переменных в рассчитываемой ячейке (in)
    flux[M] - "поток" через грань ячейки (out) */
void calc_flux_cir1( double edge_sol[M], double A[M][M], double flux[M] ) {

    int i, j;

    for ( i = 0; i < M; i++ ) {
        flux[i] = 0.0;
        for ( j = 0; j < M; j++ )
            flux[i] += A[i][j] * edge_sol[j];
    }

}

	