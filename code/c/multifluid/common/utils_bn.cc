// utils_bn.cc
// Рабочие функции для численных методов характеристического типа для решения уравнений БН
// вне зависимости от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 26 августа 2018 г.

#include "utils_bn.h"
#include "utils.h"

// Расчет матрицы одномерной системы уравнений Баера-Нунзиато в "консервативных" переменных
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_cons[M] - вектор "консервативных" переменных (in)
// A_cons[M][M] - матрица одномерной системы уравнений Баера-Нунзиато (out)
void calc_A_cons( struct ParametersCommon *paramsc, double v_cons[M], double A[M][M] ) {

    double v_ncons[M];

    double c1, c2; // скорости звука в дисперсной и газовой фазах
    double v1, v2; // скорости дисперсной и газовой фаз
    double p1, p2; // давления дисперсной и газовой фаз
    double H1, H2; // удельная энтальпия дисперсной и газовой фаз
    double g1, g2; // показатели адибаты дисперсной и газовой фаз без единицы

    convert_cons_to_noncons( paramsc, v_cons, v_ncons, 0 );

    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );

    /* введение обозначений для удобства записи компонент матрицы */
    v1 = v_ncons[V_DISP];
    p1 = v_ncons[P_DISP];
    v2 = v_ncons[V_GAS];
    p2 = v_ncons[P_GAS];
    g1 = paramsc->g1 - 1.0;
    g2 = paramsc->g2 - 1.0;
    H1 = 0.5 * pow( v1, 2.0 ) + pow( c1, 2.0 ) / g1;
    H2 = 0.5 * pow( v2, 2.0 ) + pow( c2, 2.0 ) / g2;
    
    for ( int i = 0; i < M; i++ )
        for ( int j = 0; j < M; j++ )
            A[i][j] = 0.0;

    // первая строка
    A[0][0] = v1;

    // вторая строка
    A[1][2] = 1.0;

    // третья строка
    A[2][0] = - p2 - paramsc->g1 * paramsc->p01;
    A[2][1] = g1 * H1 - pow( v1, 2.0 ) - pow( c1, 2.0 );
    A[2][2] = ( 3.0 - paramsc->g1 ) * v1;
    A[2][3] = g1;

    // четвертая строка
    A[3][0] = - v1 * A[3][0];
    A[3][1] = v1 * ( - H1 + 0.5 * g1 * pow( v1, 2.0 ) );
    A[3][2] = H1 - g1 * pow( v1, 2.0 );
    A[3][3] = v1 * paramsc->g1;

    // пятая строка
    A[4][5] = 1.0;

    // шестая строка
    A[5][0] = p2;
    A[5][4] = g2 * H2 - pow( v2, 2.0 ) - pow( c2, 2.0 );
    A[5][5] = ( 3.0 - paramsc->g2 ) * v2;
    A[5][6] = g2;

    // седьмая строка
    A[6][0] = p2 * v1;
    A[6][4] = v2 * ( - H2 + 0.5 * g2 * pow( v2, 2.0 ) );
    A[6][5] = H2 - g2 * pow( v2, 2.0 );
    A[6][6] = v2 * paramsc->g2;

}

// Расчет матрицы одномерной системы уравнений Баера-Нунзиато в примитивных переменных
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// A_ncons[M][M] - матрица одномерной системы уравнений Баера-Нунзиато в примитивных переменных (out)
void calc_A_ncons( struct ParametersCommon *paramsc, double v_ncons[M], double A_ncons[M][M] ) {

    double b1, b2;      /* объемные доли дисперсной и газовой фаз */
    double c1, c2;      /* скорости звука в дисперсной и газовой фазах */
    double v1, v2;      /* скорости дисперсной и газовой фаз */
    double p1, p2;      /* давления дисперсной и газовой фаз */
    double r1, r2;      /* плотности дисперсной и газовой фаз */

    double c2_sq;       /* квадрат скорости звука в газовой фазе */

    check_params_correctness( paramsc, v_ncons );

    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );
    c2_sq = pow( c2, 2.0 );

    /* введение обозначений для удобства записи компонент матрицы */
    
    b1 = v_ncons[B_DISP];
    r1 = v_ncons[R_DISP];
    v1 = v_ncons[V_DISP];
    p1 = v_ncons[P_DISP];

    b2 = 1.0 - b1;
    r2 = v_ncons[R_GAS];
    v2 = v_ncons[V_GAS];
    p2 = v_ncons[P_GAS];
    
    for ( int i = 0; i < M; i++ )
        for ( int j = 0; j < M; j++ )
            A_ncons[i][j] = 0.0;

    /* первая строка */
    A_ncons[0][0] = v1;

    /* вторая строка */
    A_ncons[1][1] = v1;
    A_ncons[1][2] = r1;

    /* третья строка */
    A_ncons[2][0] = ( p1 - p2 ) / b1 / r1;
    A_ncons[2][2] = v1;
    A_ncons[2][3] = 1 / r1;
    
    /* четвертая строка */
    A_ncons[3][2] = r1 * pow( c1, 2.0 );
    A_ncons[3][3] = v1;

    /* пятая строка */
    A_ncons[4][0] = ( v1 - v2 ) * r2 / b2;
    A_ncons[4][4] = v2;
    A_ncons[4][5] = r2;

    /* шестая строка */
    A_ncons[5][5] = v2;
    A_ncons[5][6] = 1 / r2;

    /* седьмая строка */
    A_ncons[6][0] = A_ncons[4][0] * c2_sq;
    A_ncons[6][5] = r2 * c2_sq;
    A_ncons[6][6] = v2;

}

void calc_omega_cons( const struct ParametersCommon* paramsc, const double v_cons[M], double omega[M][M] ) {

    double v_ncons[M];
    convert_cons_to_noncons( paramsc, v_cons, v_ncons, 0 );

    double c1, c2; // sound velocities in dispersed and gas phases
    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );

    // notation for simplicity
    double v1 = v_ncons[V_DISP];
    double v2 = v_ncons[V_GAS];
    double p2 = v_ncons[P_GAS];
    double g1 = paramsc->g1 - 1.0;
    double g2 = paramsc->g2 - 1.0;
    double H1 = 0.5 * pow( v1, 2.0 ) + pow( c1, 2.0 ) / g1;
    double H2 = 0.5 * pow( v2, 2.0 ) + pow( c2, 2.0 ) / g2;
     
    // the first eigenvector, corresponds to v1
    double A = pow( c2, 2.0 ) - pow( v1 - v2, 2.0 );
    omega[0][0] = - 0.5 * g1 * pow( v1, 2.0 ) * A;
    omega[1][0] = ( p2 + paramsc->g1 * paramsc->p01 ) * A;
    omega[2][0] = v1 * omega[1][0];
    omega[4][0] = 0.5 * paramsc->g2 * g1 * p2 * pow( v1, 2.0 );
    omega[5][0] = v1 * omega[4][0];
    omega[6][0] = 0.5 * p2 * pow( v1, 2.0 ) * g1 / g2 * 
        ( pow( c2, 2.0 ) + g2 * ( pow( v1, 2.0 ) + ( paramsc->g2 - 2.0 ) * v1 * v2 - 0.5 * ( paramsc->g2 - 2.0 ) * pow( v2, 2.0 ) ) );

    // the second eigenvector, corresponds to v1
    omega[1][1] = 1.0;
    omega[2][1] = v1;
    omega[3][1] = 0.5 * pow( v1, 2.0 );
    
    // the third eigenvector, corresponds to v1 - c1
    omega[1][2] = 1.0;
    omega[2][2] = v1 - c1;
    omega[3][2] = H1 - v1 * c1;
    
    // the fourth eigenvector, corresponds to v1 + c1
    omega[1][3] = 1.0;
    omega[2][3] = v1 + c1;
    omega[3][3] = H1 + v1 * c1;

    // the fifth eigenvector, corresponds to v2
    omega[4][4] = 1.0;
    omega[5][4] = v2;
    omega[6][4] = 0.5 * pow( v2, 2.0 );

    // the sixth eigenvector, corresponds to v2 - c2
    omega[4][5] = 1.0;
    omega[5][5] = v2 - c2;
    omega[6][5] = H2 - v2 * c2;

    // the seventh eigenvector, corresponds to v2 + c2
    omega[4][6] = 1.0;
    omega[5][6] = v2 + c2;
    omega[6][6] = H2 + v2 * c2;
}

// Расчет матрицы из собственных векторов одномерной системы уравнений Баера-Нунзиато в примитивных переменных
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// omega[M][M] - матрица из собственных векторов одномерной системы уравнений Баера-Нунзиато в примитивных переменных (out)
void calc_omega_ncons( struct ParametersCommon *paramsc, double v_ncons[M], double omega[M][M] ) {

    double b1, b2;          /* объемные доли дисперсной и газовой фаз */
    double c1, c2;          /* скорости звука в дисперсной и газовой фазах */
    double v1, v2;          /* скорости дисперсной и газовой фаз */
    double p1, p2;          /* давления дисперсной и газовой фаз */
    double r1, r2;          /* плотности дисперсной и газовой фаз */

    double c1_sq, c2_sq;    /* квадраты скоростей звука в дисперсной и газовой фазах */
    double dv;              /* v1 - v2 */
    double A;               /* для промежуточных вычислений */

    check_params_correctness( paramsc, v_ncons );

    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );
    c1_sq = pow( c1, 2.0 );
    c2_sq = pow( c2, 2.0 );

    /* введение обозначений для удобства записи компонент матрицы */
    
    b1 = v_ncons[B_DISP];
    r1 = v_ncons[R_DISP];
    v1 = v_ncons[V_DISP];
    p1 = v_ncons[P_DISP];

    b2 = 1.0 - b1;
    r2 = v_ncons[R_GAS];
    v2 = v_ncons[V_GAS];
    p2 = v_ncons[P_GAS];

    dv = v1 - v2;
    A = b2 * ( c2_sq - pow( dv, 2.0 ) );

    for ( int i = 0; i < M; i++ )
        for ( int j = 0; j < M; j++ )
            omega[i][j] = 0.0;

    /* первый собственный вектор, соответствующий собственному числу v1 */
    omega[0][0] = 1.0;
    omega[3][0] = ( p2 - p1 ) / b1;
    omega[4][0] = - r2 * pow( dv, 2.0 ) / A;
    omega[5][0] = - c2_sq * dv / A;
    omega[6][0] = omega[4][0] * c2_sq;

    /* второй собственный вектор, соответствующий собственному числу v1 */
    omega[1][1] = 1.0;

    /* третий собственный вектор, соответствующий собственному числу v1 - c1 */
    omega[1][2] = r1;
    omega[2][2] = - c1;
    omega[3][2] = r1 * c1_sq;

    /* четвертый собственный вектор, соответствующий собственному числу v1 + c1 */
    omega[1][3] = r1;
    omega[2][3] = c1;
    omega[3][3] = omega[3][2];

    /* пятый собственный вектор, соответствующий собственному числу v2 */
    omega[4][4] = 1.0;

    /* шестой собственный вектор, соответствующий собственному числу v2 - c2 */
    omega[4][5] = r2;
    omega[5][5] = - c2;
    omega[6][5] = r2 * c2_sq;

    /* седьмой собственный вектор, соответствующий собственному числу v2 + c2 */
    omega[4][6] = r2;
    omega[5][6] = c2;
    omega[6][6] = omega[6][5];

}

// Расчет матрицы, обратной матрице из собственных векторов одномерной системы уравнений Баера-Нунзиато
// в примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// omega_inverse[M][M] - матрица, обратная матрице из собственных векторов одномерной системы уравнений Баера-Нунзиато
// в примитивных переменных (out)
void calc_omega_ncons_inverse( struct ParametersCommon *paramsc, double v_ncons[M], double omega_inverse[M][M] ) {

    double b1, b2; // объемные доли дисперсной и газовой фаз
    double c1, c2; // скорости звука в дисперсной и газовой фазах
    double v1, v2; // скорости дисперсной и газовой фаз
    double p1, p2; // давления дисперсной и газовой фаз
    double r1, r2; // плотности дисперсной и газовой фаз

    double c1_sq, c2_sq; // квадраты скоростей звука в дисперсной и газовой фазах
    double dv; // v1 - v2
    double A; // для промежуточных вычислений

    check_params_correctness( paramsc, v_ncons );

    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );
    c1_sq = pow( c1, 2.0 );
    c2_sq = pow( c2, 2.0 );

    // введение обозначений для удобства записи компонент матрицы
    
    b1 = v_ncons[B_DISP];
    r1 = v_ncons[R_DISP];
    v1 = v_ncons[V_DISP];
    p1 = v_ncons[P_DISP];

    b2 = 1.0 - b1;
    r2 = v_ncons[R_GAS];
    v2 = v_ncons[V_GAS];
    p2 = v_ncons[P_GAS];

    dv = v1 - v2;
    A = 0.5 * dv / b2;

    for ( int i = 0; i < M; i++ )
        for ( int j = 0; j < M; j++ )
            omega_inverse[i][j] = 0.0;

    // первая строка
    omega_inverse[0][0] = 1.0;

    // вторая строка
    omega_inverse[1][0] = ( p2 - p1 ) / b1 / c1_sq;
    omega_inverse[1][1] = 1.0;
    omega_inverse[1][3] = - 1.0 / c1_sq;

    // третья строка
    omega_inverse[2][0] = - 0.5 * omega_inverse[1][0] / r1;
    omega_inverse[2][2] = - 0.5 / c1;
    omega_inverse[2][3] = - 0.5 * omega_inverse[1][3] / r1;

    // четвертая строка
    omega_inverse[3][0] = omega_inverse[2][0];
    omega_inverse[3][2] = - omega_inverse[2][2];
    omega_inverse[3][3] = omega_inverse[2][3];

    // пятая строка
    omega_inverse[4][4] = 1.0;
    omega_inverse[4][6] = - 1.0 / c2_sq;

    // шестая строка
    omega_inverse[5][0] = - A / ( c2 + dv );
    omega_inverse[5][5] = - 0.5 / c2;
    omega_inverse[5][6] = - 0.5 * omega_inverse[4][6] / r2;

    // седьмая строка
    omega_inverse[6][0] = A / ( c2 - dv );
    omega_inverse[6][5] =  - omega_inverse[5][5];
    omega_inverse[6][6] = omega_inverse[5][6];

}

// Расчет диагональной матрицы со знаками собственных чисел системы уравнений Баера-Нунзиато на диагонали
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// sign_lambda[M][M] - диагональная матрица со знаками собственных чисел системы уравнений Баера-Нунзиато на диагонали (out)
void calc_sign_lambda( struct ParametersCommon *paramsc, double v_ncons[M], double sign_lambda[M][M] ) {

    double c1, c2;      /* скорости звука в дисперсной и газовой фазах */
    double v1, v2;      /* скорости дисперсной и газовой фаз */

    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );

    /* введение обозначений для удобства записи компонент матрицы */
    v1 = v_ncons[V_DISP];
    v2 = v_ncons[V_GAS];

    for ( int i = 0; i < M; i++ )
        for ( int j = 0; j < M; j++ )
            sign_lambda[i][j] = 0.0;

    sign_lambda[0][0] = sign( v1, paramsc->eps_general );
    sign_lambda[1][1] = sign( v1, paramsc->eps_general );
    sign_lambda[2][2] = sign( v1 - c1, paramsc->eps_general );
    sign_lambda[3][3] = sign( v1 + c1, paramsc->eps_general );
    sign_lambda[4][4] = sign( v2, paramsc->eps_general );
    sign_lambda[5][5] = sign( v2 - c2, paramsc->eps_general );
    sign_lambda[6][6] = sign( v2 + c2, paramsc->eps_general );

}

void calc_abs_lambda( const struct ParametersCommon* paramsc, const double v_cons[M], double lambda[M][M] ) {

    double v_ncons[M];
    convert_cons_to_noncons( paramsc, v_cons, v_ncons, 0 );
    
    double c1, c2; // sound velocities in dispersed and gas phases
    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );

    // notations for simplicity
    double v1 = v_ncons[V_DISP];
    double v2 = v_ncons[V_GAS];

    lambda[0][0] = fabs( v1 );
    lambda[1][1] = fabs( v1 );
    lambda[2][2] = fabs( v1 - c1 );
    lambda[3][3] = fabs( v1 + c1 );
    lambda[4][4] = fabs( v2 );
    lambda[5][5] = fabs( v2 - c2 );
    lambda[6][6] = fabs( v2 +c2 );

}

/* Расчет диагональной матрицы с собственными числами системы уравнений Баера-Нунзиато на диагонали

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons[M] - вектор примитивных переменных (in)

   lambda[M][M] - диагональная матрица с собственными числами системы уравнений Баера-Нунзиато на диагонали (out) */
void calc_lambda( struct ParametersCommon *paramsc, double v_ncons[M], double lambda[M][M] ) {

    double c1, c2; // скорости звука в дисперсной и газовой фазах
    double v1, v2; // скорости дисперсной и газовой фаз

    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );

    // введение обозначений для удобства записи компонент матрицы
    v1 = v_ncons[V_DISP];
    v2 = v_ncons[V_GAS];

    for ( int i = 0; i < M; i++ )
        for ( int j = 0; j < M; j++ )
            lambda[i][j] = 0.0;

    lambda[0][0] = v1;
    lambda[1][1] = v1;
    lambda[2][2] = v1 - c1;
    lambda[3][3] = v1 + c1;
    lambda[4][4] = v2;
    lambda[5][5] = v2 - c2;
    lambda[6][6] = v2 + c2;

}

// Проверка корректности значений неконсервативных пераметров
// params - структура с основными параметрами вычислительного эксперимента (in) 
// v_ncons[M] - вектор неконсервативных переменных (in)
void check_params_correctness( struct ParametersCommon *paramsc, double v_ncons[M] ) {

    double b1, b2; // объемные доли дисперсной и газовой фаз
    double c1, c2; // скорости звука в дисперсной и газовой фазах
    double r1, r2; // плотности дисперсной и газовой фаз
    double v1, v2; // скорости дисперсной и газовой фаз

    double dv; // v1 - v2

    calc_sound_velocity( paramsc, v_ncons, &c1, &c2 );
    
    // введение обозначений для удобства записи
    
    b1 = v_ncons[V_DISP];
    r1 = v_ncons[R_DISP];
    v1 = v_ncons[V_DISP];

    b2 = 1.0 - b1;
    r2 = v_ncons[R_GAS];
    v2 = v_ncons[V_GAS];

    dv = v1 - v2;

    // проверки на то, что вещество в обеих фазах присутствует

    if ( fabs( r1 ) < paramsc->eps_general ) {
        printf( "\ncheck_Parameters1d_correctness -> the density of the dispersed phase is equal to zero.\n" );
        exit( EXIT_FAILURE );
    }
    if ( fabs( b1 ) < paramsc->eps_general ) {
        printf( "\ncheck_Parameters1d_correctness -> the volume fraction of dispersed phase is equal to zero.\n" );
        exit( EXIT_FAILURE );
    }
    if ( fabs( r2 ) < paramsc->eps_general ) {
        printf( "\ncheck_Parameters1d_correctness -> the density of the gaseous phase is equal to zero.\n" );
        exit( EXIT_FAILURE );
    }
    if ( fabs( b2 ) < paramsc->eps_general ) {
        printf( "\ncheck_Parameters1d_correctness -> the volume fraction of the gaseous phase is equal to zero.\n" );
        exit( EXIT_FAILURE );
    }

    // проверка на гиперболичность системы уравнений

    if ( fabs( c2 - dv ) < paramsc->eps_general || fabs( c2 + dv ) < paramsc->eps_general ) {
        printf( "\ncheck_Parameters1d_correctness -> the system of equations is not hyperbolic.\n" );
        exit( EXIT_FAILURE );
    }

    // проверки на странные ошибки

    if ( fabs( c1 ) < paramsc->eps_general ) {
        printf( "\ncheck_Parameters1d_correctness -> the sound velocity for the dispersed phase is equal to zero.\n" );
        exit( EXIT_FAILURE );
    }
    if ( fabs( c2 ) < paramsc->eps_general ) {
        printf( "\ncheck_Parameters1d_correctness -> the sound velocity for the gaseous phase is equal to zero.\n" );
        exit( EXIT_FAILURE );
    }
    
}

// Проверка, равна ли матрица системы произведению матриц омега, лямбда и омега а минус первой
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор неконсервативных переменных (in)
// Возвращает true, если декомпозияи выполняется верно; false - иначе
bool check_matrix_decomposition( struct ParametersCommon *paramsc, double v_ncons[M] ) {

    double A[M][M];             /* матрица системы */
    double omega[M][M];         /* матрица из собственных векторов */
    double omega_inverse[M][M]; /* матрица, обратная матрице из собственных векторов */
    double lambda[M][M];        /* диагональная матрица с собственными числами на диагонали */

    double m1[M][M], m2[M][M];  /* матрицы для хранения промежуточных результатов */

    calc_A_ncons( paramsc, v_ncons, A );
    calc_omega_ncons( paramsc, v_ncons, omega );
    calc_omega_ncons_inverse( paramsc, v_ncons, omega_inverse );
    calc_lambda( paramsc, v_ncons, lambda );
        
    mult_matrixes( omega, lambda, m1, M );
    mult_matrixes( m1, omega_inverse, m2, M );

    for ( int i = 0; i < M; i++ ) {
        for ( int j = 0; j < M; j++ ) {
            if ( fabs( A[i][j] - m2[i][j] ) > paramsc->eps_general )
                return false;
        }
    }

    return true;

}