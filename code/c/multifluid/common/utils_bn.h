// utils_bn.h
// Рабочие функции, специфичные для рассматриваемого класса систем уравнений типа Baer-Nunziato и Saurel-Abgrall
// (c) Уткин Павел, 2013 - 2015
// Создан: 16 февраля 2013 г.

#ifndef __UTILS_BN_H_
#define __UTILS_BN_H_

#include "struct.h"
#include "math_utils.h"

// Расчет матрицы одномерной системы уравнений Баера-Нунзиато в "консервативных" переменных
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_cons[M] - вектор "консервативных" переменных (in)
// A_cons[M][M] - матрица одномерной системы уравнений Баера-Нунзиато (out)
void calc_A_cons( struct ParametersCommon *paramsc, double v_cons[M], double A[M][M] );

// Расчет матрицы одномерной системы уравнений Баера-Нунзиато в примитивных переменных
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// A_ncons[M][M] - матрица одномерной системы уравнений Баера-Нунзиато в примитивных переменных (out)
void calc_A_ncons( struct ParametersCommon *paramsc, double v_ncons[M], double A_ncons[M][M] );

void calc_omega_cons( const struct ParametersCommon* paramsc, const double v_cons[M], double omega[M][M] );

// Расчет матрицы из собственных векторов одномерной системы уравнений Баера-Нунзиато в примитивных переменных
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// omega[M][M] - матрица из собственных векторов одномерной системы уравнений Баера-Нунзиато в примитивных переменных (out)
void calc_omega_ncons( struct ParametersCommon *paramsc, double v_ncons[M], double omega[M][M] );

// Расчет матрицы, обратной матрице из собственных векторов одномерной системы уравнений Баера-Нунзиато
// в примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// omega_inverse[M][M] - матрица, обратная матрице из собственных векторов одномерной системы уравнений Баера-Нунзиато
// в примитивных переменных (out)
void calc_omega_ncons_inverse( struct ParametersCommon *paramsc, double v_ncons[M], double omega_inverse[M][M] );

// Расчет диагональной матрицы со знаками собственных чисел системы уравнений Баера-Нунзиато на диагонали
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных (in)
// sign_lambda[M][M] - диагональная матрица со знаками собственных чисел системы уравнений Баера-Нунзиато на диагонали (out)
void calc_sign_lambda( struct ParametersCommon *paramsc, double v_ncons[M], double sign_lambda[M][M] );

void calc_abs_lambda( const struct ParametersCommon* paramsc, const double v_cons[M], double lambda[M][M] );

/* Расчет диагональной матрицы с собственными числами системы уравнений Баера-Нунзиато на диагонали

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons[M] - вектор примитивных переменных (in)

   lambda[M][M] - диагональная матрица с собственными числами системы уравнений Баера-Нунзиато на диагонали (out) */
void calc_lambda( struct ParametersCommon *paramsc, double v_ncons[M], double lambda[M][M] );

// Проверка корректности значений неконсервативных пераметров
// params - структура с основными параметрами вычислительного эксперимента (in) 
// v_ncons[M] - вектор неконсервативных переменных (in)
void check_params_correctness( struct ParametersCommon *paramsc, double v_ncons[M] );

// Проверка, равна ли матрица системы произведению матриц омега, лямбда и омега а минус первой
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор неконсервативных переменных (in)
// Возвращает true, если декомпозияи выполняется верно; false - иначе
bool check_matrix_decomposition( struct ParametersCommon *paramsc, double v_ncons[M] );

#endif // __UTILS_BN_H_