/*
 * physics_solver.h
 *
 * Учет всей "физики" в ячейке, включая межфазное взаимодействие и межгранулярные напряжения.
 *
 * (c) Уткин Павел, 2014
 *
 * Создан: 7 января 2014 г.
 *
 */

#ifndef __PHYSICS_SOLVER_H
#define __PHYSICS_SOLVER_H

#include "struct.h"

#include "memory.h"

// Решение системы обыкновенных дифференциальных уравнений, описывающих межфазное взаимодействие и иную "физику"
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// dt - текущий шаг по времени (in)
// v_ncons[M] - вектор примитивных переменных в ячейке, на выходе - измененный в результате учета "физики" (in/out)
void physics_solver( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double dt, double *v_ncons, double *v_ncons_next, double *v_ncons_prev, int *number_of_block, double *beta, double *B, double *C, int n, int ignition_flag, int initial_ignition, int cell_number, int left_ghost  ) ;

// Расчет вектора правых частей
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// right_hand_side_terms[M] - текущий вектор правых частей (out)
void calc_right_hand_side_terms( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double *v_ncons, double *right_hand_side_terms, int *number_of_block, int n, int ignition_flag, int initial_ignition ) ;

// Решение СОДУ методом Рунге-Кутты (предиктор-корректор) 2-го порядка
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// v_cons[M] - текущий вектор консервативных переменных в ячейке (in / out)
// right_hand_side_terms[M] - текущий вектор правых частей (in)
// sub_dt - шаг по времени (in)
// i_step_current - номер текущего шага по пространству (in)
void runge_kutta_predictor_2_nd_order (struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double v_cons[M], double right_hand_side_terms[M], double sub_dt, int *number_of_block, int ignition_flag, int initial_ignition );

// Решение полной системы уравнений (7) по неявной схеме Эйлера методом Ньютона - задача 1d2phc
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами 1dphc (in)
// v_cons[M1D] - текущий вектор консервативных переменных в ячейке (in / out)
// sub_dt - шаг по времени (in)
// number_of_block - номер блока (in)
// n - реальный размер векторов (in)
void impicit_euler_full_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double v_cons[M], double sub_dt, int number_of_block, int n );

// Решение полной системы уравнений (9) по неявной схеме Эйлера методом Ньютона - задача 2d2phc
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// params2d - структура с дополнительными параметрами 2d2phc (in)
// v_cons[M] - текущий вектор консервативных переменных в ячейке (in / out)
// sub_dt - шаг по времени (in)
// number_of_block - номер блока (in)
// n - реальный размер векторов
void impicit_euler_full_2d( struct ParametersCommon *paramsc, struct Parameters2d *params2d, double v_cons[M], double sub_dt, int *number_of_block, int n );

#endif /* __PHYSICS_SOLVER_H */