// utils.cc
// Рабочие функции, специфичные для рассматриваемого класса систем уравнений типа Baer-Nunziato и Saurel-Abgrall
// (c) Уткин Павел, 2013 - 2015
// Создан: 16 февраля 2013 г.

#ifndef __UTILS_1D2PHC_H_
#define __UTILS_1D2PHC_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "math_utils.h"
#include "io_1d.h"
#include "utils.h"


// Расчет шага интегрирования по времени
// paramsc - структура с основными параметрами вычислительного эксперимента
// params1d - структура с параметрами одномерной задачи
// *x - координаты узлов сетки
// **v_ncons - вектора примитивных переменных в центрах ячеек
// Возвращает шаг интегрирования по времени
double get_time_step( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, const double* x, double **v_ncons, int print ) ;

// Проверка сохранения массы в системе
// params1d - структура с параметрами одномерной задачи
// *x - координаты узлов сетки
// **v_ncons - вектора примитивных переменных в центрах ячеек
// initial_total_mass - исходная суммарная масса вещества в системе
// mass_diff - относительная погрешность суммарной массы вещества в системе относительно исходной
void check_mass( const struct Parameters1d* params1d, const double* x, double** v_ncons, const double initial_total_mass, double* mass_diff );

// Расчет полной массы вещества в системе
// params1d - структура с параметрами одномерной задачи
// *x - координаты узлов сетки
// **v_cons - вектора "консервативных" переменных в центрах ячеек
// Возвращает суммарную массу газа и дисперсной фазы в системе
double get_total_mass( const struct Parameters1d* params1d, const double* x, double** v_cons );

void calculate_DCJ(struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *Dcj);


#endif // __UTILS_1D2PHC_H_