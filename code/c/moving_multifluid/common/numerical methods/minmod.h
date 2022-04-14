// minmod.h
// Повышение порядка аппроксимации метода за счет кусочно-линейного восполнения сеточных консервативных переменных
// с использованием ограничителя minmod.
// (c) Уткин Павел, 2016
// Создан: 22 августа 2016 г.

#ifndef __MINMOD_H_
#define __MINMOD_H_

#include "utils.h"
#include "math_utils.h"

// Реконструкция консервативных переменных с использованием ограничителя minmod
// paramsc - структура с основными параметрами вычислительного эксперимента
// params1d - структура с параметрами одномерной задачи
// u_prev - вектора примитивных переменных во всех ячейках расчетной области
// xc - массив координат центров ячеек сетки
// slopes - вектора наклонов во всех внутренних ячейках расчетной области
void reconstruction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, double **u_prev, double *xc, double **slopes, int number_of_scalars );

#endif // __MINMOD_H_