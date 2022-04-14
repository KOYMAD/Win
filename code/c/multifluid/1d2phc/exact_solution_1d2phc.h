/*
 * exact_solution_1d2phc.h
 *
 * Построение точного решения задачи о распаде разрыва для системы уравнений Баера - Нунзиато.
 *
 * (c) Уткин Павел, 2013
 *
 * Создан: 24 августа 2013 г.
 *
 */

#ifndef __EXACT_SOLUTION_1D2PHC_H
#define __EXACT_SOLUTION_1D2PHC_H

#include "godunov_bn.h"
#include "utils.h"
#include "memory.h"
#include "grid.h"

// точное решение изначально строится на отрезке [-0.5;0.5], потом сдвигается, иходя из геометрии расчетной области
#define LEFT_BOUN_EX_SOL    -0.5
#define RIGHT_BOUN_EX_SOL   0.5

// Построение точного решения задачи о распаде разрыва для системы уравнений Баера-Нунзиато
// output_directory - директория, куда будет записан файл с точным решением (in) 
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с параметрами одномерной задачи (in)
// debug_info - структура с отладочной информацией (in)
void build_exact_sol( char *output_directory, struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, int n );

#endif /* __EXACT_SOLUTION_1D2PHC_H */