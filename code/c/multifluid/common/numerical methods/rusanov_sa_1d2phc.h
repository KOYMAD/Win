// rusanov_sa_1d2phc.h
// Метод Русанова численного интегрирования уравнений Saurel-Abgrall, 1D случай
// Реализовано по: Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) Уткин Павел, 2018
// Создан: 4 марта 2018 г.

#ifndef __RUSANOV_SA_1D2PHC_H_
#define __RUSANOV_SA_1D2PHC_H_

#include "relaxation.h"

// Метод Русанова численного интегрирования уравнений Saurel-Abgrall, 1D случай
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге
// configuration_pressure - конфигурационное давление
void rusanov_1d( struct ParametersCommon* paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
                 double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                 double dt, double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure ) ;

// Гиперболический оператор метода Русанова, 1D случай
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
void Lh_rusanov_1d( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
                    const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                    const double dt, const double h, double solution_ncons[M] );

// Метод Русанова расчета потоков в двухфазной среде, 1D случай
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных слева от разрыва
// right_ncons - вектор примитивных переменных справа от разрыва
// flux - рассчитанный вектор потока
// smax - оценка для скорости S+
void rusanov_flux_1d( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M], array1D* flux, double *splus );

// Расчет оценки для скорости волны S+ с учетом обеих фаз
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons[M] - вектор примитивных переменных слева от разрыва
// right_ncons[M] - вектор примитивных переменных справа от разрыва
// splus - оценка для скорости волны S+
void calc_splus_rusanov( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M], double* splus );

#endif // __RUSANOV_SA_1D2PHC_H_