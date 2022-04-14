// hll.h
// Метод HLL из Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) Уткин Павел, 2015
// Создан: 8 июля 2015 г.

#ifndef __HLL_H_
#define __HLL_H_

#include "struct.h"
#include "utils.h"
#include "eos.h"
#include "relaxation.h"

// Метод HLL численного интегрирования уравнений Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией (in)
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге
// step_number - номер текущего шага по пространству
// n - реальный размер векторов
// is_pressure_relaxation_after_this_step - true, если после решения одномерной задачи по данному направлению нужно применять релаксацию давлений; false - иначе
//                              если в данной задачу вообще не нужно применять релаксацию давлений, то true не повлияет на результат, так как данная опция контролируется в parameters.dat
// number_of_scalars - число лагранжевых скаляров
// curr_time - текущий момент времени
// configuration_presure - конфигурационное давление
void hll( struct ParametersCommon* paramsc, struct Parameters1d* params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
          double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure, double S_center, double S_left, double S_right, int l );

// Гиперболический оператор метода HLL
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
// n - реальный размер векторов без учета лагранжевых скаляров
// number_of_scalars - количество дополнительных уравнений, оно же количество лагранжевых скаляров
void Lh_HLL( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars, double S_center,  double S_left, double S_right   );

// Метод Harten - Lax - van Leer (HLL) расчета потоков в двухфазной среде
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных слева от разрыва
// right_ncons - вектор примитивных переменных справа от разрыва
// flux - рассчитанный вектор потока
// splus - оценка для скорости волны S+
// sminus - оценка для скорости волны S-
// n - реальный размер векторов
void hll_flux( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M],
               array1D* flux, double* splus, double* sminus, int n, int number_of_scalars  ) ;

// Расчет оценок для скоростей волн S+ и S-
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons[M] - вектор примитивных переменных слева от разрыва
// right_ncons[M] - вектор примитивных переменных справа от разрыва
// phase - фаза, для которой оцениваются скорости
// splus - оценка для скорости волны S+
// sminus - оценка для скорости волны S-
void calc_splus_sminus( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M],
                        const Phase phase, double* splus, double* sminus );

// Аппроксимация дифференциала dB2 в правой части
// left_minus - левая грань, вектор консервативных переменных слева от разрыва
// left_plus - левая грань, вектор консервативных переменных справа от разрыва
// right_minus - правая грань, вектор консервативных переменных слева от разрыва
// right_plus - правая грань, вектор консервативных переменных справа от разрыва
// splus_l - оценка для скорости волны S+ на левом ребре
// sminus_l - оценка для скорости волны S- на левом ребре
// splus_r - оценка для скорости волны S+ на правом ребре
// sminus_l - оценка для скорости волны S- на правом ребре
// dB2dx - искомая аппроксимация
void calc_dB2( const double left_minus[M], const double left_plus[M], const double right_minus[M], const double right_plus[M],
               double splus_l, double sminus_l, double splus_r, double sminus_r, double* dB2dx, double S_center, double S_left, double S_right );

/*
// Релаксационный оператор метода HLL
// params - структура с параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия релаксационного оператора
void Lr( const struct Parameters1d* params, const double left_ncons[M], const double center_ncons[M],
         const double right_ncons[M], const double dt, const double h, double solution_ncons[M] );

// Оператор релаксации скорости метода HLL
// params - структура с параметрами вычислительного эксперимента
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия источникого релаксационного оператора скорости
void Lrv( const struct Parameters1d* params, const double center_ncons[M], const double dt, const double h, double solution_ncons[M] );

// Оператор релаксации давления метода HLL
// ---
// Реализация по статье:
// Иванов И.Э. Численное моделирование многофазных течений с большим содержанием дисперсной фазы // Вестник МАИ. – 2009. – Т. 16, № 2. – С. 62 – 70.
// Алгоритм детально описан в \science\utkin\docs\Релаксация давления.docx
// ---
// params - структура с параметрами вычислительного эксперимента
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора и оператора релаксации скорости
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического и полного релаксационного операторов
void Lrp_Ivanov( const struct Parameters1d* params, const double center_ncons[M], double solution_ncons[M] );
*/
#endif // __HLL_H