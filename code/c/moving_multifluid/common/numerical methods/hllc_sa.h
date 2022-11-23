// hllc_sa.h
// Метод HLLC численного интегрирования уравнений Saurel-Abgrall
// Реализовано по: Li Q. et al. Difference scheme for two-phase flow // Applied Mathematics and Mechanics. - 2004. - V. 25, No. 5. - P. 536 - 545
// и
// Liang S. et al. Solving seven-equation model for compressible two-phase flow using multiple GPUs // Computers & Fluids. - 2014. - V. 99. - P. 156 - 171
// (c) Уткин Павел, 2017
// Создан: 29 марта 2017 г.

#ifndef __HLLC_H_
#define __HLLC_H_

#include "relaxation.h"
#include "hll.h"

// Метод HLLC численного интегрирования уравнений Saurel-Abgrall
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
// n - реальный размер векторов
void hllc_1d(struct ParametersCommon* paramsc, struct Parameters1d* params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
          double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure,  double body_velocity, int i, int *status,double cont_left[M], double cont_right[M] ) ;

// Гиперболический оператор метода HLLC
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
void Lh_HLLC( struct DebugInfo *debug_info, struct Parameters1d* params1d,  struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars,  int i, int *status, double body_velocity,double cont_left[M], double cont_right[M]  );

// Метод Harten - Lax - van Leer - Contact (HLLC) расчета потоков в двухфазной среде
// Liang S. et al. Solving seven-equation model for compressible two-phase flow using multiple GPUs //
// Computers & Fluids. - 2014. - V. 99. - P. 156 - 171.
// paramsc - структура с основными параметрами вычислительного эксперимента
// left_ncons - вектор примитивных переменных слева от разрыва
// right_ncons - вектор примитивных переменных справа от разрыва
// flux - рассчитанный вектор потока
// phi - объемная доля дисперсной фазы на ребре
// s_cont - оценка для скорости контактного разрыва
void hllc_flux( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M], array1D* flux, double* phi, double* s_cont );

// Расчет вектора - оценки состояния слева или справа от контактного разрыва в методе HLLC
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных слева или справа от разрыва
// s - оценка скорости левой или правой волны
// s_star - оценка для s со звездочкой
// q_star - оценки состояния слева или справа от контактного разрыва
void calc_contact_vector( const struct ParametersCommon *paramsc, const double v_ncons[M], const double s, const double s_star, array1D* q_star );

#endif // __HLLC_H_