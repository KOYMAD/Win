// relaxation.h
// Релаксация скоростей и давлений фаз на межфазной границе

// Релаксация скоростей реализована по статье:
// Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467.

// Реализация давлений реализована по статье:
// Иванов И.Э. Численное моделирование многофазных течений с большим содержанием дисперсной фазы // Вестник МАИ. – 2009. – Т. 16, № 2. – С. 62 – 70.
// Алгоритм детально описан в \science\utkin\docs\Релаксация давления.docx

// (c) Уткин Павел, 2017
// Создан: 29 марта 2017 г.

#ifndef __RELAXATION_H_
#define __RELAXATION_H_

#include "struct.h"
#include "utils.h"
#include "eos.h"
#include "io.h"
#include "source_terms.h"

// Релаксационный оператор для решения системы уравнений типа Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// left_ncons - вектор примитивных переменных в ячейке слева от рассчитываемой
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
// right_ncons - вектор примитивных переменных в ячейке справа от рассчитываемой
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия релаксационного оператора
// n - реальный размер вектора
// curr_time - текущий момент времени
// configuration_pressure - конфигурационное давление
void Lr( const struct ParametersCommon* paramsс, const struct Parameters1d* params1d, struct DebugInfo *debug_info, const double left_ncons[M], const double center_ncons[M],
         const double right_ncons[M], const double dt, const double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure);

// Оператор релаксации скорости для решения системы уравнений типа Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора
// dt - временной шаг
// h - пространственный шаг
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия релаксационного оператора скорости
void Lrv( const struct ParametersCommon* paramsc, const double center_ncons[M], const double dt, const double h, double solution_ncons[M], int n );

// Оператор релаксации давления для решения системы уравнений типа Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора и оператора релаксации скорости
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического и полного релаксационного операторов
void Lrp_Ivanov( const struct ParametersCommon* paramsc, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int n );

// Оператор релаксации давления для решения системы уравнений типа Baer-Nunziato с учетом компактирования Schwendeman
// paramsc - структура с основными параметрами вычислительного эксперимента
// params1d - структура c параметрами вычислительного эксперимента, присущими 1d случаю
// debug_info - структура с отладочной информацией
// center_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического оператора и оператора релаксации скорости
// solution_ncons - вектор примитивных переменных в рассчитываемой ячейке после действия гиперболического и полного релаксационного операторов
// step_number - номер текущей ячейки
void Lrp_compaction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int step_number, int n, double *configuration_pressure);

#endif // __RELAXATION_H_