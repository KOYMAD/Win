// utils.h
// Рабочие функции, специфичные для рассматриваемого класса систем уравнений типа Baer-Nunziato и Saurel-Abgrall
// вне зависимости от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 26 августа 2018 г.

#ifndef __UTILS_H_
#define __UTILS_H_

#include "struct.h"
#include "utils_bn.h"
#include "io_1d.h"

// различные разностные схемы для решения уравнений БН
#include "cir_1_bn.h"
#include "cir_2_bn.h"
#include "cir_3_bn.h"
#include "cir_4_bn.h"
#include "godunov_bn.h"

// методы решения уравнений СА
#include "hll.h"
#include "rusanov_sa_1d2phc.h"
#include "hllc_sa.h"

void initiate_status( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double left, double right, int *status, int *left_ghost, int *right_ghost );

void calc_status (  struct ParametersCommon *paramsc, struct Parameters1d *params1d, double coordinate_of_left_boundary_of_body, double coordinate_of_right_boundary_of_body,
    int *status, int *index_of_left_ghost_cell, int *index_of_right_ghost_cell, double **v_ncons, double v_left[M], double v_right[M]);

void convert_2D_to_1D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector);

void convert_1D_to_2D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector);

// Инициализация вектора-решения
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// files_directory - директория, в которой находятся все внешние файлы, требуемые для расчета (in)
// time_mom - структура, определяющая текущий момент времени (out)
// **initial_solution - массив векторов в примитивных переменных - начальных условий в каждой ячейке (out)
void init_solution_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, char *files_directory, struct TimeMoment *time_mom, double **initial_solution );

// Обработка граничного условия
// params1d - структура с параметрами одномерной задачи (in)
// v_ncons[M] - вектор примитивных переменных (in)
// boun_type - тип граничного условия (in)
// boun_v[M] - вектор примитивных переменных в фиктивной ячейке (out)
// curr_t - текущий момент времени (in)
// curr_inflow_parameters[M] - текущие значения параметров втекания(in)
void boundary( struct Parameters1d *params1d, struct Parameters2d *params2d, double v_ncons[M], int boun_type, double boun_v[M], double curr_t, double curr_inflow_parameters[M] ) ;

void current_inflow_parameters( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, int boun_direction, double *curr_inflow_params);

// Расчет одного шага по времени в одной ячейке
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// u_left - вектор примитивных переменных в ячейке слева от рассчитываемой
// u_center - вектор примитивных переменных в рассчитываемой ячейке
// u_right - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - шаг интегрирования по времени
// h - размер ячейки
// u_next - вектор примитивных переменных в рассчитываемой ячейке на следующем временном шаге
// n - реальный размер векторов
// is_pressure_relaxation_now - true, если после решения одномерной задачи по данному направлению нужно применять релаксацию давлений; false - иначе
//                              если в данной задачу вообще не нужно применять релаксацию давлений, то true не повлияет на результат, так как данная опция контролируется в parameters.dat
// number_of_scalars - число лагранжевых скаляров в данной задаче
// curr_time - текущий момент времени
// configuration_pressure - конфигурационное давление
void calc_step_in_cell( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double u_left[M], double u_center[M],
                        double u_right[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                        double dt, double h, double u_next[M], int step_number, int n, bool is_pressure_relaxation_now, int number_of_scalars, double curr_time, double *configuration_pressure, double body_velocity, int i, int *status, double cont_left[M], double cont_right[M] ) ;

// Инициализация вектора-решения
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// files_directory - директория, в которой находятся все внешние файлы, требуемые для расчета (in)
// time_mom - структура, определяющая текущий момент времени (out)
// **initial_solution - массив векторов в примитивных переменных - начальных условий в каждой ячейке (out)
void init_solution_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, char *files_directory, struct TimeMoment *time_mom, double **initial_solution );

// Инициализация структуры для хранения и передачи отладочной информации
// output_file_directory - директория, в которую должны быть записаны файлы
// debug_info - структура с отладочной информацией
void init_debug_info( const char *output_file_directory, struct DebugInfo *debug_info );

// Расчет усредненного межфазного давления
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор неконсервативных переменных в рассматриваемой ячейке
// Возвращает искомое давление
double calc_p_i( const struct ParametersCommon* paramsc, const double v_ncons[M] );

// Расчет усредненной межфазной плотности, имеет смысл только для модели Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор неконсервативных переменных в рассматриваемой ячейке
// Возвращает искомую плотность
double calc_r_i( const struct ParametersCommon* paramsc, const double v_ncons[M] );

// Расчет скоростей звука для обеих фаз по полному вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// c1 - скорость звука в дисперсной фазе
// c2 - скорость звука в газовой фазе
void calc_sound_velocity( const struct ParametersCommon* paramsc, const double v_ncons[M], double* c1, double* c2 );

// Расчет скоростей звука для заданной фазы по вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// phase - идентификатор фазы, для которой рассчитывается скорость звука - газовая или дисперсная
// Возвращает скорость звука
double calc_sound_velocity_one_phase( const struct ParametersCommon* paramsc, const double v_ncons[M], const Phase phase );

// Расчет скорости звука для одной фазы по однофазному вектору переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// phase - идентификатор фазы, для которой рассчитывается скорость звука - газовая или дисперсная
// Возвращает скорость звука
double calc_sound_velocity_reduced( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase );

// Расчет усредненной межфазной скорости
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор неконсервативных переменных в рассматриваемой ячейке
// Возвращает искомую скорость
double calc_u_i( const struct ParametersCommon* paramsc, const double v_ncons[M] );

// Преобразование вектора "консервативных" переменных в вектор примитивных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// v_ncons - вектор примитивных переменных
// number_of_scalars - количество лагранжевых скаляров
// Возвращает: SUCCEESS            давление в обеих фазах положительно
//             NEGATIVE_PRESSURE   давление хотя бы в одной из фаз отрицательно
ReturnCodes convert_cons_to_noncons( const struct ParametersCommon* paramsc, const double v_cons[M], double v_ncons[M], int number_of_scalars  );

/* Преобразование редуцированного однофазного вектора "консервативных" переменных в вектор примитивных

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_cons[M_REDUCTION] - вектор "консервативных" переменных (in)
   phase - идентификатор фазы, для которой осуществляется преобразование - газовая или дисперсная (in)

   v_ncons[M_REDUCTION] - вектор примитивных переменных (out) */
void convert_cons_to_noncons_reduction( struct ParametersCommon *paramsc, double v_cons[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] );

/* Формирование однофазного вектора по полному двухфазному вектору
 
   full_vector[M] - полный двухфазный вектор (in)
   phase - идентификатор фазы, для которой по полному вектору строится редуцированный - газовая или дисперсная (in)
 
   reduced_vector[M_REDUCTION] - редуцированный однофазный вектор (out) */
void convert_full_to_reduced( double full_vector[M], Phase phase, double reduced_vector[M_REDUCTION] );

// Преобразование вектора примитивных переменных в вектор "консервативных"
// params - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// v_cons - вектор "консервативных" переменных
// number_of_scalars - количество лагранжевых скаляров
void convert_noncons_to_cons( const struct ParametersCommon *paramsc, const double v_ncons[M], double v_cons[M], int number_of_scalars );

/* Преобразование однофазного вектора примитивных переменных в вектор "консервативных"

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons[M_REDUCTION] - вектор примитивных переменных (in)
   phase - идентификатор фазы, для которой осуществляется преобразование - газовая или дисперсная (in)

   v_cons[M_REDUCTION] - вектор "консервативных" переменных (out) */
void convert_noncons_to_cons_reduction( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase, double v_cons[M_REDUCTION] );

// Формирование части полного двухфазного вектора по сокращенному однофазному
// reduced_vector[M_REDUCTION] - редуцированный однофазный вектор (in)
// phase - идентификатор фазы, для которой по редуцированному вектору строится часть полного - газовая или дисперсная (in)
// full_vector[M] - полный двухфазный вектор (out)
void convert_reduced_to_full( double reduced_vector[M], Phase phase, double full_vector[M_REDUCTION] );

// Расчет вектора дифференциального "потока" по вектору "консервативных" переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// flux - рассчитываемый вектор дифференциального "потока"
void diff_flux_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* flux, int number_of_scalars ) ;

// Расчет вектора дифференциального "потока" по вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// flux - рассчитываемый вектор дифференциального "потока"
void diff_flux_ncons( const struct ParametersCommon* paramsc, const double v_ncons[M], array1D* flux, int number_of_scalars ) ;

/* Расчет однофазного вектора дифференциального "потока" по сокращенному однофазному вектору примитивных переменных

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons - сокращенный однофазный вектор примитивных переменных (in)
   phase - фаза - газовая или дисперсная, для которой рассчитывается "поток" (in)

   flux - рассчитываемый сокращенный однофазный вектор дифференциального "потока" (out) */
void diff_flux_ncons_reduced( struct ParametersCommon *paramsc, double *v_ncons, Phase phase, double *flux );

// Расчет неконсервативного вектора правых частей по вектору "консервативных" переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// rhst - неконсервативный вектор правых частей
void rhst_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* rhst, int number_of_scalars );

// Расчет неконсервативного вектора правых частей по вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// rhst - неконсервативный вектор правых частей
void rhst_ncons( const struct ParametersCommon* paramsс, const double v_ncons[M], array1D* rhst, int number_of_scalars );

// Отладочная печать компонент вектора в ячейке
// debug_info - структура с отладочной информацией
void vector_debug_print( const struct DebugInfo *debug_info );

// Печать полной информации о проблемном месте в случае аварийной остановки
// debug_info - структура с отладочной информацией (in)
// index - индекс проблемного элемента в массиве давлений (in)
void debug_print( struct DebugInfo *debug_info, int index );

// Приведение параметров задачи к безразмерному виду
// paramsc - структура с основными параметрами вычислительного эксперимента, часть полей проходит процедуру приведения к безразмерному виду (in/out)
// params1d - структура с параметрами одномерной задачи, часть полей проходит процедуру приведения к безразмерному виду (in/out)
void dimensionalization( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d ) ;

void current_block_number(struct ParametersCommon *params, struct Parameters1d *params1d, struct Parameters2d *params2d, int step_number_x, int step_number_y, int *number_of_block);

#endif // __UTILS_H_