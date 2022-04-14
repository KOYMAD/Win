/*
 * godunov.h
 *
 * Метод Годунова расчета потоков для одномерной системы уравнений Баера-Нунзцато.
 *
 * Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of
 * compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526.
 * 
 * (c) Уткин Павел, 2013
 *
 * Создан: 4 июля 2013 г.
 *
 */

#ifndef __GODUNOV_H_
#define __GODUNOV_H_

#include <string.h>
#include "utils.h"
#include "relaxation.h"


// Метод Годунова интегрирования одномерной системы уравнений Баера-Нунциато,
// расчет одного шага по времени в одной ячейке
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_params[M] - вектор примитивных переменных в ячейке слева от рассчитываемой (in)
// center_params[M] - вектор примитивных переменных в рассчитываемой ячейке (in)
// right_params[M] - вектор примитивных переменных в ячейке справа от рассчитываемой (in)
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой (in)
// slopes_center - вектор наклонов в рассчитываемой ячейке (in)
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге (out)
// n - реальный размер векторов
// configuration_pressure - конфигурационное давление
void godunov( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_params[M], double center_params[M],
              double right_params[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
              double dt, double h, double solution[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, double curr_time, double *configuration_pressure );

// Анализ перепада значений объемной доли дисперсной фазы слева и справа от разрыва
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// left_beta - величина объемной доли дисперсной фазы слева от разрыва (in)
// right_beta - величина объемной доли дисперсной фазы справа от разрыва (in)
// Возвращает одну из констант перечисления Disp_phase_cases
Disp_phase_cases what_case( struct ParametersCommon *paramsc, double left_beta, double right_beta );


// Если в ячейке фактически отсутствует дисперсная фаза, то установить в ней фоновые параметры
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке (in/out)
void set_background_state( struct ParametersCommon *paramsc, double solution[M] );

// Расчет полного "потока" через грань ячейки методом Годунова
// params - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_params[M] - вектор примитивных переменных слева от разрыва (in)
// right_params[M] - вектор примитивных переменных справа от разрыва (in)
// solver_part - константа, которая определяет, какую из четырех частей солвера использовать для расчета "потока" (in)
// edge - идентификатор грани, через которую считается "поток" - LEFT или RIGHT (in)
// curr_cell_beta - значение объемной доли в рассчитываемой ячейке (in)
// godunov_flux[M] - полный "поток" Годунова (out)
// n - реальный размер вектора left_params, right_params и godunov_dlux (in)
// Возвращает: SUCCESS          в случае успеха
//             GODUNOV_FAILS    в случае невозможности построить решение задачи Римана
ReturnCodes godunov_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                          Disp_phase_cases solver_part, Direction edge, double curr_cell_beta, double godunov_flux[M], int n ) ;


// Получение решения на следующем шаге по времени для случая отсутствия неконсервативного члена в правой части - 
// полное расщепление фаз
// paramsс - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// center_params[M] - вектор примитивных переменных в рассчитываемой ячейке (in)
// left_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грании ячейки слева от рассчитываемой (in)
// left_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани рассчитываемой ячейки (in)
// right_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грани рассчитываемой ячейки (in)
// right_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани ячейки справа от рассчитываемой (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// solution[M] - вектор примитивных переменных в рассчитываемой ячейке на следующем шаге (out)
// n (in)
void full_decouple_case_sol( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double center_params[M], double left_minus_ncons[M], double left_plus_ncons[M],
                             double right_minus_ncons[M], double right_plus_ncons[M], double dt, double h, double solution[M], int n ) ;


/* Тебуется расщепить систему уравнений и рассчитать сокращенные вектора потоков, но затем
   сформировать полный "поток", используя в качестве объемной доли значение в текущей ячейке.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   left_params_full[M] - полный вектор примитивных переменных слева от разрыва (in)
   right_params_full[M] - полный вектор примитивных переменных справа от разрыва (in)
   beta - значение объемной доли в рассчитываемой ячейке (in)
   
   flux_full[M] - полный "поток" (out) */
void full_decouple_case_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params_full[M], double right_params_full[M],
                              double beta, double flux_full[M] );


/*  Функция расчета консервативной составляющей "потока" методом Годунова
    для одномерной системы уравнений Баера-Нунциато

    paramsc - структура с основными параметрами вычислительного эксперимента (in)
    debug_info - структура с отладочной информацией (in)
    left_params[M] - вектор примитивных переменных слева от разрыва (in)
    right_params[M] - вектор примитивных переменных справа от разрыва (in)
    s - значение автомодельной переменной x/t (in)
    solver_part - константа, которая определяет, какую из четырех частей солвера использовать для расчета "потока" (in)

    flux[M] - консервативная составляющая вектора "потока" (out)
    v_ncons_res[M] - вектор-решение задачи о распаде разрыва, нужен в выходных параметрах для построения точного решения (out)
    v_solid_cont - скорость контактного разрыва в дисперсной фазе (out)
    p_solid_cont_l - давление в дисперсной фазе слева от контактного разрыва в дисперсной фазе (out)
    p_solid_cont_r - давление в дисперсной фазе справа от контактного разрыва в дисперсной фазе (out)
    n - реальный размер вектора left_params, right_params и godunov_dlux (in)
    Возвращает: SUCCEESS        успешная сходимость итераций
                GODUNOV_FAILS   итерации не сошлись, решение задачи Римана построить не удалось */
ReturnCodes godunov_cons_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                               double s, Disp_phase_cases solver_part, double flux[M], double v_ncons_res[M], double *v_solid_cont,
                               double *p_solid_cont_l, double *p_solid_cont_r, int n, double cont_ncons[M] );

// Построение начального приближения - давления и скорости на контактном разрыве в дисперсной фазе для обеих фаз
// для итерационного процесса
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_params[M] - вектор примитивных переменных слева от разрыва (in)
// right_params[M] - вектор примитивных переменных справа от разрыва (in)
// c - вектор со скоростями звука фаз по разные стороны от разрыва (in)
// p_cont_gas - начальное приближение для давлений в газе на контактном разрыве в дисперсной фазе (out)
// v_cont_gas - начальное приближение для скорости газа на контактном разрыве в дисперсной фазе (out)
// p_cont_disp - начальное приближение для давлений в дисперсной фазе на контактном разрыве в дисперсной фазе (out)
// v_cont_disp - начальное приближение для скорости дисперсной фазы на контактном разрыве в дисперсной фазе (out)
// solver_part - константа, которая определяет, какая из конфигураций объемной доли дисперсной фазы реализуется (in)
// Возвращает: SUCCEESS         начальное приближение построено успешно
//             GODUNOV_FAILS    начальное приближение построить не удалось
ReturnCodes calc_p_v_initial_guess( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M],
                                    double right_params[M], double c[K_GENERAL_CASE], Disp_phase_cases solver_part,
                                    double *p_cont_gas, double *v_cont_gas, double *p_cont_disp, double *v_cont_disp );

/* Расчет давлений и скоростей обеих фаз на контактном разрыве в дисперсной фазе
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   sound_velocities[K_GENERAL_CASE] - вектор со скоростями звука фаз по разные стороны от разрыва - 
                                      sound_velocities[GAS_LEFT] - в газе слева, sound_velocities[GAS_RIGHT] - в газе справа,
                                      sound_velocities[DISP_LEFT] - в дисперсной фазе слева, sound_velocities[DISP_RIGHT] - в дисперсной фазе справа (in)
   solver_part - константа, которая определяет, какая из конфигураций объемной доли дисперсной фазы реализуется (in)
   
   solid_discontinuity_pressures[K_GENERAL_CASE] - давление газа и дисперсной фазы по разные стороны от контактного разрыва
                                                   в дисперсной фазе, на входе - начальное приближение, на выходе - результат;
                                                   [GAS_LEFT] - давление слева в газе, [GAS_RIGHT] - давление справа в газе,
                                                   [DISP_LEFT] - давление слева в дисперсной фазе, [DISP_RIGHT] - давление справа в дисперсной фазе (in/out)
   v1 - скорость контактного разрыва в дисперсной фазе (out)
   v2 - скорость контактного разрыва в газовой фазе (out)
   v21 - скорость газа слева от контактного разрыва дисперсной фазы (out)
   v22 - скорость газа справа от контактного разрыва дисперсной фазы (out)
   
   Возвращает: SUCCEESS         успешная сходимость итераций
               GODUNOV_FAILS    итерации не сошлись, решение задачи Римана построить не удалось */ 
ReturnCodes calc_p_v( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                      double sound_velocities[K_GENERAL_CASE], Disp_phase_cases solver_part,
                      double solid_discontinuity_pressures[K_GENERAL_CASE], double *v1, double *v2, double *v21, double *v22 );

/* Рабочая функция заполнения вспомогательного массива со скоростями звука фаз слева и справа от разрыва

   c2l - скорость звука в газовой фазе слева от разрыва (in)
   c2r - скорость звука в газовой фазе справа от разрыва (in)
   c1l - скорость звука в дисперсной фазе слева от разрыва (in)
   c1r - скорость звука в дисперсной фазе справа от разрыва (in)

   sound_velocities[K_GENERAL_CASE] - вектор скоростей звука (out) */
void fill_sound_velocities( double c2l, double c2r, double c1l, double c1r, double sound_velocities[K_GENERAL_CASE] );

/* Рабочая функция инициализации вектора давлений слева и справа от разрыва в дисперсной фазе

   p_cont_gas - давление газа на контактном разрыве без учета разрыва объемной доли дисперсной фазы (in)
   p_cont_solid - давление дисперсной фазы на контактном разрыве без учета разрыва объемной доли дисперсной фазы (in)

   solid_discontinuity_pressures[K_GENERAL_CASE] - вектор начальных давлений слева и справа от разрыва объемной доли
   дисперсной фазы (out) */
void init_solid_discontinuity_pressures( double p_cont_gas, double p_cont_solid,
                                         double solid_discontinuity_pressures[K_GENERAL_CASE] );

/* Итерационная процедура расчета давления и скорости на контактном разрыве в среде без разрыва объемной доли дисперсной фазы

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 155. - Subroutine STARPU.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   v_ncons_l - вектор примитивных переменных слева от разрыва (in)
   v_ncons_r - вектор примитивных переменных справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   phase - идентификатор фазы, для которой ищется начальное приближение - газовая или дисперсная (in)

   p_cont - давление на контактном разрыве (out)
   v_cont - скорость на контактном разрыве (out)
   
   Возвращает: SUCCEESS         успешная сходимость итераций
               GODUNOV_FAILS    итерации не сошлись, решение задачи Римана построить не удалось */
int calc_contact_pressure_velocity( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double *v_ncons_l, double *v_ncons_r,
                                    int vector_size, double cl, double cr, Phase phase, double *p_cont, double *v_cont );

/* Определение начального приближения для расчета давления на контактном разрыве в среде без разрыва объемной доли дисперсной фазы

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l - вектор примитивных переменных слева от разрыва (in)
   v_ncons_r - вектор примитивных переменных справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   phase - идентификатор фазы, для которой ищется начальное приближение - газовая или дисперсная (in)

   Возвращает искомое начальное приближения */
double pressure_initial_guess( struct ParametersCommon *paramsc, double *v_ncons_l, double *v_ncons_r, int vector_size,
                               double cl, double cr, int phase );

/* Расчет функции F, определяющей скорость среды на контактном разрыве в среде без разрыва объемной доли дисперсной фазы, и ее производной по давлению среды DF

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + поправка на двучленное уравнение состояния:
   Годунов С.К. и др. Численное решение многомерных задач газовой динамики. - 
   М.: Наука, 1976. - С. 110 - 111. - Формулы (13.16), (13.17).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   curr_press - давление с предыдущей итерации (in)
   v_ncons - вектор примитивных переменных (in)
   vector_size - размер вектора переменных (in)
   c - скорость звука в невозмущенной среде (in)
   phase - идентификатор фазы, для которой ищется начальное приближение - газовая или дисперсная (in)

   F - значение функции (out)
   DF - значение производной (out) */
void calc_F_and_DF( struct ParametersCommon *paramsc, double curr_press, double *v_ncons, int vector_size, double c, int phase, double *F, double *DF );

/* Расчет функции G, определяющей плотность газа по разные стороны от контактного разрыва в среде без разрыва объемной доли дисперсной фазы,
   и ее производной по давлению DG

   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - Функция G определяется формулами (7).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   curr_press - давление с предыдущей итерации (in)
   v_ncons[M] - вектор примитивных переменных (in)
   
   G - значение функции (out)
   DG - значение производной (out) */
void calc_G_and_DG( struct ParametersCommon *paramsc, double curr_press, double v_ncons[M], double *G, double *DG );

/* Итерационная процедура для поиска давлений в газовой и дисперсной фазах на контактном разрыве в дисперсной фазе
   методом Ньютона-Рафсона
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   debug_info - структура с отладочной информацией (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   dsolver_part - константа, которая определяет, какая из конфигураций объемной доли дисперсной фазы реализуется (in)
   sound_velocities[K_GENERAL_CASE] - вектор со скоростями звука фаз по разные стороны от разрыва - 
                                      sound_velocities[GAS_LEFT] - в газе слева, sound_velocities[GAS_RIGHT] - в газе справа,
                                      sound_velocities[DISP_LEFT] - в дисперсной фазе слева,
                                      sound_velocities[DISP_RIGHT] - в дисперсной фазе справа (in)

   solid_discontinuity_pressures[K_GENERAL_CASE] - давление газа и дисперсной фазы по разные стороны от контактного разрыва
                                                   в дисперсной фазе, на входе - начальное приближение, на выходе - результат;
                                                   [GAS_LEFT] - давление слева в газе, [GAS_RIGHT] - давление справа в газе,
                                                   [DISP_LEFT] - давление слева в дисперсной фазе, [DISP_RIGHT] - давление справа в дисперсной фазе (in/out)
   v1 - скорость контактного разрыва в дисперсной фазе (out)
   v2 - скорость контактного разрыва в газовой фазе (out)
   v21 - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v22 - скорость газа справа от контактного разрыва в дисперсной фазе (out)
   
   Возвращает: SUCCEESS         успешная сходимость итераций
               GODUNOV_FAILS    итерации не сошлись, решение задачи Римана построить не удалось */
int calc_solid_discontinuity_pressures( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M],
                                        double right_params[M], Disp_phase_cases solver_part, double sound_velocities[K_GENERAL_CASE],
                                        double solid_discontinuity_pressures[K_GENERAL_CASE], double *v1, double *v2,
                                        double *v21, double *v22 );

/* Расчет вектор-функции для определения соотношений на контактном разрыве в дисперсной фазе, ее матрицы Якоби, а также скоростей
   контактных разрывов и скоростей газа слева и справа от контактного разрыва в дисперсной фазе.
   Общий случай наличия дисперсной фазы по обе стороны от разрыва.
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - Искомая функция определяется формулами
   (23), (25).

   Общее правило обозначения переменных с двумя индексами:
   - первый индекс - номер фазы: 1 - дисперсная фаза, 2 - газовая фаза
   - второй индекс - положение относительно контактного разрыва: 1 - слева, 2 - справа
   Индексы l и r соответствуют параметрам в невозмущенной среде слева и справа от разрыва, соответственно.
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   c[K_GENERAL_CASE] - вектор скоростей звука в дисперсной и газовой фазах по разные стороны от разрыва (in)
   curr_p[K_GENERAL_CASE] - текущий вектор давлений дисперсной и газовой фаз слева и справа от разрыва пористости (in)

   P[K_GENERAL_CASE] - искомая вектор-функция (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - искомая матрица Якоби (out)
   v_cont_disp - скорость контактного разрыва в дисперсной фазе (out)
   v_cont_gas - скорость контактного разрыва в газовой фазе (out)
   v_gas_left - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v_gas_right - скорость газа справа от контактного разрыва в дисперсной фазе (out) */
void calc_P_and_DP_general_case( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                 double c[M], double curr_p[M], double P[M],
                                 double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                 double *v_gas_left, double *v_gas_right );

/* Функция отбора решения в дисперсной фазе

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
   + поправка на двучленное уравнение состояния: Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann
   problem and a high-resolution Godunov method for a model of compressible two-phase flow // Journal of Computational Physics.
   - 2006. - V. 212. - P. 490 - 526. - Формулы (9) - (11).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор неконсервативных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор неконсервативных переменных справа от разрыва (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   p1 - давление слева от контактного разрыва в дисперсной фазе (in)
   p2 - давление справа от контактного разрыва в дисперсной фазе (in)
   v_cont - скорость на контактном разрыве (in)
   s - значение x/t, для которого отбирается решение (in)

   v_ncons_res[M] - вектор неконсервативных переменных с отобранными компонентами для дисперсной фазы (out) */
void sample_solid_solution( struct ParametersCommon *paramsc, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr,
                            double p1, double p2, double v_cont, double s, double v_ncons_res[M], double cont_ncons[M] );

/* Функция отбора решения в газовой фазе

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор неконсервативных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор неконсервативных переменных справа от разрыва (in)
   cl - скорость звука в газе слева от разрыва (in)
   cr - скорость звука в газе справа от разрыва (in)
   p1 - давление газа слева от контактного разрыва в дисперсной фазе (in)
   p2 - давление газа справа от контактного разрыва в дисперсной фазе (in)
   v_cont_solid - скорость контактного разрыва в дисперсной фазе (in)
   v_cont_gas - скорость контактного разрыва в газовой фазе (in)
   v1 - скорость газа слева от контактного разрыва в дисперсной фазе (in)
   v2 - скорость газа справа от контактного разрыва в дисперсной фазе (in)
   s - значение x/t, для которого отбирается решение (in)

   v_ncons_res[M] - отобранный вектор неконсервативных переменных (out) */
void sample_gas_solution( struct ParametersCommon *paramsc, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr, double p1,
                          double p2, double v_cont_solid, double v_cont_gas, double v1, double v2, double s, double v_ncons_res[M], double cont_ncons[M] );

/* Расчет неконсервативной составляющей вектора "потока"
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model
   of compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - Формула (30).
   
   v_solid_cont - скорость контактного разрыва в дисперсной фазе (in)
   beta_l - объемная доля дисперсной фазы слева от разрыва (in)
   beta_r - объемная доля дисперсной фазы справа от разрыва (in)
   p_solid_l - давление в дисперсной фазе слева от контактного разрыва в дисперсной фазе (in)
   p_solid_r - давление в дисперсной фазе справа от контактного разрыва в дисперсной фазе (in)
   
   v_ncons_term[M] - неконсервативная составляющая вектора "потока" (out) */
void calc_ncons_term( double v_solid_cont, double beta_l, double beta_r, double p_solid_l, double p_solid_r,
                      double v_ncons_term[M] );

/* Расчет вектор-функции для определения соотношений на контактном разрыве в дисперсной фазе, ее матрицы Якоби, а также скоростей
   контактных разрывов и скоростей газа слева и справа от контактного разрыва в дисперсной фазе.
   ОСОБЫЙ СЛУЧАЙ ОТСУТСТВИЯ ДИСПЕРСНОЙ ФАЗЫ СПРАВА ОТ РАЗРЫВА.
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 501 - 502.

   Общее правило обозначения переменных с двумя индексами:
   - первый индекс - номер фазы: 1 - дисперсная фаза, 2 - газовая фаза
   - второй индекс - положение относительно контактного разрыва: 1 - слева, 2 - справа
   Индексы l и r соответствуют параметрам в невозмущенной среде слева и справа от разрыва, соответственно.
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   c[K_GENERAL_CASE] - вектор скоростей звука в дисперсной и газовой фазах по разные стороны от разрыва (in)
   curr_p[K_GENERAL_CASE] - текущий вектор давлений дисперсной и газовой фаз слева и справа от разрыва пористости (in)

   curr_p[GAS_LEFT] - давление газа слева от контактного разрыва в дисперсной фазе
   curr_p[GAS_RIGHT] - давление газа справа от контактного разрыва в дисперсной фазе
   curr_p[DISP_LEFT] - давление дисперсной фазы слева от контактного разрыва в дисперсной фазе

   P[K_GENERAL_CASE] - искомая вектор-функция (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - искомая матрица Якоби (out)
   v_cont_disp - скорость контактного разрыва в дисперсной фазе (out)
   v_cont_gas - скорость контактного разрыва в газовой фазе (out)
   v_gas_left - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v_gas_right - скорость газа справа от контактного разрыва в дисперсной фазе (out) */
void calc_P_and_DP_left_disp_phase( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                    double c[K_GENERAL_CASE], double curr_p[K_GENERAL_CASE], double P[M],
                                    double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                    double *v_gas_left, double *v_gas_right );

/* Расчет вектор-функции для определения соотношений на контактном разрыве в дисперсной фазе, ее матрицы Якоби, а также скоростей
   контактных разрывов и скоростей газа слева и справа от контактного разрыва в дисперсной фазе.
   ОСОБЫЙ СЛУЧАЙ ОТСУТСТВИЯ ДИСПЕРСНОЙ ФАЗЫ СЛЕВА ОТ РАЗРЫВА.
   
   Формулы приведены в диссертации.

   Общее правило обозначения переменных с двумя индексами:
   - первый индекс - номер фазы: 1 - дисперсная фаза, 2 - газовая фаза
   - второй индекс - положение относительно контактного разрыва: 1 - слева, 2 - справа
   Индексы l и r соответствуют параметрам в невозмущенной среде слева и справа от разрыва, соответственно.
   
   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   left_params[M] - вектор примитивных переменных слева от разрыва (in)
   right_params[M] - вектор примитивных переменных справа от разрыва (in)
   c[K_GENERAL_CASE] - вектор скоростей звука в дисперсной и газовой фазах по разные стороны от разрыва (in)
   curr_p[K_GENERAL_CASE] - текущий вектор давлений дисперсной и газовой фаз слева и справа от разрыва пористости (in)
   
   curr_p[GAS_LEFT] - давление газа слева от контактного разрыва в дисперсной фазе
   curr_p[GAS_RIGHT] - давление газа справа от контактного разрыва в дисперсной фазе
   curr_p[DISP_RIGHT] - давление дисперсной фазы справа от контактного разрыва в дисперсной фазе

   P[K_GENERAL_CASE] - искомая вектор-функция (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - искомая матрица Якоби - последние строка и столбец не имеют смысла (out)
   v_cont_disp - скорость контактного разрыва в дисперсной фазе (out)
   v_cont_gas - скорость контактного разрыва в газовой фазе (out)
   v_gas_left - скорость газа слева от контактного разрыва в дисперсной фазе (out)
   v_gas_right - скорость газа справа от контактного разрыва в дисперсной фазе (out) */
void calc_P_and_DP_right_disp_phase( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                     double c[M], double curr_p[M], double P[M],
                                     double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                     double *v_gas_left, double *v_gas_right );

// Метод Годунова интегрирования одномерной системы уравнений Баера-Нунциато при отсутствия градиента объемной доли
// дисперсной фазы, расчет одного шага по времени в одной ячейке для одной фазы
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// center_params_full[M] - вектор примитивных переменных в рассчитываемой ячейке для полной системы (in)
// left_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грании ячейки слева от рассчитываемой для полной системы (in)
// left_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани рассчитываемой ячейки для полной системы (in)
// right_minus_ncons[M] - реконструированный вектор примитивных переменных на правой грани рассчитываемой ячейки для полной системы (in)
// right_plus_ncons[M] - реконструированный вектор примитивных переменных на левой грани ячейки справа от рассчитываемой для полной системы (in)
// phase - идентификатор фазы, для которой интегрируется система уравнений - газовая или дисперсная (in)
// dt - временной шаг (in)
// h - пространственный шаг (in)
// v_ncons_res (out)
// solution_full[M] - вектор примитивных переменных в рассчитываемой ячейке для полной системы на следующем шаге (out)
// n - реальный размер векторов center_params_full, left_minus_ncons_full, left_plus_ncons_full, right_minus_ncons_full, right_plus_ncons_full, solution_full
void godunov_classical_one_phase( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double center_params_full[M],
                                  double left_minus_ncons_full[M], double left_plus_ncons_full[M], double right_minus_ncons_full[M],
                                  double right_plus_ncons_full[M], Phase phase, double dt, double h, double v_ncons_res_left[M], 
                                  double v_ncons_res_right[M], double solution_full[M], int n );

// Расчет классического однофазного потока С.К. Годунова
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных слева от разрыва (in)
// right_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных справа от разрыва (in)
// phase - идентификатор фазы, для которой рассчитывается поток - газовая или дисперсная (in)
// v_ncons_res[M_REDUCTION] - рассчитанный вектор неконсервативных переменных (in)
// flux[M_REDUCTION] - рассчитанный вектор потока (out)
void godunov_flux_classical( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                             double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons_res[M_REDUCTION], double flux[M_REDUCTION], double cont_ncons_red[M_REDUCTION] );

// Построение решение классической однофазной задачи Римана
// params - структура с параметрами вычислительного эксперимента (in)
// debug_info - структура с отладочной информацией (in)
// left_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных слева от разрыва (in)
// right_ncons_params[M_REDUCTION] - однофазный  вектор примитивных переменных справа от разрыва (in)
// phase - идентификатор фазы, для которой рассчитывается поток - газовая или дисперсная (in)
// v_ncons[M_REDUCTION] - рассчитанный вектор-решение (out)
void get_classical_Riemann_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                                     double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] );

/* Функция отбора решения в одной из фаз по сокращенному однофазному вектору

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
   + поправка на двучленное уравнение состояния: Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann
   problem and a high-resolution Godunov method for a model of compressible two-phase flow // Journal of Computational Physics.
   - 2006. - V. 212. - P. 490 - 526. - Формулы (9) - (11).

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons_l - вектор неконсервативных переменных слева от разрыва (in)
   v_ncons_r - вектор неконсервативных переменных справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   phase - фаза, для которой отбирается решение - газовая или дисперсная (in)
   p_cont - давление на контактном разрыве (in)
   v_cont - скорость на контактном разрыве (in)
   s - значение x/t, для которого отбирается решение (in)

   v_ncons_res - вектор неконсервативных переменных с отобранными компонентами (out) */
void sample_reduced( struct ParametersCommon *paramsc, double *v_ncons_l, double *v_ncons_r, int vector_size, double cl, double cr,
                     Phase phase, double p_cont, double v_cont, double s, double *v_ncons_res, double cont_ncons_red[M_REDUCTION]  );

/* Функция для построения массива точек, соответствующих множеству состояний, куда можно перевести однофазную среду
   слева или справа от разрыва, пустив по ней ударную волну или волну разрежения

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons - вектор неконсервативных переменных слева или справа от разрыва (in)
   vector_size - размер вектора переменных (in)
   c - скорость звука слева или справа от разрыва (in)
   phase - фаза, для которой строится адиабата - газовая или дисперсная (in)
   dir - направление, для которого строится адиабата - слева или справа от разрыва (in) */
void draw_adiabatic_curve( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *v_ncons, int vector_size, double c, Phase phase, Direction dir );

// Возвращает значение объемной доли твердой фазы, используемое в full_decouple_case_flux (случай NO_GRAD) при расчете потока
// case_beta - номер рассматриваемого варианта ( 0 - берется значение в центре ячейки, 1 - берется значение слева 
//             или справа от ребра, через которое считается поток, в зависимости от знака скорости дисперсной фазы) (in)
// left_edge_beta - значение объемной доли дисперсной фазы слева от ребра, через которое считается поток (in)
// center_beta - значение объемной доли дисперсной фазы в центре рассчитываемой ячейки (in)
// right_edge_beta - значение объемной доли дисперсной фазы справа от ребра, через которое считается поток (in)
// solid_velocity - значение скорости дисперсной фазы, получаемое в результате решения задачи Римана о распаде разрыва (in)
double full_decouple_case_flux_volume_fraction ( int case_beta, double left_edge_beta, double center_beta, double right_edge_beta, double solid_velocity );

#endif /* __GODUNOV_H_ */