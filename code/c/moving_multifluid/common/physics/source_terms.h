#ifndef __SOURCE_TERMS_H_
#define __SOURCE_TERMS_H_

#include "eos.h"

// Возвращает в виде вектора начальные параметры в конкретном блоке
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами для 1d2phc случая (in)
// params2d - структура с дополнительными параметрами для 2d2phc случая (in)
// number_of_block[2] - массив с номером блока. В 1d случае имеет значение только number_of_block[0] (in)
// v_ncons[M] - вектор наачльных примитивных переменных в блоке (out)
void calc_initial_parameters_values(const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, int *number_of_block, double v_ncons[M]);

// Расчет силы межфазного трения
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает силу трения
double calc_friction_force( struct ParametersCommon *paramsc, double v_ncons[M] );

// Расчет коэффициента сопротивления
// Houim R.W., Oran E.S. A multiphase model for compressible granular-gaseous flows: formulation and initial tests // J. Fluid Mech. - 2016. - V. 789. - P. 166 - 220.
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает коэффициент сопротивления
double calc_Cd( struct ParametersCommon *paramsc, double v_ncons[M] );

// Расчет скорости компактирования
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с дополнительными параметрами для 1d2phc случая (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает скорость компактирования
double calc_compaction_rate( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block );

// Расчет configuration pressure Beta
// params1d - структура с дополнительными параметрами для 1d2phc случая (in)
// v_ncons[M] - текущий вектор примитивных переменных в ячейке (in)
// Возвращает configuration pressure Beta
double calc_configuration_pressure ( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block );

double compaction_energy( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block );

// Расчет химической реакции
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных в центре ячейки (in)
double calc_chemical_reaction(const struct ParametersCommon *paramsc, double v_ncons[M], int ignition_flag);

// Расчет теплопереноса
// paramsc - структура с общими параметрами вычислительного эксперимента (in)
// v_ncons[M] - вектор примитивных переменных в центре ячейки (in)
double calc_heat_transfer(const struct ParametersCommon *paramsc, double v_ncons[M]);

#endif // __SOURCE_TERMS_H_