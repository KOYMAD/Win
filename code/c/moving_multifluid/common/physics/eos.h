// eos.h
// Уравнения состояния фаз
// (c) Уткин Павел, 2015
// Создан: 1 октября 2015 г.

#ifndef __EOS_H
#define __EOS_H

#include "struct.h"

// Возвращает внутреннюю энергию газовой фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double e_gas( const struct ParametersCommon* paramsc, const double p, const double r );

// Возвращает давление газовой фазы как функцию внутренней энергии и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// e - внутренняя энергия
// r - плотность
double p_gas( const struct ParametersCommon* paramsc, const double e, const double r );

// Возвращает внутреннюю энергию дисперсной фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double e_disp( const struct ParametersCommon* paramsc, const double p, const double r );

// Возвращает давление дисперсной фазы как функцию внутренней энергии и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// e - внутренняя энергия
// r - плотность
double p_disp( const struct ParametersCommon* paramsc, const double e, const double r );

// Возвращает температуру дисперсной фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double T_disp( const struct ParametersCommon* paramsc, const double p, const double r );

// Возвращает температуру газовой фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double T_gas( const struct ParametersCommon* paramsc, const double p, const double r );

#endif // __EOS_H