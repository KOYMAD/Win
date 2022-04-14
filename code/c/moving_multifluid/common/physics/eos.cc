// eos.cc
// Уравнения состояния фаз
// (c) Уткин Павел, 2015
// Создан: 1 октября 2015 г.

#include "eos.h"

// Возвращает внутреннюю энергию газовой фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double e_gas( const struct ParametersCommon* paramsc, const double p, const double r ) {

    return ( ( p + paramsc->g2 * paramsc->p02 )  / r / ( paramsc->g2 - 1.0 )  / (1.0 + paramsc->b_virial * r) );

}

// Возвращает давление газовой фазы как функцию внутренней энергии и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// e - внутренняя энергия
// r - плотность
double p_gas( const struct ParametersCommon* paramsc, const double e, const double r ) {

    double pressure_gas = ( e * r * ( paramsc->g2 - 1.0 ) * (1.0 + paramsc->b_virial * r) - paramsc->g2 * paramsc->p02);
    if (pressure_gas < 0.0){
   //     printf( "\np_gas -> Pressure of gas phase should not be negative.\n\n" );
   //     exit( EXIT_FAILURE );
    }
    return pressure_gas;

}

// Возвращает внутреннюю энергию дисперсной фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double e_disp( const struct ParametersCommon* paramsc, const double p, const double r ) {

    return ( ( p + paramsc->g1 * paramsc->p01 )  / r / ( paramsc->g1 - 1.0 ) );

}

// Возвращает давление дисперсной фазы как функцию внутренней энергии и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// e - внутренняя энергия
// r - плотность
double p_disp( const struct ParametersCommon* paramsc, const double e, const double r ) {
    
    double pressure_disp = ( e * r * ( paramsc->g1 - 1.0 ) - paramsc->g1 * paramsc->p01 );

    if (pressure_disp < 0.0){
      //  printf( "\np_disp -> Pressure of solid phase should not be negative.\n\n" );
      //  exit( EXIT_FAILURE );
    }
    
    return pressure_disp;

}

// Возвращает температуру дисперсной фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double T_disp( const struct ParametersCommon* paramsc, const double p, const double r ) {

    return ( ( p + paramsc->g1*paramsc->p01 )  / r / ( paramsc->g1 - 1.0 ) / paramsc->specific_heat_cv1);

}

// Возвращает температуру газовой фазы как функцию давления и плотности
// paramsc - структура с основными параметрами вычислительного эксперимента
// p - давление
// r - плотность
double T_gas( const struct ParametersCommon* paramsc, const double p, const double r ) {

    return ( ( p + paramsc->p02 )  / r / ( paramsc->g2 - 1.0 )  / (1.0  + paramsc->b_virial * r) / paramsc->specific_heat_cv2);

}