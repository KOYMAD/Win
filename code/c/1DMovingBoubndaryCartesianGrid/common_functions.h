#ifndef __COMMON_FUNCTIONS_H
#define __COMMON_FUNCTIONS_H

#include "structures.h"
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

void calc_conservative_variables ( struct Parameters *params, struct Primitive_vector primitive, struct Conservative_vector *conservative );

void calc_primitive_variables ( struct Parameters *params, struct Conservative_vector conservative, struct Primitive_vector *primitive );

void calc_sound_velocity ( const struct Parameters* params, struct Primitive_vector primitive, double* c );

double max( double a, double b );

double min( double a, double b );

#endif /* __COMMON_FUNCTIONS_H */