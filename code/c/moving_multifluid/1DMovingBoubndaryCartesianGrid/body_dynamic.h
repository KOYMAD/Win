#ifndef __BODY_DYNAMIC_H
#define __BODY_DYNAMIC_H

#include "structures.h"
#include <stdio.h>
#include <stdlib.h>
#include "common_functions.h"
//#include "math.h"
//#include "Constants.h"

void calc_body_velosity ( struct Parameters *params, double *v, double time_step, int *status,
struct Conservative_vector *conservative_vectors, struct TimeMoment time_mom );

void moving_the_body ( double body_velocity, double dt, double *coordinate_of_left_boundary_of_body,
                      double *coordinate_of_right_boundary_of_body );

void calc_pressure_for_movement ( struct Parameters *params, struct Primitive_vector *primitive,
                                 int dir, double body_velocity );

#endif /* __BODY_DYNAMIC_H */