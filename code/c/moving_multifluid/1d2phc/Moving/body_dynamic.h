#ifndef __BODY_DYNAMIC_H
#define __BODY_DYNAMIC_H

#include "structures.h"
#include <stdio.h>
#include <stdlib.h>
#include "struct.h"
#include "common_functions.h"
//#include "math.h"
//#include "Constants.h"

double calc_body_velosity (struct ParametersCommon *paramsc, struct Parameters1d *params1d, double body_velocity, double time_step, int *status,
double cont_left[M], double cont_right[M]);

void moving_the_body (struct Parameters1d *params1d, double body_velocity, double dt, double *coordinate_of_left_boundary_of_body,
                      double *coordinate_of_right_boundary_of_body );

void calc_pressure_for_movement (const struct ParametersCommon *paramsc, struct Parameters1d *params1d,double primitive[M],
                                 int dir, double body_velocity );

#endif /* __BODY_DYNAMIC_H */