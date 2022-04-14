#ifndef __INITIAL_AND_BOUNDARY_DATA_H
#define __INITIAL_AND_BOUNDARY_DATA_H

#include <stdio.h>
#include "structures.h"
#include "common_functions.h"
#include "body_dynamic.h"

void initiate_data ( struct Parameters *params, struct Conservative_vector *conservative,
                    int *status, struct TimeMoment *time_mom, int *zones_count );

void initiate_status( struct Parameters params, double left, double right, int *status );

void boundary_conditions ( struct Parameters *params, struct Conservative_vector *conservative, 
struct Conservative_vector *behind_left_boundary_conservative_vector, struct Conservative_vector *behind_right_boundary_conservative_vector );

#endif /* __INITIAL_AND_BOUNDARY_DATA_H */