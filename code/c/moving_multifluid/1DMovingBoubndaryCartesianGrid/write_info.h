#ifndef __WRITE_INFO_H
#define __WRITE_INFO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structures.h"
#include "common_functions.h"


void write_results ( struct Parameters *params, int *status, struct Conservative_vector *conservative, struct TimeMoment time_mom );

void write_piston_trajectory ( struct Parameters *params, double coordinate_of_left_boundary_of_body, struct TimeMoment time_mom );

#endif /* __WRITE_INFO_H */