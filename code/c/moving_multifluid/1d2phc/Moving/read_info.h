#ifndef __READ_INFO_H_
#define __READ_INFO_H_

#include "string.h"
#include <stdio.h>
#include "constants.h"
#include <stdlib.h>
#include "structures.h"
#include "math.h"

void read_parameters( struct Parameters *params );

void dimless_parameters( struct Parameters *params );

void read_last_time_moment( struct Parameters *params, struct TimeMoment *time_mom );

void read_solution ( struct Parameters *params, struct Conservative_vector *conservative, struct TimeMoment *time_mom );

#endif /* __READ_INFO_H_ */