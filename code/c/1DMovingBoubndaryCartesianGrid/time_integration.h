#ifndef __TIME_INTEGRATION_H
#define __TIME_INTEGRATION_H

#include "structures.h"
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "common_functions.h"
#include "flux.h"
#include "initial_and_boundary_data.h"
#include "memory.h"
#include "body_dynamic.h"

void is_time_evolution( const struct Parameters* params, struct TimeMoment time_mom, bool *condition_to_enter_the_cycle );

void calc_time_step ( struct Parameters *params, int *status, struct Conservative_vector *conservative, double body_velocity, double *time_step );

void time_integration ( struct Parameters *params, struct Conservative_vector *conservative, double time_step,
                       int *status, double body_velocity, double *coordinate_of_left_boundary_of_body,
                      double *coordinate_of_right_boundary_of_body );

void calc_flux_through_edge_of_ghost_cell ( struct Parameters *params, struct Conservative_vector conservative, 
                                           double body_velocity, enum Direction dir, struct Flux_vector *flux,
                                           struct Conservative_vector *Riemann_conservative_solution );

void calc_conservative_on_next_time_moment ( struct Parameters *params, struct Conservative_vector *conservative,
                                    struct Flux_vector left_flux, struct Flux_vector right_flux, double time_step );

void calc_status ( Parameters *params, double coordinate_of_left_boundary_of_body, double coordinate_of_right_boundary_of_body,
    struct Conservative_vector contact_discontinuity_conservative_on_left_border,
    struct Conservative_vector contact_discontinuity_conservative_on_right_border,
    int *status, struct Conservative_vector *conservative, int index_of_left_ghost_cell, int index_of_right_ghost_cell );

#endif /* __TIME_INTEGRATION_H */