#ifndef __FLUX_H
#define __FLUX_H

#include "structures.h"
#include "common_functions.h"

void solve_Riemann_problem ( struct Parameters *params, struct Conservative_vector left_conserative,
struct Conservative_vector right_conserative, struct Primitive_vector *solution_primitive,
struct Primitive_vector *contact_discontinuity_primitive );

void godunov_flux( struct Parameters *params, struct Conservative_vector left_conserative, struct Conservative_vector right_conserative, struct Flux_vector *flux );

void calc_contact_pressure_velocity( struct Parameters *params, struct Primitive_vector left_primitive,
struct Primitive_vector right_primitive, double cl, double cr, double *p_cont, double *v_cont );

void pressure_initial_guess( struct Parameters *params, struct Primitive_vector left_primitive,
struct Primitive_vector right_primitive, double cl, double cr, double *p_guess );

void calc_function_and_derivative( struct Parameters *params, double curr_press, struct Primitive_vector primitive, double c, double *f, double *deriv );

void sample( struct Parameters *params, struct Primitive_vector left_primitive, struct Primitive_vector right_primitive,
			double cl, double cr, double p_cont, double v_cont, double *contact_discontinuity_density,
                        double s, struct Primitive_vector *solution_primitive );

void diff_flux_ncons( struct Parameters *params, struct Primitive_vector primitive, struct Flux_vector *flux );

#endif /* __FLUX_H */