#include "body_dynamic.h"

void calc_body_velosity ( struct Parameters *params, double *v, double time_step, int *status,
struct Conservative_vector *conservative_vectors, struct TimeMoment time_mom )
{
    struct Primitive_vector primitive_vector;
    double left_pressure = 0.0, right_pressure = 0.0;
    int dir; // направление внутренней нормали к границе тела ( 1 - слева направо, -1 - справа налево)
    for ( int i = 1 ; i < params->number_of_cells - 1 ; i++ )
    {
        if ( BOUNDARY == status[i] && INNER == status[i - 1] )
        {
            calc_primitive_variables ( params, conservative_vectors[i], &primitive_vector );
            calc_pressure_for_movement ( params, &primitive_vector, 1, *v );
            left_pressure = primitive_vector.pressure;
        }
        if ( BOUNDARY == status[i] && INNER == status[i + 1] )
        {
            calc_primitive_variables ( params, conservative_vectors[i], &primitive_vector );
            calc_pressure_for_movement ( params, &primitive_vector, -1, *v );
            right_pressure = primitive_vector.pressure;
        }
    }
    //*v += time_step * params->body_cross_section_devided_to_mass * ( left_pressure - right_pressure );
}

// определение координат границы тела по прошествии шага по времени
void moving_the_body ( double body_velocity, double dt, double *coordinate_of_left_boundary_of_body,
                      double *coordinate_of_right_boundary_of_body )
{
    *coordinate_of_left_boundary_of_body += body_velocity * dt;
    *coordinate_of_right_boundary_of_body += body_velocity * dt;
}

void calc_pressure_for_movement ( struct Parameters *params, struct Primitive_vector *primitive,
                                 int dir, double body_velocity )
{
    double c;
    calc_sound_velocity (params, *primitive, &c);
    double vel_diff = body_velocity - primitive->velosity;
    if ( fabs(vel_diff) < params->eps_general ) {}
    else if ( ( vel_diff < 0.0 && 1 == dir ) || ( vel_diff > 0.0 && -1 == dir ) ) // две волны разрежения
    {
        if ( 1 == dir ) vel_diff *= -1.0;
        primitive->pressure *= pow( 1.0 - 0.5 * ( params->g - 1.0 ) * vel_diff / c, ( 2.0 * params->g / ( params->g - 1.0 ) ) );
    }
    else // две ударные волны
    {
        primitive->pressure += 0.25 * ( params->g + 1.0 ) * primitive->density * vel_diff * vel_diff *
            ( 1.0 + sqrt( 1.0 + 16.0 * c * c / ( vel_diff * vel_diff * ( params->g + 1.0 ) * ( params->g + 1.0 ) ) ) );
    }
}