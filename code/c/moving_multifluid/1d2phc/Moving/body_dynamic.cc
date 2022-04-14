#include "body_dynamic.h"
#include "utils.h"
double calc_body_velosity (struct ParametersCommon *paramsc,  struct Parameters1d *params1d, double v, double time_step, int *status,
double cont_left[M], double cont_right[M])
{
    double body_velocity;
    body_velocity = v;
    double primitive[M];
    double left_pressure = 0.0, right_pressure = 0.0;
    left_pressure = cont_left[P_GAS] * (1 - cont_left[B_DISP]) + cont_left[P_DISP] * cont_left[B_DISP];
    right_pressure = cont_right[P_GAS];
    if ((left_pressure - right_pressure) > params1d->resistive_pressure)
        body_velocity += time_step * params1d->body_cross_section_divided_to_mass * ( left_pressure - right_pressure - params1d->resistive_pressure  ); 
     printf("\n %lf %lf %lf %lf %lf", left_pressure, right_pressure, params1d->body_cross_section_divided_to_mass, body_velocity, primitive[V_GAS]);
   return body_velocity;
}

// определение координат границы тела по прошествии шага по времени
void moving_the_body (struct Parameters1d *params1d, double body_velocity, double dt, double *coordinate_of_left_boundary_of_body,
                      double *coordinate_of_right_boundary_of_body )
{printf("\n %lf", body_velocity);
    *coordinate_of_left_boundary_of_body += body_velocity * dt;
    *coordinate_of_right_boundary_of_body += body_velocity * dt;
    if (*coordinate_of_right_boundary_of_body >= params1d->right_boundary_x){
        printf("\n end of calculation region");
        exit(0);
    }
    printf("\n%lf", *coordinate_of_left_boundary_of_body);
}

//void calc_pressure_for_movement (const struct ParametersCommon *paramsc, struct Parameters1d *params1d,double primitive[M],
//                                 int dir, double body_velocity )
//{
//    
//    
//    double c1;
//    c1 = calc_sound_velocity_one_phase(paramsc, primitive, GAS_PHASE);   
//    double vel_diff = body_velocity - primitive[V_GAS];
//    printf("\n%lf dif %lf", vel_diff, c);
//    if ( fabs(vel_diff) < paramsc->eps_general ) {}
//    else if ( ( vel_diff > 0.0 && 1 == dir ) || ( vel_diff < 0.0 && -1 == dir ) ) // две волны разряжения
//    {
//        printf("check");
//        if ( -1 == dir ) vel_diff *= -1.0;
//        primitive[P_GAS] *= pow( 1.0 - 0.5 * ( paramsc->g2 - 1.0 ) * vel_diff / c, ( 2.0 * paramsc->g2 / ( paramsc->g2 - 1.0 ) ) );
//    }
//    else // две ударные волны
//    {
//        primitive[P_GAS] += 0.25 * ( paramsc->g2 + 1.0 ) * primitive[R_GAS] * vel_diff * vel_diff *
//            ( 1.0 + sqrt( 1.0 + 16.0 * c * c / ( vel_diff * vel_diff * ( paramsc->g2 + 1.0 ) * ( paramsc->g2 + 1.0 ) ) ) );
//    }
//        double c2;
//    c2 = calc_sound_velocity_one_phase(paramsc, primitive, DISPERSED_PHASE);   
//    double vel_diff = body_velocity - primitive[V_GAS];
//    printf("\n%lf dif %lf", vel_diff, c);
//    if ( fabs(vel_diff) < paramsc->eps_general ) {}
//    else if ( ( vel_diff > 0.0 && 1 == dir ) || ( vel_diff < 0.0 && -1 == dir ) ) // две волны разряжения
//    {
//        printf("check");
//        if ( -1 == dir ) vel_diff *= -1.0;
//        primitive[P_GAS] *= pow( 1.0 - 0.5 * ( paramsc->g2 - 1.0 ) * vel_diff / c, ( 2.0 * paramsc->g2 / ( paramsc->g2 - 1.0 ) ) );
//    }
//    else // две ударные волны
//    {
//        primitive[P_GAS] += 0.25 * ( paramsc->g2 + 1.0 ) * primitive[R_GAS] * vel_diff * vel_diff *
//            ( 1.0 + sqrt( 1.0 + 16.0 * c * c / ( vel_diff * vel_diff * ( paramsc->g2 + 1.0 ) * ( paramsc->g2 + 1.0 ) ) ) );
//    }
//}

