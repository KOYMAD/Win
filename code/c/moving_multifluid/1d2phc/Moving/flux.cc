#include "flux.h"

// ��������� �������� - ����������� ���������� �� ���������� ������� ����� ��� ������
void solve_Riemann_problem ( struct Parameters *params, struct Conservative_vector left_conserative, 
struct Conservative_vector right_conserative, struct Primitive_vector *solution_primitive,
struct Primitive_vector *contact_discontinuity_primitive )
{
    double cl, cr;              /* �������� ����� ����� � ������ �� ������� */

    struct Primitive_vector left_primitive, right_primitive;
    calc_primitive_variables ( params, left_conserative, &left_primitive );
    calc_primitive_variables ( params, right_conserative, &right_primitive );
    calc_sound_velocity ( params, left_primitive, &cl );
    calc_sound_velocity ( params, right_primitive, &cr );

    if ( 2.0 * ( cl + cr ) / ( params->g - 1.0 ) <= right_primitive.velosity - left_primitive.velosity )
    {
        // ������ ������������� �������
        printf( "Error: vacuum is generated.\n" );
	system ( "Pause" );
    }

    // ������������ ��������� ������� �������� � �������� ���� �� ���������� �������
    calc_contact_pressure_velocity( params, left_primitive, right_primitive, cl, cr,
        &contact_discontinuity_primitive->pressure, &contact_discontinuity_primitive->velosity );

    // ����� �������
    sample( params, left_primitive, right_primitive, cl, cr, contact_discontinuity_primitive->pressure,
        contact_discontinuity_primitive->velosity, &contact_discontinuity_primitive->density, 0.0,
        solution_primitive );
}

/*  ����� �������� ������� �������
    left_params[M] - ������ �������������� ���������� ����� �� ������� (in)
    right_params[M] - ������ �������������� ���������� ������ �� ������� (in)
    flux[M] - ������������ ������ ������ (out) */
void godunov_flux( struct Parameters *params, struct Conservative_vector left_conserative, struct Conservative_vector right_conserative, struct Flux_vector *flux )
{
    double cl, cr;                      /* �������� ����� ����� � ������ �� ������� */
    double p_cont, v_cont;              /* �������� � �������� �� ���������� ������� */

    struct Primitive_vector left_primitive, right_primitive, solution_primitive;
    calc_primitive_variables ( params, left_conserative, &left_primitive );
    calc_primitive_variables ( params, right_conserative, &right_primitive );
    calc_sound_velocity ( params, left_primitive, &cl );
    calc_sound_velocity ( params, right_primitive, &cr );

    if ( 2.0 * ( cl + cr ) / ( params->g - 1.0 ) <= right_primitive.velosity - left_primitive.velosity )
    {
        // ������ ������������� �������
        printf( "Error: vacuum is generated.\n" );
	system ( "Pause" );
    }

    // ������������ ��������� ������� �������� � �������� ���� �� ���������� �������
    calc_contact_pressure_velocity( params, left_primitive, right_primitive, cl, cr, &p_cont, &v_cont );

    double contact_discontinuity_density; // ��������� ����� ��� ������ �� ����������� �������

    // ����� �������
    sample( params, left_primitive, right_primitive, cl, cr, p_cont, v_cont,
        &contact_discontinuity_density, 0.0, &solution_primitive );

    // ������ ������ �������� �� ������� ���������������� ����������
    diff_flux_ncons( params, solution_primitive, flux );
}

/* ������������ ��������� ������� �������� � �������� ���� �� ���������� �������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 155. - Subroutine STARPU.

   v_ncons_l[M] - ������ ����������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ����������� ���������� ������ �� ������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   p_cont - �������� �� ���������� ������� (out)
   v_cont - �������� �� ���������� ������� (out) */
void calc_contact_pressure_velocity( struct Parameters *params, struct Primitive_vector left_primitive,
struct Primitive_vector right_primitive, double cl, double cr, double *p_cont, double *v_cont )
{
    double vl, vr;      /* �������� ����� � ������ �� ������� */
    double p_old;       /* �������� �������� �� ���������� �������� */
    double fl, fr;      /* �������� ������� */
    double fld, frd;    /* �������� ����������� */
    int iter_num = 0;   /* ���������� ����������� �������� */
    double criteria;    /* ���������� ��� ����������� ���������� */

    vl = left_primitive.velosity;
    vr = right_primitive.velosity;

    /* ������ ���������� ����������� ��� �������� */
    pressure_initial_guess( params, left_primitive, right_primitive, cl, cr, &p_old );

    /* ������� ����������� ��������� ��� ���������� �������� �� ���������� ������� ������� �������-������� */
    do {
        calc_function_and_derivative( params, p_old, left_primitive, cl, &fl, &fld );
        calc_function_and_derivative( params, p_old, right_primitive, cr, &fr, &frd );
        *p_cont = p_old - ( fl + fr + vr - vl ) / ( fld + frd );
        criteria = 2.0 * fabs( ( *p_cont - p_old ) / ( *p_cont + p_old ) );
        iter_num++;
	if ( iter_num > params->max_iter_num )
	{
            printf( "Error in Godunov flux function calculation. Number of iterations exceeds the maximum value.\n" );
            system ( "Pause" );
        }
        if ( *p_cont < 0.0 )
	{
            printf( "Error in Godunov flux function calculation. Pressure is negative.\n" );
            system ( "Pause" );         
        }
        p_old = *p_cont;
	} while ( criteria > params->eps_general );

    // �������� ����������� �������
    *v_cont = 0.5 * ( vl + vr + fr - fl );
}

/* ����������� ���������� ����������� ��� ������� ��������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   v_ncons_l[M] - ������ ����������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ����������� ���������� ������ �� ������� (in)
   p_guess - ������� ��������� ����������� (out) */
void pressure_initial_guess( struct Parameters *params, struct Primitive_vector left_primitive,
struct Primitive_vector right_primitive, double cl, double cr, double *p_guess )
{
    double rl, vl, pl;                  /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;                  /* ����������� ���������� ������ �� ������� */
    /* ��������� �����������, ������������ �� ����������� ������������ ��������������� �������
       � ����������� ���������� */
    double p_lin;
    double p_min, p_max;                /* ����������� � ������������ �������� ����� � ������ �� ������� */
    double p_ratio;                     /* ������� �� �������� ����� � ������ �� ������� */
    double plr, v_star, p1, p2, g1, g2; /* ��������������� ���������� ��� ������������� �������� */

    /* �������� ����������� ��� �������� */
    rl = left_primitive.density;
    vl = left_primitive.velosity;
    pl = left_primitive.pressure;
    rr = right_primitive.density;
    vr = right_primitive.velosity;
    pr = right_primitive.pressure;

    /* ��������� ����������� �� �������� ������
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = min( pl, pr );
    p_max = max( pl, pr );
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= params->p_max_ratio ) && ( p_min <= p_lin && p_lin <= p_max ) )
    {
        /* ��������� ����������� �� ��������������� ������ */
        *p_guess = p_lin;
    }
    else
    {
        if ( p_lin < p_min )
	{
            /* ��������� ����������� �� ���� ������ ����������
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 302. - Formulas (9.35), (9.36). */
	    plr = pow( pl / pr, 0.5 * ( params->g - 1.0 ) / params->g );
            v_star = ( plr * vl / cl + vr / cr + 2.0 * ( plr - 1.0 ) / ( params->g - 1.0 ) ) / ( plr / cl + 1.0 / cr );
            g1 = 0.5 * ( params->g - 1.0 );
            p1 = 1.0 + g1 * ( vl - v_star ) / cl;
            p2 = 1.0 + g1 * ( v_star - vr ) / cr;
            g2 = params->g / g1;
            *p_guess = 0.5 * ( pl * pow( p1, g2 ) + pr * pow( p2, g2 ) );
        }
	else
	{
            /* ��������� ����������� �� ���� ������� ������
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48).*/
            g1 = 2.0 / ( params->g + 1.0 );
            g2 = ( params->g - 1.0 ) / ( params->g + 1.0 );
            p1 = sqrt( g1 / rl / ( g2 * pl + p_lin ) );
            p2 = sqrt( g1 / rr / ( g2 * pr + p_lin ) );
            *p_guess = ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
        }
    }
}

/* ������ ������� � �� ����������� ��� ������������� ��������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.

   curr_press - �������� � ���������� �������� (in)
   v_ncons[M] - ������ ����������� ���������� (in)
   c - �������� �����
   f - �������� ������� (out)
   deriv - �������� ����������� (out) */
void calc_function_and_derivative( struct Parameters *params, double curr_press, struct Primitive_vector primitive, double c, double *f, double *deriv )
{
    double r, v, p;                 /* ����������� ���������� */
    double p_ratio, fg, a, b, q;    /* ��������������� ���������� */

    /* �������� ����������� ��� �������� */
	r = primitive.density;
	v = primitive.velosity;
	p = primitive.pressure;
    
    if ( curr_press <= p )
	{
        /* ����� ���������� */
        p_ratio = curr_press / p;
		fg = 2.0 / ( params->g - 1.0 );
        *f = fg * c * ( pow( p_ratio, 1.0 / fg / params->g ) - 1.0 );
        *deriv = ( 1.0 / r / c ) * pow( p_ratio, - 0.5 * ( params->g + 1.0 ) / params->g );
    }
    else
	{
        /* ������� ����� */
        a = 2.0 / r / ( params->g + 1.0 );
        b = ( params->g - 1.0 ) * p / ( params->g + 1.0 );
        q = sqrt( a / ( b + curr_press ) );
        *f = ( curr_press - p ) * q;
        *deriv = ( 1.0 - 0.5 * ( curr_press - p ) / ( b + curr_press ) ) * q;
    }
}

/* ������� ������ �������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE. 

   v_ncons_l[M] - ������ ���������������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ���������������� ���������� ������ �� ������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   p_cont - �������� �� ���������� ������� (in)
   v_cont - �������� �� ���������� ������� (in)
   s - �������� x/t, ��� �������� ���������� ������� (in)
   v_ncons_res[M] - ���������� ������ ���������������� ���������� (out) */
void sample( struct Parameters *params, struct Primitive_vector left_primitive, struct Primitive_vector right_primitive,
			double cl, double cr, double p_cont, double v_cont, double *contact_discontinuity_density,
                        double s, struct Primitive_vector *solution_primitive )
{
    double rl, vl, pl;  /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;  /* ����������� ���������� ������ �� ������� */

    /* �������� ����� ���� */
    double shl, stl;    /* �������� "������" � "������" ����� ����� ���������� */
    double sl;          /* �������� ����� ������� ����� */

    /* �������� ������ ���� */
    double shr, str;    /* �������� "������" � "������" ������ ����� ���������� */
    double sr;          /* �������� ������ ������� ����� */

    double cml, cmr;    /* �������� ����� ����� � ������ �� ����������� ������� */
    double c;           /* ��������� �������� ����� ������ ����� ���������� */
    double p_ratio;
    double r, v, p;     /* ���������� �������� ���������, �������� � �������� */

    /* �������� ����������� ��� �������� ������ */
    rl = left_primitive.density;
    vl = left_primitive.velosity;
    pl = left_primitive.pressure;
    rr = right_primitive.density;
    vr = right_primitive.velosity;
    pr = right_primitive.pressure;

    double GAMMA = params->g;
    double G1 = 0.5 * ( GAMMA - 1.0 ) / GAMMA;
    double G2 = 0.5 * ( GAMMA + 1.0 ) / GAMMA;
    double G3 = 2.0 * GAMMA / ( GAMMA - 1.0 );
    double G4 = 2.0 / ( GAMMA - 1.0 );
    double G5 = 2.0 / ( GAMMA + 1.0 );
    double G6 = ( GAMMA - 1.0 ) / ( GAMMA + 1.0 );
    double G7 = 0.5 * ( GAMMA - 1.0 );
    double G8 = GAMMA - 1.0;

    if ( s <= v_cont ) {
        /* ��������������� ����� - ����� �� ����������� ������� */
        if ( p_cont <= pl ) {
            /* ����� ����� ���������� */
            /* ��������� ����� �� ����������� ������� */
            *contact_discontinuity_density = rl * pow( p_cont / pl, 1.0 / GAMMA );
            shl = vl - cl;
            if ( s <= shl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( p_cont / pl, G1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* ��������� ����� �� ����������� ������� */
                    r = rl * pow( p_cont / pl, 1.0 / GAMMA );
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* ��������� ������ ����� ����� ���������� */
                    v = G5 * ( cl + G7 * vl + s );
                    c = G5 * ( cl + G7 * ( vl - s ) );
                    r = rl * pow( c / cl, G4 );
                    p = pl * pow( c / cl, G3 );
                }
            }
        }
        else {
            /* ����� ������� ����� */
            p_ratio = p_cont / pl;
            /* ��������� �� ����� ������� ������ */
            *contact_discontinuity_density = rl * ( p_ratio + G6 ) / ( p_ratio * G6 + 1.0 );
            sl = vl - cl * sqrt( G2 * p_ratio + G1 );
            if ( s <= sl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* ��������� �� ����� ������� ������ */
                r = rl * ( p_ratio + G6 ) / ( p_ratio * G6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        /* ��������������� ����� - ������ �� ����������� ������� */
        if ( p_cont > pr ) {
            /* ������ ������� ����� */
            p_ratio = p_cont / pr;
            /* ��������� �� ������ ������� ������ */
            *contact_discontinuity_density = rr * ( p_ratio + G6 ) / ( p_ratio * G6 + 1.0 );
            sr = vr + cr * sqrt( G2 * p_ratio + G1 );
            if ( s >= sr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* ��������� �� ������ ������� ������ */
                r = rr * ( p_ratio + G6 ) / ( p_ratio * G6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* ������ ����� ���������� */
            /* ��������� ������ �� ����������� ������� */
            *contact_discontinuity_density = rr * pow( p_cont / pr, 1.0 / GAMMA );
            shr = vr + cr;
            if ( s >= shr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( p_cont / pr, G1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* ��������� ������ �� ����������� ������� */
                   r = rr * pow( p_cont / pr, 1.0 / GAMMA );
                   v = v_cont;
                   p = p_cont;
               }
               else {
                    /* ��������� ������ ������ ����� ���������� */
                    v = G5 * ( -cr + G7 * vr + s );
                    c = G5 * ( cr - G7 * ( vr - s ) );
                    r = rr * pow( c / cr, G4 );
                    p = pr * pow( c / cr, G3 );
               }
            }
        }
    }
    
    /* ������������ ��������� ������� � ����������� */
    solution_primitive->density = r;
    solution_primitive->velosity = v;
    solution_primitive->pressure = p;
}// ������� ������ �������

/* ������ ������� ����������������� ������ �� ������� ����������� ����������
   v_ncons[M] - ������ ����������� ���������� (in)
   flux[M] - �������������� ������ ����������������� ������ (out) */
void diff_flux_ncons( struct Parameters *params, struct Primitive_vector primitive, struct Flux_vector *flux )
{
    struct Conservative_vector conservative;
    calc_conservative_variables ( params, primitive, &conservative );

    flux->mass_flux = conservative.specific_momentum;                                    /* ����� */
    flux->momentum_flux = conservative.specific_momentum * primitive.velosity + primitive.pressure;          /* ������� */
    flux->energy_flux = ( conservative.specific_energy + primitive.pressure ) * primitive.velosity;      /* ������ ������� */
}