#include "time_integration.h"

// ����������, ���������� ������ ��� ���������� ��������� �����/���������� �����
// params - ��������� � ����������� ��������������� ������������
// time_iter_num - ������� ���������� ����� �� �������
// curr_time - ������� ������ �������
// ���������� true, ���� ����� ���������� ������; false - �����
void is_time_evolution( const struct Parameters* params, struct TimeMoment time_mom, bool *condition_to_enter_the_cycle )
{
    switch( params->exit_time_cycle )
    {
        case ITERATIONS:
            // �������� �������� - �������� ���������� ����� �� �������
            if ( time_mom.steps_num < params->time_cycle_iterations )
                *condition_to_enter_the_cycle = true;
            else
                *condition_to_enter_the_cycle = false;
            break;
        case FINAL_TIME:
            // �������� �������� - �������� ������ �������
            if ( time_mom.curr_t < params->t_fin )
                *condition_to_enter_the_cycle = true;
            else
                *condition_to_enter_the_cycle = false;
            break;
        default:
            printf( "\nis_time_evolution -> wrong exit_time_cycle value\n\n" );
            system ( "Pause" );
    }
}

//������� �������� ��� �� ������� � ���������� ���
//params - ��������� � ����������� ��������������� ������������ (in)
//*conservative - ������ �������������� ���������� (out)
// status - ������ �������� ��������� �����
void calc_time_step ( struct Parameters *params, int *status, struct Conservative_vector *conservative, double body_velocity, double *time_step )
{
    double char_time; // ����� �� ������� ���������� �� ������ ���������������� �� ���� ������
    switch ( params->time_step_method )
    {
        case DYNAMIC_TIME_STEP:
        {
            double delta = ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary ) / params->number_of_cells;
            double c; // �������� ����� � ������� ������
            for ( int i = 0 ; i < params->number_of_cells; i++ )
            {
                struct Primitive_vector primitive;
                if ( INNER == status[i] || BOUNDARY == status[i] )
                {
                    calc_primitive_variables( params, conservative[i], &primitive );
                    calc_sound_velocity ( params, primitive, &c ); // ������ �������� ����� � ������� ������
                    char_time = ( delta / ( abs ( primitive.velosity ) + c ) );	
                    if ( char_time < *time_step ) *time_step = char_time;
                }
            }
            if ( (*time_step) * params->cfl > delta / abs ( body_velocity ) )
            {
	        printf( "\ncalc_time_step -> too big time step or body velocity\n" );
	        system ( "Pause" );
	    }
            break;
        }
        case CONSTANT_TIME_STEP:
            *time_step = params->dt;
            break;
        default:
            system ( "Pause" );
    }
    *time_step *= params->cfl;
}

void time_integration ( struct Parameters *params, struct Conservative_vector *conservative,
                       double time_step, int *status, double body_velocity, double *coordinate_of_left_boundary_of_body,
                      double *coordinate_of_right_boundary_of_body )
{
    struct Flux_vector *array_of_fluxes; // ������ ������� ����� ����
    get_memory_for_1D_flux_vector_array ( params->number_of_cells + 1, &array_of_fluxes );

    // �������������� ���������� ����� � ������ �� ������ ��������� �������
    struct Conservative_vector behind_left_boundary_conservative_vector;
    struct Conservative_vector behind_right_boundary_conservative_vector;
    // ���� ��������� ������� - ������ �������������� ���������� ����� � ������ �� ��������� �������
    boundary_conditions ( params, conservative, &behind_left_boundary_conservative_vector,
        &behind_right_boundary_conservative_vector );
    // ������ ������� �� ����� � ������ �������� ��������� �������
    godunov_flux ( params, behind_left_boundary_conservative_vector, conservative[0], &(array_of_fluxes[0]) );
    godunov_flux ( params, conservative[params->number_of_cells - 1], behind_right_boundary_conservative_vector,
        &(array_of_fluxes[params->number_of_cells]) );
    // ������ �������������� ���������� �� ���������� ������� �� ����� ������� ���������� ����
    struct Conservative_vector contact_discontinuity_conservative_on_left_border;
    // ������ �������������� ���������� �� ���������� ������� �� ������ ������� ���������� ����
    struct Conservative_vector contact_discontinuity_conservative_on_right_border;
    int index_of_left_ghost_cell, index_of_right_ghost_cell; // ������� ����� � ������ ����������� �����

    // ������ �������
    for ( int i = 0 ; i < params->number_of_cells - 1 ; i++ ) // ���� �� ������ ����� ����� ������ �� �������������
    {
        if ( INNER == status[i] ) // ���� ������ ����������
            // ��������� ������� ������ - ��������� � ����� � ������ ������� ��������������
            godunov_flux ( params, conservative[i], conservative[i+1], &array_of_fluxes[i+1] );
        else if ( BOUNDARY == status[i] ) // ���� ������ ���������
        {
            if ( INNER == status[i+1] ) // ���� ������ ������ ����������
                // ��������� ������� ������ - ��������� � ����� � ������ ������� ��������������
                godunov_flux ( params, conservative[i], conservative[i+1], &array_of_fluxes[i+1] );
            else if ( GHOST == status[i+1] ) // ���� ������ ������ �����������
            {
                index_of_left_ghost_cell = i + 1;
                // � �������� ���������� ������ ���������� ��������, ��������� �� �������� ������� ����
	        calc_flux_through_edge_of_ghost_cell ( params, conservative[i], body_velocity, FROM_LEFT_TO_RIGHT,
                    &array_of_fluxes[i+1], &contact_discontinuity_conservative_on_left_border );
            }
        }
        else if ( GHOST == status[i] ) // ���� ������ �����������
        {
            if ( BOUNDARY == status[i+1] ) // ���� ������ ������ ���������
            {
                index_of_right_ghost_cell = i;
                // � �������� ���������� ����� ���������� ��������, ��������� �� �������� ������� ����
                calc_flux_through_edge_of_ghost_cell ( params, conservative[i+1], body_velocity, FROM_RIGHT_TO_LEFT,
                    &array_of_fluxes[i+1], &contact_discontinuity_conservative_on_right_border );
            }
            else continue;
        }
	else if ( OUTER == status[i] ) {} // ���� ������ �������, ������ �� � ����� �� ��������������
	else
	{
	    printf( "\ntime_integration -> wrong status value\n\n" );
	    system ( "Pause" );
	}
    }

    // ���������� �������������� ���������� �� ����� ���� �� �������
    for ( int i = 0 ; i < params->number_of_cells ; i++ ) // ���� �� �������
    {
        if ( INNER == status[i] || BOUNDARY == status[i] ) // ���� ������ ���������� ��� ���������
            calc_conservative_on_next_time_moment ( params, &(conservative[i]), array_of_fluxes[i],
            array_of_fluxes[i+1], time_step );
        else continue;
    }

    // ������������ ���� � ����������� ����� ��������� ��� ����� �������
    moving_the_body ( body_velocity, time_step, coordinate_of_left_boundary_of_body,
        coordinate_of_right_boundary_of_body );

    // �������� �������� ����� � ������������� ����� � ����� ��������
    calc_status ( params, *coordinate_of_left_boundary_of_body, *coordinate_of_right_boundary_of_body,
        contact_discontinuity_conservative_on_left_border, contact_discontinuity_conservative_on_right_border,
        status, conservative, index_of_left_ghost_cell, index_of_right_ghost_cell );
}

// ������ ������ ����� �����, ����������� ��������� ������ � �����������, ���������� � ���� ������� ���������� ����
void calc_flux_through_edge_of_ghost_cell ( struct Parameters *params, struct Conservative_vector conservative,
                                           double body_velocity, enum Direction dir, struct Flux_vector *flux, 
                                           struct Conservative_vector *contact_discontinuity_conservative )
{
    struct Primitive_vector Riemann_primitive_solution; // ������ ����������� ���������� ������� ������ ������
    struct Conservative_vector ghost_conservative; // ������ �������������� ���������� � ����������� ������
    struct Primitive_vector ghost_primitive; // ������ ����������� ���������� � ����������� ������
    struct Primitive_vector current_primitive; // ������ ����������� ���������� � ������� ��������� ������
    calc_primitive_variables ( params, conservative, &current_primitive );
    ghost_primitive = current_primitive;
    ghost_primitive.velosity = 2 * body_velocity - current_primitive.velosity;
    calc_conservative_variables ( params, ghost_primitive, &ghost_conservative );
    struct Primitive_vector contact_discontinuity_primitive; // ����������� ���������� �� ���������� �������
    // ��������������� ������� ������ ������ � ����������� �������� �� ���������� �������
    if ( FROM_LEFT_TO_RIGHT == dir )
    {
        solve_Riemann_problem ( params, conservative, ghost_conservative, &Riemann_primitive_solution,
            &contact_discontinuity_primitive );
    }
    else if ( FROM_RIGHT_TO_LEFT == dir )
    {
        solve_Riemann_problem ( params, ghost_conservative, conservative, &Riemann_primitive_solution,
            &contact_discontinuity_primitive );
    }
    // ������ �������������� ���������� �� ���������� �������
    calc_conservative_variables ( params, contact_discontinuity_primitive, contact_discontinuity_conservative );
    //calc_conservative_variables ( params, Riemann_primitive_solution, Riemann_conservative_solution );
    // �������� ����� �� ���������� �������
    double contact_discontinuity_sound_velocity;
    calc_sound_velocity ( params, contact_discontinuity_primitive, &contact_discontinuity_sound_velocity );
    // ���� ���� �������� ����� ������� �������� ����� �� ���������� �������
    if ( body_velocity < -contact_discontinuity_sound_velocity )
    {
        if ( FROM_LEFT_TO_RIGHT == dir ) diff_flux_ncons( params, Riemann_primitive_solution, flux );
        else if ( FROM_RIGHT_TO_LEFT == dir ) diff_flux_ncons( params, current_primitive, flux );
    }
    // ���� ���� �������� ��������� �������� ����� �� ���������� �������
    else if ( body_velocity < contact_discontinuity_sound_velocity )
    {
        if ( FROM_LEFT_TO_RIGHT == dir ) godunov_flux ( params, conservative, *contact_discontinuity_conservative, flux );
        else if ( FROM_RIGHT_TO_LEFT == dir ) godunov_flux ( params, *contact_discontinuity_conservative, conservative, flux );
        else
        {
            printf( "\ncalc_flux_through_edge_of_ghost_cell -> wrong Direction value\n\n" );
            system ( "Pause" );
        }
    }
    // ���� ���� �������� ������ ������� �������� ����� �� ���������� �������
    else
    {
        if ( FROM_LEFT_TO_RIGHT == dir ) diff_flux_ncons( params, current_primitive, flux );
        else if ( FROM_RIGHT_TO_LEFT == dir ) diff_flux_ncons( params, Riemann_primitive_solution, flux );
    }
}

void calc_conservative_on_next_time_moment ( struct Parameters *params, struct Conservative_vector *conservative,
                                    struct Flux_vector left_flux, struct Flux_vector right_flux, double time_step )
{
    double grid_step = ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary ) /
        params->number_of_cells; // ��� �����

    conservative->specific_mass -= ( right_flux.mass_flux - left_flux.mass_flux ) * time_step / grid_step;
    conservative->specific_momentum -= ( right_flux.momentum_flux - left_flux.momentum_flux ) * time_step / grid_step;
    conservative->specific_energy -= ( right_flux.energy_flux - left_flux.energy_flux ) * time_step / grid_step;
}

// int i (in) - ������ ����������� ������ (���������� ������� ����)
//void calc_status ( Parameters *params, double coordinate_of_left_boundary_of_body, double coordinate_of_right_boundary_of_body,
//    struct Conservative_vector contact_discontinuity_conservative_on_left_border,
//    struct Conservative_vector contact_discontinuity_conservative_on_right_border,
//    int *status, struct Conservative_vector *conservative, int index_of_left_ghost_cell, int index_of_right_ghost_cell )
//{
//    struct Conservative_vector null_conservative = {0.0, 0.0, 0.0};
//    double grid_step = ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary )
//        / params->number_of_cells;
//
//    // ����� ������� ���� ���������� ������� � ������ i + 1
//    if ( coordinate_of_left_boundary_of_body >= ( index_of_left_ghost_cell + 1 ) * grid_step )
//    {
//        status[index_of_left_ghost_cell + 1] = GHOST;
//        status[index_of_left_ghost_cell] = BOUNDARY;
//        conservative[index_of_left_ghost_cell] = contact_discontinuity_conservative_on_left_border;
//        status[index_of_left_ghost_cell - 1] = INNER;
//    }
//    else if ( coordinate_of_left_boundary_of_body >= index_of_left_ghost_cell * grid_step ) // �������� � ��� �� ������
//    {}
//    else // ���������� ������ � ������ i - 1
//    {
//        status[index_of_left_ghost_cell - 2] = BOUNDARY;
//        status[index_of_left_ghost_cell - 1] = GHOST;
//        status[index_of_left_ghost_cell] = OUTER;
//        conservative[index_of_left_ghost_cell] = null_conservative;
//    }
//
//    // ������ ������� ���� ���������� ������� � ������ i + 1
//    if ( coordinate_of_right_boundary_of_body >= ( index_of_right_ghost_cell + 1 ) * grid_step ) // ���������� ������� � ������ i + 1
//    {
//        status[index_of_right_ghost_cell + 2] = BOUNDARY;
//        status[index_of_right_ghost_cell + 1] = GHOST;
//        status[index_of_right_ghost_cell] = OUTER;
//        conservative[index_of_right_ghost_cell] = null_conservative;
//    }
//    else if ( coordinate_of_right_boundary_of_body >= index_of_right_ghost_cell * grid_step ) // �������� � ��� �� ������
//    {}
//    else // ���������� ������ � ������ i - 1
//    {
//        status[index_of_right_ghost_cell - 1] = GHOST;
//        status[index_of_right_ghost_cell] = BOUNDARY;
//        conservative[index_of_right_ghost_cell] = contact_discontinuity_conservative_on_right_border;
//        status[index_of_right_ghost_cell + 1] = INNER;
//    }
//}