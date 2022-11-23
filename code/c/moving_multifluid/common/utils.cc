// utils.cc
// ������� �������, ����������� ��� ���������������� ������ ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// ��� ����������� �� ����������� ������
// (c) ����� �����, 2018
// ������: 26 ������� 2018 �.

#include "utils.h"
void initiate_status( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double left, double right, int *status, int *left_ghost, int *right_ghost )
{
    double grid_step = ( params1d->right_boundary_x - params1d->left_boundary_x ) /
        params1d->cells_number; // ��� �����
    for ( int i = 0 ; i < params1d->cells_number ; i++ )
    {
        if ( left >= i * grid_step ) // ���� ����� ������� ���� ������ ����� ������� ������
        {
            if ( left >= ( i + 1 ) * grid_step ) // ���� ����� ������� ���� ������ ������ ������� ������
            {
                if ( right >= i * grid_step ) // ���� ������ ������� ���� ������ ����� ������� ������
                {
                    if ( right >= ( i + 1 ) * grid_step ) // ���� ������ ������� ���� ������ ������ ������� ������
                        status[i] = INNER;
                    else // ���� ������ ������� ���� ����� ������ ������� ������
                    {
                        printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                        system ( "Pause" );
                    }
                }
                else // ���� ������ ������� ���� ����� ����� ������� ������
                {
                    printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                    system ( "Pause" );
                }
            }
            else // ���� ����� ������� ���� ����� ������ ������� ������
            {
                if ( right >= i * grid_step ) // ���� ������ ������� ���� ������ ����� ������� ������
                {
                    if ( right >= ( i + 1 ) * grid_step ){ // ���� ������ ������� ���� ������ ������ ������� ������
                        status[i] = GHOST;
                        *left_ghost = i;
                    }
                    else // ���� ������ ������� ���� ����� ������ ������� ������
                    {
                        printf( "\ninitiate_status -> all body is inside current cell\n" );
                        system ( "Pause" );
                    }
                }
                else // ���� ������ ������� ���� ����� ����� ������� ������
                {
                    printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                    system ( "Pause" );
                }
            }
        }
        else // ���� ����� ������� ���� ����� ����� ������� ������
        {
            if ( left >= ( i + 1 ) * grid_step ) // ���� ����� ������� ���� ������ ������ ������� ������
            {
                printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                system ( "Pause" );
            }
            else // ���� ����� ������� ���� ����� ������ ������� ������
            {
                if ( right >= i * grid_step ) // ���� ������ ������� ���� ������ ����� ������� ������
                {
                    if ( right >= ( i + 1 ) * grid_step ) // ���� ������ ������� ���� ������ ������ ������� ������
                        status[i] = OUTER;
                    else{ // ���� ������ ������� ���� ����� ������ ������� ������
                        status[i] = GHOST;
                        *right_ghost = i;
                    }

                }
                else // ���� ������ ������� ���� ����� ����� ������� ������
                {
                    if ( right >= ( i + 1 ) * grid_step ) // ���� ������ ������� ���� ������ ������ ������� ������
                    {
                        printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                        system ( "Pause" );
                    }
                    else // ���� ������ ������� ���� ����� ������ ������� ������
                        status[i] = INNER;
                }
            }
        }
    }
    for ( int i = 0 ; i < params1d->cells_number ; i++ )
    {
        if ( GHOST == status[i] )
        {
            if ( i != 0 && INNER == status[i - 1] && OUTER == status[i + 1] )
                status[i - 1] = BOUNDARY;
            else if ( i != params1d->cells_number - 1 && INNER == status[i + 1] 
            && OUTER == status[i - 1] )
                status[i + 1] = BOUNDARY;
            else
            {
                printf( "\ndefine_status -> wrong status definition\n" );
                system ( "Pause" );
            }
        }
    }
}

void calc_status ( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double coordinate_of_left_boundary_of_body, double coordinate_of_right_boundary_of_body,
    int *status, int *index_of_left_ghost_cell, int *index_of_right_ghost_cell, double **v_ncons,double v_left[M], double v_right[M] )
{
    double grid_step = ( params1d->right_boundary_x - params1d->left_boundary_x )
        / params1d->cells_number;

    // ����� ������� ���� ���������� ������� � ������ i + 1
    if ( coordinate_of_left_boundary_of_body >= ( *index_of_left_ghost_cell + 1 ) * grid_step )
    {
        printf("\n right \n %d %lf %lf %lf", *index_of_left_ghost_cell, grid_step, coordinate_of_left_boundary_of_body, v_left[P_GAS]);

        status[*index_of_left_ghost_cell + 1] = GHOST;
        status[*index_of_left_ghost_cell] = BOUNDARY;
        status[*index_of_left_ghost_cell - 1] = INNER;
        v_ncons[ *index_of_left_ghost_cell ][P_GAS] = v_left[P_GAS];
        v_ncons[ *index_of_left_ghost_cell ][V_GAS] = v_left[V_GAS];
        v_ncons[ *index_of_left_ghost_cell ][R_GAS] = v_left[R_GAS];
        v_ncons[ *index_of_left_ghost_cell ][P_DISP] = v_left[P_DISP];
        v_ncons[ *index_of_left_ghost_cell ][V_DISP] = v_left[V_DISP];
        v_ncons[ *index_of_left_ghost_cell ][R_DISP] = v_left[R_DISP];
        v_ncons[ *index_of_left_ghost_cell ][B_DISP] =  0.8 * v_left[B_DISP];
        *index_of_left_ghost_cell += 1;
    }
    else if ( coordinate_of_left_boundary_of_body >= *index_of_left_ghost_cell * grid_step ) // �������� � ��� �� ������
    {}
    else // ���������� ������ � ������ i - 1
    {
        status[*index_of_left_ghost_cell - 2] = BOUNDARY;
        status[*index_of_left_ghost_cell - 1] = GHOST;
        status[*index_of_left_ghost_cell] = OUTER;
        v_ncons[ *index_of_left_ghost_cell ][P_GAS] = 0;
        v_ncons[ *index_of_left_ghost_cell ][V_GAS] = 0;
        v_ncons[ *index_of_left_ghost_cell ][R_GAS] = 0;
        *index_of_left_ghost_cell -= 1;
    }

    // ������ ������� ���� ���������� ������� � ������ i + 1
    if ( coordinate_of_right_boundary_of_body >= ( *index_of_right_ghost_cell + 1 ) * grid_step ) // ���������� ������� � ������ i + 1
    {
        status[*index_of_right_ghost_cell + 2] = BOUNDARY;
        status[*index_of_right_ghost_cell + 1] = GHOST;
        status[*index_of_right_ghost_cell] = OUTER;
        v_ncons[ *index_of_right_ghost_cell ][P_GAS] = 1e5;
        v_ncons[ *index_of_right_ghost_cell ][V_GAS] = 0;
        v_ncons[ *index_of_right_ghost_cell ][R_GAS] = 1.2;
        v_ncons[ *index_of_right_ghost_cell ][P_DISP] = 1e5;
        v_ncons[ *index_of_right_ghost_cell ][V_DISP] = 0;
        v_ncons[ *index_of_right_ghost_cell ][R_DISP] = 1578;
        *index_of_right_ghost_cell += 1;
    }
    else if ( coordinate_of_right_boundary_of_body >= *index_of_right_ghost_cell * grid_step ) // �������� � ��� �� ������
    {}
    else // ���������� ������ � ������ i - 1
    {
        status[*index_of_right_ghost_cell - 1] = GHOST;
        status[*index_of_right_ghost_cell] = BOUNDARY;
        status[*index_of_right_ghost_cell + 1] = INNER;
        v_ncons[ *index_of_right_ghost_cell ][P_GAS] = v_right[P_GAS];
        v_ncons[ *index_of_right_ghost_cell ][V_GAS] = v_right[V_GAS];
        v_ncons[ *index_of_right_ghost_cell ][R_GAS] = v_right[R_GAS];
        v_ncons[ *index_of_right_ghost_cell ][P_DISP] = v_right[P_DISP];
        v_ncons[ *index_of_right_ghost_cell ][V_DISP] = v_right[V_DISP];
        v_ncons[ *index_of_right_ghost_cell ][R_DISP] = v_right[R_DISP];
        v_ncons[ *index_of_right_ghost_cell ][B_DISP] = v_right[B_DISP];
        *index_of_right_ghost_cell -= 1;
    }
}
void convert_2D_to_1D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector){
    u_1D_vector[B_DISP] = u_2D_vector[B_DISP_2D];
    u_1D_vector[R_DISP] = u_2D_vector[R_DISP_2D];
    u_1D_vector[P_DISP] = u_2D_vector[P_DISP_2D];
    u_1D_vector[R_GAS] = u_2D_vector[R_GAS_2D];
    u_1D_vector[P_GAS] = u_2D_vector[P_GAS_2D];
    switch (dir){
        case X_DIRECTION:
            u_1D_vector[V_DISP] = u_2D_vector[V_DISP_2D];
            u_1D_vector[V_GAS] = u_2D_vector[V_GAS_2D];
            break;
        case Y_DIRECTION:
            u_1D_vector[V_DISP] = u_2D_vector[U_DISP_2D];
            u_1D_vector[V_GAS] = u_2D_vector[U_GAS_2D];
            break;
        default:
            printf( "\nconvert_2D_to_1D -> wrong direction - only x- or y- are possible\n\n" );
            exit( EXIT_FAILURE );
    }
}

void convert_1D_to_2D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector){
    u_2D_vector[B_DISP_2D] = u_1D_vector[B_DISP];
    u_2D_vector[R_DISP_2D] = u_1D_vector[R_DISP];
    u_2D_vector[P_DISP_2D] = u_1D_vector[P_DISP];
    u_2D_vector[R_GAS_2D] = u_1D_vector[R_GAS];
    u_2D_vector[P_GAS_2D] = u_1D_vector[P_GAS];
    switch (dir){
        case X_DIRECTION:
            u_2D_vector[V_DISP_2D] = u_1D_vector[V_DISP];
            u_2D_vector[V_GAS_2D] = u_1D_vector[V_GAS];
            break;
        case Y_DIRECTION:
            u_2D_vector[U_DISP_2D] = u_1D_vector[V_DISP];
            u_2D_vector[U_GAS_2D] = u_1D_vector[V_GAS];
            break;
        default:
            printf( "\nconvert_2D_to_1D -> wrong direction - only x- or y- are possible\n\n" );
            exit( EXIT_FAILURE );
    }
}


// ������������� �������-������� ��� ���������� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// files_directory - ����������, � ������� ��������� ��� ������� �����, ��������� ��� ������� (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out)
void init_solution_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, char *files_directory, struct TimeMoment *time_mom, double **initial_solution ) {

    FILE *restart_info_file = NULL; // ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������
    int i_cell_global = 0; // ���������� ����� ������
    
    /* ������� ������� ���� � ����������� � ��������, ������������ �������� params->isContinue */
/*    try_to_restart( paramsc, files_directory, time_mom, restart_info_file );

    if ( paramsc->isContinue ) {
        // ����������� ������ �������� ������������ �����������
        read_solution( paramsc, params1d, restart_info_file, initial_solution );
    }
    else { */
        /* ������������� ����� � ��������� ������� */
        for ( int i_block = 0; i_block < params1d->ic_blocks_number; i_block++ ) { // ���� �� ������
            for ( int i_cell = params1d->cell_begin[i_block]; i_cell <= params1d->cell_end[i_block]; i_cell++ ) { // ���� �� ������� ������ �����
                initial_solution[i_cell_global][B_DISP] = params1d->block_values[i_block][B_DISP];
                initial_solution[i_cell_global][R_DISP] = params1d->block_values[i_block][R_DISP];
                initial_solution[i_cell_global][V_DISP] = params1d->block_values[i_block][V_DISP];
                initial_solution[i_cell_global][P_DISP] = params1d->block_values[i_block][P_DISP];
                initial_solution[i_cell_global][R_GAS] = params1d->block_values[i_block][R_GAS];
                initial_solution[i_cell_global][V_GAS] = params1d->block_values[i_block][V_GAS];
                initial_solution[i_cell_global][P_GAS] = params1d->block_values[i_block][P_GAS];
                for (int i_scalar = 0; i_scalar < params1d->number_of_scalars; i_scalar++)
                    initial_solution[i_cell_global][Z0 + i_scalar] = params1d->scalar_values[i_block][i_scalar];
                i_cell_global++;
            }
        }
 //   }
    
}


// ��������� ���������� �������
// params1d - ��������� � ����������� ���������� ������ (in)
// v_ncons[M] - ������ ����������� ���������� (in)
// boun_type - ��� ���������� ������� (in)
// boun_v[M] - ������ ����������� ���������� � ��������� ������ (out)
// curr_t - ������� ������ ������� (in)
// curr_inflow_parameters[M] - ������� �������� ���������� ��������(in)
void boundary( struct Parameters1d *params1d, struct Parameters2d *params2d, double v_ncons[M], int boun_type, double boun_v[M], double curr_t, double curr_inflow_parameters[M] ) {

    switch ( boun_type ) {
        case WALL:
            boun_v[B_DISP] = v_ncons[B_DISP]; // �������� ���� ���������� ����
	    boun_v[R_DISP] = v_ncons[R_DISP]; // ��������� ���������� ����
	    boun_v[V_DISP] = - v_ncons[V_DISP]; // �������� ���������� ����
            boun_v[P_DISP] = v_ncons[P_DISP]; // �������� ���������� ����
            boun_v[R_GAS] = v_ncons[R_GAS]; // ��������� ������� ����
            boun_v[V_GAS] = - v_ncons[V_GAS]; // �������� ������� ����
            boun_v[P_GAS] = v_ncons[P_GAS]; // �������� ������� ����
	    break;
        case FREE:
            boun_v[B_DISP] = v_ncons[B_DISP];   /* �������� ���� ���������� ���� */
	    boun_v[R_DISP] = v_ncons[R_DISP];   /* ��������� ���������� ���� */ 
	    boun_v[V_DISP] = v_ncons[V_DISP];   /* �������� ���������� ���� */
            boun_v[P_DISP] = v_ncons[P_DISP];   /* �������� ���������� ���� */
            boun_v[R_GAS] = v_ncons[R_GAS];   /* ��������� ������� ���� */
            boun_v[V_GAS] = v_ncons[V_GAS];   /* �������� ������� ���� */
            boun_v[P_GAS] = v_ncons[P_GAS];   /* �������� ������� ���� */
            break;
        case INFLOW:
            boun_v[B_DISP] = curr_inflow_parameters[B_DISP]; /* �������� ���� ���������� ���� */
	    boun_v[R_DISP] = curr_inflow_parameters[R_DISP]; /* ��������� ���������� ���� */ 
	    boun_v[V_DISP] = curr_inflow_parameters[V_DISP]; /* �������� ���������� ���� */
            boun_v[P_DISP] = curr_inflow_parameters[P_DISP]; /* �������� ���������� ���� */
            boun_v[R_GAS] = curr_inflow_parameters[R_GAS]; /* ��������� ������� ���� */
            boun_v[V_GAS] = curr_inflow_parameters[V_GAS]; /* �������� ������� ���� */
            boun_v[P_GAS] = curr_inflow_parameters[P_GAS]; /* �������� ������� ���� */
            break;
        case COMPLEX_INFLOW:
            if ( curr_t < params1d->change_inflow_time ) {
                boun_v[B_DISP] = params1d->inflow_values[B_DISP]; /* �������� ���� ���������� ���� */
	        boun_v[R_DISP] = params1d->inflow_values[R_DISP]; /* ��������� ���������� ���� */ 
	        boun_v[V_DISP] = params1d->inflow_values[V_DISP]; /* �������� ���������� ���� */
                boun_v[P_DISP] = params1d->inflow_values[P_DISP]; /* �������� ���������� ���� */
                boun_v[R_GAS] = params1d->inflow_values[R_GAS]; /* ��������� ������� ���� */
                boun_v[V_GAS] = params1d->inflow_values[V_GAS]; /* �������� ������� ���� */
                boun_v[P_GAS] = params1d->inflow_values[P_GAS]; /* �������� ������� ���� */
            }
            else {
                boun_v[B_DISP] = params1d->changed_inflow_values[B_DISP]; /* �������� ���� ���������� ���� */
	        boun_v[R_DISP] = params1d->changed_inflow_values[R_DISP]; /* ��������� ���������� ���� */ 
	        boun_v[V_DISP] = params1d->changed_inflow_values[V_DISP]; /* �������� ���������� ���� */
                boun_v[P_DISP] = params1d->changed_inflow_values[P_DISP]; /* �������� ���������� ���� */
                boun_v[R_GAS] = params1d->changed_inflow_values[R_GAS]; /* ��������� ������� ���� */
                boun_v[V_GAS] = params1d->changed_inflow_values[V_GAS]; /* �������� ������� ���� */
                boun_v[P_GAS] = params1d->changed_inflow_values[P_GAS]; /* �������� ������� ���� */
            }
            break;
        default:
            printf( "\nboundary -> wrong boundary condition.\n" );
            exit( EXIT_FAILURE );
            break;
    }
    for (int i = 0; i < params1d->number_of_scalars; i++){
        if ( boun_type == INFLOW )
            boun_v[Z0 + i] = curr_inflow_parameters[Z0 + i];
        else
            boun_v[Z0 + i] = v_ncons[Z0 + i];
    }
}

void current_inflow_parameters( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, int boun_direction, double *curr_inflow_params){

    if ( paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC ){
        
        for (int i = 0; i < M1D + params1d->number_of_scalars; i++)
            curr_inflow_params[i] = params1d->inflow_values[i];

    }
    else if ( paramsc->program_name == TWOD2PHC ){
        if (boun_direction == LEFT_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_left_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_left_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_left_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_left_bc[P_DISP_2D];
            curr_inflow_params[R_GAS] = params2d->inflow_values_left_bc[R_GAS_2D];
            curr_inflow_params[V_GAS] = params2d->inflow_values_left_bc[V_GAS_2D];
            curr_inflow_params[P_GAS] = params2d->inflow_values_left_bc[P_GAS_2D]; 
        }
        else if (boun_direction == RIGHT_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_right_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_right_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_right_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_right_bc[P_DISP_2D];
            curr_inflow_params[R_GAS] = params2d->inflow_values_right_bc[R_GAS_2D];
            curr_inflow_params[V_GAS] = params2d->inflow_values_right_bc[V_GAS_2D];
            curr_inflow_params[P_GAS] = params2d->inflow_values_right_bc[P_GAS_2D];    
        }
        else if (boun_direction == UP_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_up_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_up_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_up_bc[U_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_up_bc[P_DISP_2D];
            curr_inflow_params[R_GAS] = params2d->inflow_values_up_bc[R_GAS_2D];
            curr_inflow_params[V_GAS] = params2d->inflow_values_up_bc[U_GAS_2D];
            curr_inflow_params[P_GAS] = params2d->inflow_values_up_bc[P_GAS_2D];
        }
        else if (boun_direction == DOWN_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_down_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_down_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_down_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_down_bc[P_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_down_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_down_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_down_bc[P_DISP_2D]; 
        }
        else{
            
            printf_s("Wrong boun_direction number - %d. Only %d, %d and %d, %d are allowed.\n", boun_direction, LEFT_BOUNDARY, RIGHT_BOUNDARY, UP_BOUNDARY, DOWN_BOUNDARY);
            exit(EXIT_FAILURE);

        }
    }
    else{
        
        printf_s("Wrong program name number - %d. Only %d, %d and %d are allowed.\n", paramsc->program_name, ONED2PHC, ONED3PHC, TWOD2PHC);
        exit(EXIT_FAILURE);

    }

}

// ������ ������ ���� �� ������� � ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� �����������
// u_left - ������ ����������� ���������� � ������ ����� �� ��������������
// u_center - ������ ����������� ���������� � �������������� ������
// u_right - ������ ����������� ���������� � ������ ������ �� ��������������
// slopes_left - ������ �������� � ������ ����� �� ��������������
// slopes_center - ������ �������� � �������������� ������
// slopes_right - ������ �������� � ������ ������ �� ��������������
// dt - ��� �������������� �� �������
// h - ������ ������
// u_next - ������ ����������� ���������� � �������������� ������ �� ��������� ��������� ����
// n - �������� ������ ��������
// is_pressure_relaxation_now - true, ���� ����� ������� ���������� ������ �� ������� ����������� ����� ��������� ���������� ��������; false - �����
//                              ���� � ������ ������ ������ �� ����� ��������� ���������� ��������, �� true �� �������� �� ���������, ��� ��� ������ ����� �������������� � parameters.dat
// number_of_scalars - ����� ����������� �������� � ������ ������
// curr_time - ������� ������ �������
// configuration_pressure - ���������������� ��������
void calc_step_in_cell( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double u_left[M], double u_center[M],
                        double u_right[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                        double dt, double h, double u_next[M], int step_number, int n, bool is_pressure_relaxation_now, int number_of_scalars, double curr_time, double *configuration_pressure, double body_velocity, int i, int *status, double cont_left[M], double cont_right[M] ) {
    
    switch ( paramsc->numerical_method ) {
        case CIR1:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR_DISP scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_1( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case CIR2:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR_GAS scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_2( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case CIR3:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR3 scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_3( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case CIR4:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR4 scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_4( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case GODUNOV:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> Godunov scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            godunov( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, is_pressure_relaxation_now, curr_time, configuration_pressure  );
            break;
        case HLL:
            hll( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, is_pressure_relaxation_now, params1d->number_of_scalars, curr_time, configuration_pressure, body_velocity, i, status,  cont_left, cont_right );
            break;
        case RUSANOV:
            rusanov_1d( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, curr_time, configuration_pressure  );
            break;
        case HLLC:
            hllc_1d( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, is_pressure_relaxation_now, params1d->number_of_scalars, curr_time, configuration_pressure, body_velocity, i, status,  cont_left, cont_right );
            break;
        default:
            printf( "\ncalc_step_in_cell -> wrong scheme number\n" );
            exit( EXIT_FAILURE );
        }
    //printf("\n%lf %lf %lf \n", u_next[P_GAS]);
}


// ������������� ��������� ��� �������� � �������� ���������� ����������
// output_file_directory - ����������, � ������� ������ ���� �������� �����
// debug_info - ��������� � ���������� �����������
void init_debug_info( const char *output_file_directory, struct DebugInfo *debug_info ) {

    (*debug_info).current_cell = -1;
    (*debug_info).current_cell_x = -1.0;
    (*debug_info).neighbour_cell = -1;
    for ( int i_component = 0; i_component < M; i_component++ ) {
        (*debug_info).current_cell_vncons[i_component] = -1.0;
        (*debug_info).neighbour_cell_vncons[i_component] = -1.0;
    }
    
    for ( int i_relax_case = 0; i_relax_case < RELAX_CASES_NUM; i_relax_case++ )
        (*debug_info).relaxation_cases[i_relax_case] = 0;
    for ( int i_godunov_case = 0; i_godunov_case < GODUNOV_CASES_NUM; i_godunov_case++ )
        (*debug_info).godunov_cases[i_godunov_case] = 0;

    char string[MAX_STRING_SIZE];

    strcpy_s( string, output_file_directory );
    strcat_s( string, "\\relaxation.dat" );
    if ( ( fopen_s( &(*debug_info).relaxation_out, string, "wt" ) ) != 0 ) {
        printf( "\nread_Parameters1d -> can't open file %s for format writing\n\n", string );
        exit( EXIT_FAILURE );
    }

    strcpy_s( string, output_file_directory );
    strcat_s( string, "\\godunov.dat" );
    if ( ( fopen_s( &(*debug_info).godunov_out, string, "wt" ) ) != 0 ) {
        printf( "\nread_Parameters1d -> can't open file %s for format writing\n\n", string );
        exit( EXIT_FAILURE );
    }

}

// ������ ������������ ���������� ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ���������������� ���������� � ��������������� ������
// ���������� ������� ��������
double calc_p_i( const struct ParametersCommon* paramsc, const double v_ncons[M] ) {

    switch ( paramsc->media_model ) {
        case BAER_NUNZIATO:
            return v_ncons[P_GAS];
        case SAUREL_ABGRALL:
            return ( v_ncons[B_DISP] * v_ncons[P_DISP] + ( 1.0 - v_ncons[B_DISP] ) * v_ncons[P_GAS] );
        default:
            printf( "calc_p_i -> wrong value of params->media_model\n" );
            exit( EXIT_FAILURE );
    }

}

// ������ ����������� ��������� ���������, ����� ����� ������ ��� ������ Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ���������������� ���������� � ��������������� ������
// ���������� ������� ���������
double calc_r_i( const struct ParametersCommon* paramsc, const double v_ncons[M] ) {

    switch ( paramsc->media_model ) {
        case SAUREL_ABGRALL:
            return ( v_ncons[B_DISP] * v_ncons[R_DISP] + ( 1.0 - v_ncons[B_DISP] ) * v_ncons[R_GAS] );
        case BAER_NUNZIATO:
        default:
            printf( "calc_r_i -> wrong value of params->media_model\n" );
            exit( EXIT_FAILURE );
    }

}

// ������ ��������� ����� ��� ����� ��� �� ������� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// c1 - �������� ����� � ���������� ����
// c2 - �������� ����� � ������� ����
void calc_sound_velocity( const struct ParametersCommon* paramsc, const double v_ncons[M], double* c1, double* c2 ) {

    if (v_ncons[P_DISP] < 0.0 || v_ncons[P_GAS] < 0.0){
        if (v_ncons[P_DISP] < 0.0)
            printf("\ncalc_sound_velocity -> pressure of solid phase is negative\n");
            printf("\ncalc_sound_velocity -> pressure of solid phase is %lf\n", v_ncons[P_DISP]);
            printf("\ncalc_sound_velocity -> volume_fraction of solid phase is %lf\n", v_ncons[B_DISP]);
        if (v_ncons[P_GAS] < 0.0)
            printf("\ncalc_sound_velocity -> pressure of gas phase is negative\n");
        exit(EXIT_FAILURE);
    }

    *c1 = sqrt( paramsc->g1 * ( v_ncons[P_DISP] + paramsc->p01 ) / v_ncons[R_DISP] );
    *c2 = sqrt( paramsc->g2 * ( v_ncons[P_GAS] + paramsc->p02 ) / v_ncons[R_GAS] * 
        (1.0 + paramsc->b_virial * v_ncons[R_GAS] - paramsc->b_virial * v_ncons[R_GAS] * paramsc->b_virial * v_ncons[R_GAS] / paramsc->g2 / (1.0 + paramsc->b_virial * v_ncons[R_GAS]) ) );

}

// ������ ��������� ����� ��� �������� ���� �� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// phase - ������������� ����, ��� ������� �������������� �������� ����� - ������� ��� ����������
// ���������� �������� �����
double calc_sound_velocity_one_phase( const struct ParametersCommon* paramsc, const double v_ncons[M], const Phase phase ) {

    if ( phase == GAS_PHASE ){
        if ( v_ncons[P_GAS] < 0.0){
            printf("\ncalc_sound_velocity_one_phase -> pressure of gas phase is negative %lf \n", v_ncons[P_GAS]);
            exit(EXIT_FAILURE);
        }
        return sqrt( paramsc->g2 * ( v_ncons[P_GAS] + paramsc->p02 ) / v_ncons[R_GAS] * 
        (1.0 + paramsc->b_virial * v_ncons[R_GAS] - paramsc->b_virial * v_ncons[R_GAS] * paramsc->b_virial * v_ncons[R_GAS] / paramsc->g2 / (1.0 + paramsc->b_virial * v_ncons[R_GAS]) ) );
    }
    else{
        if ( v_ncons[P_DISP] < 0.0){
            printf("\ncalc_sound_velocity_one_phase -> pressure of solid phase is negative\n");
            exit(EXIT_FAILURE);
        }
        return sqrt( paramsc->g1 * ( v_ncons[P_DISP] + paramsc->p01 ) / v_ncons[R_DISP] );
    }
}

// ������ �������� ����� ��� ����� ���� �� ����������� ������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// phase - ������������� ����, ��� ������� �������������� �������� ����� - ������� ��� ����������
// ���������� �������� �����
double calc_sound_velocity_reduced( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase ) {

    if ( phase == GAS_PHASE )
        return sqrt( paramsc->g2 * v_ncons[P] / v_ncons[R] );
    else
        return sqrt( paramsc->g1 * ( v_ncons[P] + paramsc->p01 ) / v_ncons[R] );

}

// ������ ����������� ��������� ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ���������������� ���������� � ��������������� ������
// ���������� ������� ��������
double calc_u_i( const struct ParametersCommon* paramsc, const double v_ncons[M] ) {

    double sum1, sum2; 
    
    switch ( paramsc->media_model ) {
        case BAER_NUNZIATO:
            return v_ncons[V_DISP];
        case SAUREL_ABGRALL:
            sum1 = v_ncons[B_DISP] * v_ncons[R_DISP] * v_ncons[V_DISP] + ( 1.0 - v_ncons[B_DISP] ) * v_ncons[R_GAS] * v_ncons[V_GAS];
            sum2 = calc_r_i( paramsc, v_ncons );
            return ( sum1 / sum2 );
        default:
            printf( "calc_u_i -> wrong value of params->media_model\n" );
            exit( EXIT_FAILURE );
    }

}

// �������������� ������� "��������������" ���������� � ������ �����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// v_ncons - ������ ����������� ����������
// ����������: SUCCEESS            �������� � ����� ����� ������������
//             NEGATIVE_PRESSURE   �������� ���� �� � ����� �� ��� ������������
ReturnCodes convert_cons_to_noncons( const struct ParametersCommon* paramsc, const double v_cons[M], double v_ncons[M], int number_of_scalars ) {

    double b2 = 1 - v_cons[B_DISP]; // �������� ���� ������� ����
    ReturnCodes return_code = SUCCESS; // ��� �������� �������

    v_ncons[B_DISP] = v_cons[B_DISP]; // b1, �������� ���� ���������� ����
    v_ncons[R_DISP] = v_cons[R_DISP] / v_cons[B_DISP]; // r1, ��������� ���������� ����
    v_ncons[V_DISP] = v_cons[V_DISP] / v_cons[R_DISP]; // v1, �������� ���������� ����
    v_ncons[P_DISP] = ( paramsc->g1 - 1.0 ) * ( v_cons[P_DISP] / v_cons[B_DISP] - 0.5 * pow( v_cons[V_DISP], 2.0 ) / v_cons[B_DISP] / v_cons[R_DISP] ) -
        paramsc->g1 * paramsc->p01; // p1, �������� ���������� ����
    v_ncons[R_GAS] = v_cons[R_GAS] / b2; // r2, ��������� ������� ����
    v_ncons[V_GAS] = v_cons[V_GAS] / v_cons[R_GAS]; // v2, �������� ������� ����
    v_ncons[P_GAS] = ( paramsc->g2 - 1.0 ) * (1.0 + paramsc->b_virial * v_cons[R_GAS] / b2) * ( v_cons[P_GAS] / b2 - 0.5 * pow( v_cons[V_GAS], 2.0 ) / b2 / v_cons[R_GAS] ) -
        paramsc->g2 * paramsc->p02; // p2, �������� ������� ���� 

    for (int i = 0; i < number_of_scalars; i++)
        v_ncons[Z0 + i] = v_cons[Z0 + i] / v_cons[R_DISP];

    return return_code;

}

/* �������������� ��������������� ����������� ������� "��������������" ���������� � ������ �����������

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_cons[M_REDUCTION] - ������ "��������������" ���������� (in)
   phase - ������������� ����, ��� ������� �������������� �������������� - ������� ��� ���������� (in)

   v_ncons[M_REDUCTION] - ������ ����������� ���������� (out) */
void convert_cons_to_noncons_reduction( struct ParametersCommon *paramsc, double v_cons[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] ) {

    /* ��������� */
    v_ncons[R] = v_cons[R];
    
    /* �������� */
    v_ncons[V] = v_cons[V] / v_cons[R];
    
    /* �������� */
    switch ( phase ) {
        case GAS_PHASE:
            v_ncons[P] = ( paramsc->g2 - 1.0 ) * ( v_cons[P] - 0.5 * pow( v_cons[V], 2.0 ) / v_cons[R] );
            break;
        case DISPERSED_PHASE:
            v_ncons[P] = ( paramsc->g1 - 1.0 ) * ( v_cons[P] - 0.5 * pow( v_cons[V], 2.0 ) / v_cons[R] ) -
                paramsc->g1 * paramsc->p01;
            break;
        default:
            printf( "\nconvert_cons_to_noncons_reduction -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

/* ������������ ����������� ������� �� ������� ����������� �������
 
   full_vector[M] - ������ ���������� ������ (in)
   phase - ������������� ����, ��� ������� �� ������� ������� �������� �������������� - ������� ��� ���������� (in)
 
   reduced_vector[M_REDUCTION] - �������������� ���������� ������ (out) */
void convert_full_to_reduced( double full_vector[M], Phase phase, double reduced_vector[M_REDUCTION] ) {

    switch ( phase ) {
        case GAS_PHASE:
            reduced_vector[R] = full_vector[R_GAS];
            reduced_vector[V] = full_vector[V_GAS];
            reduced_vector[P] = full_vector[P_GAS];
            break;
        case DISPERSED_PHASE:
            reduced_vector[R] = full_vector[R_DISP];
            reduced_vector[V] = full_vector[V_DISP];
            reduced_vector[P] = full_vector[P_DISP];
            break;
        default:
            printf( "\nconvert_full_to_reduced -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

// �������������� ������� ����������� ���������� � ������ "��������������"
// params - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// v_cons - ������ "��������������" ����������
// number_of_scalars - ���������� ����������� ��������
void convert_noncons_to_cons( const struct ParametersCommon *paramsc, const double v_ncons[M], double v_cons[M], int number_of_scalars ) {

    double b2 = 1 - v_ncons[B_DISP]; // �������� ���� ������� ���� 

    v_cons[B_DISP] = v_ncons[B_DISP]; // b1, �������� ���� ���������� ���� 
    v_cons[R_DISP] = v_ncons[B_DISP] * v_ncons[R_DISP]; // b1 * r1
    v_cons[V_DISP] = v_cons[R_DISP] * v_ncons[V_DISP]; // b1 * r1 * v1
    v_cons[P_DISP] = v_ncons[B_DISP] * v_ncons[R_DISP] * ( 0.5 * pow( v_ncons[V_DISP], 2.0 ) + v_ncons[P_DISP] / v_ncons[R_DISP] /
        ( paramsc->g1 - 1.0 ) + paramsc->g1 * paramsc->p01 / ( paramsc->g1 - 1.0 ) / v_ncons[R_DISP] ); // b1 * r1 * E1 
    v_cons[R_GAS] = b2 * v_ncons[R_GAS]; // b2 * r2
    v_cons[V_GAS] = v_cons[R_GAS] * v_ncons[V_GAS]; // b2 * r2* v2
    v_cons[P_GAS] = v_cons[R_GAS] * ( 0.5 * pow( v_ncons[V_GAS], 2.0 ) + ( v_ncons[P_GAS] + paramsc->g2 * paramsc->p02 )/ v_ncons[R_GAS] /
		( paramsc->g2 - 1.0 ) / (1.0 + paramsc->b_virial * v_ncons[R_GAS])) ; // b2 * r2 * E2

    for (int i = 0; i < number_of_scalars; i++)
        v_cons[Z0 + i] = v_ncons[Z0 + i] * v_ncons[R_DISP] * v_ncons[B_DISP];
}


/* �������������� ����������� ������� ����������� ���������� � ������ "��������������"

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons[M_REDUCTION] - ������ ����������� ���������� (in)
   phase - ������������� ����, ��� ������� �������������� �������������� - ������� ��� ���������� (in)

   v_cons[M_REDUCTION] - ������ "��������������" ���������� (out) */
void convert_noncons_to_cons_reduction( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase, double v_cons[M_REDUCTION] ) {

    /* ����� */
    v_cons[R] = v_ncons[R];
    
    /* ������� */
    v_cons[V] = v_ncons[R] * v_ncons[V];

    /* ������� */
    switch ( phase ) {
        case GAS_PHASE:
            v_cons[P] = v_ncons[R] * ( 0.5 * pow( v_ncons[V], 2.0 ) + v_ncons[P] / v_ncons[R] / ( paramsc->g2 - 1.0 ) );
            break;
        case DISPERSED_PHASE:
            v_cons[P] = v_ncons[R] * ( 0.5 * pow( v_ncons[V], 2.0 ) + v_ncons[P] / v_ncons[R] /
                ( paramsc->g1 - 1.0 ) + paramsc->g1 * paramsc->p01 / ( paramsc->g1 - 1.0 ) / v_ncons[R] );
            break;
        default:
            printf( "\nconvert_noncons_to_cons_reduction -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

// ������������ ����� ������� ����������� ������� �� ������������ �����������
// reduced_vector[M_REDUCTION] - �������������� ���������� ������ (in)
// phase - ������������� ����, ��� ������� �� ��������������� ������� �������� ����� ������� - ������� ��� ���������� (in)
// full_vector[M] - ������ ���������� ������ (out)
void convert_reduced_to_full( double reduced_vector[M], Phase phase, double full_vector[M_REDUCTION] ) {

    switch ( phase ) {
        case GAS_PHASE:
            full_vector[R_GAS] = reduced_vector[R];
            full_vector[V_GAS] = reduced_vector[V];
            full_vector[P_GAS] = reduced_vector[P];
            break;
        case DISPERSED_PHASE:
            full_vector[R_DISP] = reduced_vector[R];
            full_vector[V_DISP] = reduced_vector[V];
            full_vector[P_DISP] = reduced_vector[P];
            break;
        default:
            printf( "\nconvert_reduced_to_full -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

// ������ ������� ����������������� "������" �� ������� "��������������" ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// flux - �������������� ������ ����������������� "������"
void diff_flux_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* flux, int number_of_scalars ) {

    double v_ncons[M]; // ������ ����������� ����������
        
    convert_cons_to_noncons( paramsc, v_cons, v_ncons, number_of_scalars );

    double b2 = 1 - v_ncons[B_DISP]; // �������� ���� ������� ����

    (*flux)[B_DISP] = 0.0; // ����� �������� ���� ���������� ����
    (*flux)[R_DISP] = v_cons[V_DISP]; // ����� ����� � ���������� ����
    (*flux)[V_DISP] = v_ncons[B_DISP] * ( v_ncons[R_DISP] * pow( v_ncons[V_DISP], 2.0 ) + v_ncons[P_DISP] ); // ����� �������� � ���������� ����
    (*flux)[P_DISP] = v_ncons[V_DISP] * ( v_cons[P_DISP] + v_ncons[B_DISP] * v_ncons[P_DISP] ); // ����� ������� � ���������� ����
    (*flux)[R_GAS] = v_cons[V_GAS]; // ����� ����� � ������� ����
    (*flux)[V_GAS] = b2 * ( v_ncons[R_GAS] * pow( v_ncons[V_GAS], 2.0 ) + v_ncons[P_GAS] ); // ����� �������� � ������� ����
    (*flux)[P_GAS] = v_ncons[V_GAS] * ( v_cons[P_GAS] + b2 * v_ncons[P_GAS] ); // ����� ������� � ������� ����

    for ( int i = 0; i < number_of_scalars; i++ )
        (*flux)[Z0 + i] = v_cons[V_DISP] * v_ncons[Z0 + i];
}

// ������ ������� ����������������� "������" �� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// flux - �������������� ������ ����������������� "������"
void diff_flux_ncons( const struct ParametersCommon* paramsc, const double v_ncons[M], array1D* flux, int number_of_scalars ) {

    double v_cons[M]; // ������ "��������������" ����������

    convert_noncons_to_cons( paramsc, v_ncons, v_cons, number_of_scalars );

    double b2 = 1 - v_ncons[B_DISP]; // �������� ���� ������� ����
    
    (*flux)[B_DISP] = 0.0; // ����� �������� ���� ���������� ����
    (*flux)[R_DISP] = v_cons[V_DISP]; // ����� ����� � ���������� ����
    (*flux)[V_DISP] = v_ncons[B_DISP] * ( v_ncons[R_DISP] * pow( v_ncons[V_DISP], 2.0 ) + v_ncons[P_DISP] ); // ����� �������� � ���������� ����
    (*flux)[P_DISP] = v_ncons[V_DISP] * ( v_cons[P_DISP] + v_ncons[B_DISP] * v_ncons[P_DISP] ); // ����� ������� � ���������� ����
    (*flux)[R_GAS] = v_cons[V_GAS]; // ����� ����� � ������� ����
    (*flux)[V_GAS] = b2 * ( v_ncons[R_GAS] * pow( v_ncons[V_GAS], 2.0 ) + v_ncons[P_GAS] ); // ����� �������� � ������� ����
    (*flux)[P_GAS] = v_ncons[V_GAS] * ( v_cons[P_GAS] + b2 * v_ncons[P_GAS] ); // ����� ������� � ������� ����

    for ( int i = 0; i < number_of_scalars; i++ )
        (*flux)[Z0 + i] = v_cons[V_DISP] * v_ncons[Z0 + i];
}

/* ������ ����������� ������� ����������������� "������" �� ������������ ����������� ������� ����������� ����������

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons - ����������� ���������� ������ ����������� ���������� (in)
   phase - ���� - ������� ��� ����������, ��� ������� �������������� "�����" (in)

   flux - �������������� ����������� ���������� ������ ����������������� "������" (out) */
void diff_flux_ncons_reduced( struct ParametersCommon *paramsc, double *v_ncons, Phase phase, double *flux ) {

    double v_cons[M_REDUCTION];     /* ����������� ���������� ������ "��������������" ���������� */
    
    /* �������������� ������� ���������������� ���������� � ������ "��������������" */
    convert_noncons_to_cons_reduction( paramsc, v_ncons, phase, v_cons );
    
    flux[R] = v_cons[V];                                        /* ����� ����� */
    flux[V] = v_ncons[R] * pow( v_ncons[V], 2.0 ) + v_ncons[P]; /* ����� �������� */
    flux[P] = v_ncons[V] * ( v_cons[P] + v_ncons[P] );          /* ����� ������� */

}

// ������ ����������������� ������� ������ ������ �� ������� "��������������" ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// rhst - ���������������� ������ ������ ������
void rhst_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* rhst, int number_of_scalars ) {

    double v_ncons[M]; // ������ ����������� ����������
        
    convert_cons_to_noncons( paramsc, v_cons, v_ncons, number_of_scalars );

    (*rhst)[B_DISP] = - v_ncons[V_DISP]; // ��������� ���������������
    (*rhst)[R_DISP] = 0.0; // ��� � ���������� ����
    (*rhst)[V_DISP] = v_ncons[P_GAS]; // ��� � ���������� ����
    (*rhst)[P_DISP] = v_ncons[P_GAS] * v_ncons[V_DISP]; // ��� � ���������� ����
    (*rhst)[R_GAS] = 0.0; // ��� � ������� ����
    (*rhst)[V_GAS] = - (*rhst)[V_DISP]; // ��� � ������� ����
    (*rhst)[P_GAS] = - (*rhst)[P_DISP]; // ��� � ������� ����

    for ( int i = 0; i < number_of_scalars; i++ )
        (*rhst)[Z0 + i] = 0.0;
}

// ������ ����������������� ������� ������ ������ �� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// rhst - ���������������� ������ ������ ������
void rhst_ncons( const struct ParametersCommon* params�, const double v_ncons[M], array1D* rhst, int number_of_scalars  ) {

    double p_i = calc_p_i( params�, v_ncons ); // ����������� ��������� ��������
    double u_i = calc_u_i( params�, v_ncons ); // ����������� ��������� ��������
    
    (*rhst)[B_DISP] = - u_i; // ��������� ���������������
    (*rhst)[R_DISP] = 0.0; // ��� � ���������� ����
    (*rhst)[V_DISP] = p_i; // ��� � ���������� ����
    (*rhst)[P_DISP] = p_i * u_i; // ��� � ���������� ����
    (*rhst)[R_GAS] = 0.0; // ��� � ������� ����
    (*rhst)[V_GAS] = - p_i; // ��� � ������� ����
    (*rhst)[P_GAS] = - (*rhst)[P_DISP]; // ��� � ������� ����

    for ( int i = 0; i < number_of_scalars; i++ )
        (*rhst)[Z0 + i] = 0.0;

}

// ���������� ������ ��������� ������� � ������
// debug_info - ��������� � ���������� �����������
void vector_debug_print( const struct DebugInfo *debug_info ) {

    printf( "Vector of primitive variables in cell:\n" );
    printf( "b1 = %.3e\n", debug_info->current_cell_vncons[B_DISP] );
    printf( "r1 = %.3e\n", debug_info->current_cell_vncons[R_DISP] );
    printf( "v1 = %.3e\n", debug_info->current_cell_vncons[V_DISP] );
    printf( "p1 = %.3e\n", debug_info->current_cell_vncons[P_DISP] );
    printf( "r2 = %.3e\n", debug_info->current_cell_vncons[R_GAS] );
    printf( "v2 = %.3e\n", debug_info->current_cell_vncons[V_GAS] );
    printf( "p2 = %.3e\n", debug_info->current_cell_vncons[P_GAS] );

}

// ������ ������ ���������� � ���������� ����� � ������ ��������� ��������� � ������ ��������
// debug_info - ��������� � ���������� ����������� (in)
// index - ������ ����������� �������� � ������� �������� (in)
void debug_print( struct DebugInfo *debug_info, int index ) {

    // �������� � ��������� � ����� ����?
    switch ( index ) {
        case GAS_LEFT:
        case GAS_RIGHT:
            // ����������, ���������� ��������� � ������� ����
            printf( "for the gas phase\n" );
            break;
        case DISP_LEFT:
        case DISP_RIGHT:
            // ����������, ���������� ��������� � ���������� ����
            printf( "for the dispersed phase\n" );
            break;
        }

    printf( "\nDebug report: current cell %d (x = %f), neighbour cell %d\n", debug_info->current_cell, debug_info->current_cell_x,
        debug_info->neighbour_cell );

    printf( "\nCurrent cell primitive vector:\n\n" );
    printf( "b1 = %.3e\n", debug_info->current_cell_vncons[B_DISP] );
    printf( "r1 = %.3e\n", debug_info->current_cell_vncons[R_DISP] );
    printf( "v1 = %.3e\n", debug_info->current_cell_vncons[V_DISP] );
    printf( "p1 = %.3e\n", debug_info->current_cell_vncons[P_DISP] );
    printf( "r2 = %.3e\n", debug_info->current_cell_vncons[R_GAS] );
    printf( "v2 = %.3e\n", debug_info->current_cell_vncons[V_GAS] );
    printf( "p2 = %.3e\n", debug_info->current_cell_vncons[P_GAS] );

    printf( "\nNeighbour cell primitive vector:\n\n" );
    printf( "b1 = %.3e\n", debug_info->neighbour_cell_vncons[B_DISP] );
    printf( "r1 = %.3e\n", debug_info->neighbour_cell_vncons[R_DISP] );
    printf( "v1 = %.3e\n", debug_info->neighbour_cell_vncons[V_DISP] );
    printf( "p1 = %.3e\n", debug_info->neighbour_cell_vncons[P_DISP] );
    printf( "r2 = %.3e\n", debug_info->neighbour_cell_vncons[R_GAS] );
    printf( "v2 = %.3e\n", debug_info->neighbour_cell_vncons[V_GAS] );
    printf( "p2 = %.3e\n\n", debug_info->neighbour_cell_vncons[P_GAS] );

}

// ���������� ���������� ������ � ������������� ����
// paramsc - ��������� � ��������� ����������� ��������������� ������������, ����� ����� �������� ��������� ���������� � ������������� ���� (in/out)
// params1d - ��������� � ����������� ���������� ������, ����� ����� �������� ��������� ���������� � ������������� ���� (in/out)
void dimensionalization( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d ) {
	
    double density_scale = paramsc->mass_scale / pow( paramsc->length_scale, 3.0 ); // ����������� ������� ���������
    double velocity_scale = paramsc->length_scale / paramsc->time_scale; // ����������� ������� ��������
    double pressure_scale = paramsc->mass_scale / paramsc->length_scale / pow( paramsc->time_scale, 2.0 ); // ����������� ������� ��������

    double energy_scale = pressure_scale / density_scale;

    // ��� �� �������
    if ( paramsc->constant_time_step )
        paramsc->dt /= paramsc->time_scale;

    // ���
    paramsc->p01 /= pressure_scale; // ��������� � ��������� ��������� ���������� ����
    paramsc->p02 /= pressure_scale; // ��������� � ��������� ��������� ������� ����
    // ����������� ������������ �������� � ����
    paramsc->mu2 /= ( paramsc->mass_scale / ( paramsc->length_scale * paramsc->time_scale ) );

    // ����� �����������
    paramsc->stop_time /= paramsc->time_scale;

    // ������� ��������� ��� ������� ������ ���������� ���������� ���� �� ���� �� ������ �� �������
    paramsc->background_density /= density_scale; // ��������� ���������� ����
    paramsc->background_velocity /= velocity_scale; // �������� ���������� ����
    paramsc->background_pressure /= pressure_scale; // �������� ���������� ����

    // "������"
    paramsc->particle_diameter /= paramsc->length_scale; // ������� ������� ���������� ����
    paramsc->compaction_viscosity /= (pressure_scale * paramsc->time_scale);
    paramsc->interface_drag_coef /= (density_scale / paramsc->time_scale);
	
    paramsc->tau_parameter /= (energy_scale / paramsc->mass_scale);

    if (paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC){

        // ������� ��������� �������
        params1d->left_boundary_x /= paramsc->length_scale;
        params1d->right_boundary_x /= paramsc->length_scale;

        // ��������� �������
        for ( int i_block = 0; i_block < params1d->ic_blocks_number; i_block++ ) {
            params1d->block_values[i_block][R_DISP] /= density_scale; // ��������� ���������� ����
            params1d->block_values[i_block][V_DISP] /= velocity_scale; // �������� ���������� ����
            params1d->block_values[i_block][P_DISP] /= pressure_scale; // �������� ���������� ����
            params1d->block_values[i_block][R_GAS] /= density_scale; // ��������� ������� ����
            params1d->block_values[i_block][V_GAS] /= velocity_scale; // �������� ������� ����
            params1d->block_values[i_block][P_GAS] /= pressure_scale; // �������� ������� ����
        }
    
        // ��������� �����
        if ( params1d->left_bc == INFLOW || params1d->right_bc == INFLOW || params1d->left_bc == COMPLEX_INFLOW || params1d->right_bc == COMPLEX_INFLOW ) {
            params1d->inflow_values[R_DISP] /= density_scale; // ��������� ���������� ����
            params1d->inflow_values[V_DISP] /= velocity_scale; // �������� ���������� ����
            params1d->inflow_values[P_DISP] /= pressure_scale; // �������� ���������� ����
            params1d->inflow_values[R_GAS] /= density_scale; // ��������� ������� ����
            params1d->inflow_values[V_GAS] /= velocity_scale; // �������� ������� ����
            params1d->inflow_values[P_GAS] /= pressure_scale; // �������� ������� ����
        }

        // ��������� �������� �����
        if ( params1d->left_bc == COMPLEX_INFLOW || params1d->right_bc == COMPLEX_INFLOW ) {
            params1d->change_inflow_time /= paramsc->time_scale;
            params1d->changed_inflow_values[R_DISP] /= density_scale; // ��������� ���������� ����
            params1d->changed_inflow_values[V_DISP] /= velocity_scale; // �������� ���������� ����
            params1d->changed_inflow_values[P_DISP] /= pressure_scale; // �������� ���������� ����
            params1d->changed_inflow_values[R_GAS] /= density_scale; // ��������� ������� ����
            params1d->changed_inflow_values[V_GAS] /= velocity_scale; // �������� ������� ����
            params1d->changed_inflow_values[P_GAS] /= pressure_scale; // �������� ������� ����
        }
    }
    else if (paramsc->program_name == TWOD2PHC){

        double velocity_scale1 = paramsc->length_scale / paramsc->time_scale;    /* ����������� ������� �������� */
        double velocity_scale2 = paramsc->length_scale / paramsc->time_scale;   // (S) ��� 2-�� ������� ��� �������� ������ ��� �������� (������� �������� � �� �����)

        /* ������� ��������� ������� */
        params2d->left_boundary_x /= paramsc->length_scale;
        params2d->right_boundary_x /= paramsc->length_scale;

        /* ��������� ������� */
        for ( int i_block = 0; i_block < params2d->x_blocks_number; i_block++ ) 
        {
            for ( int j_block = 0; j_block < params2d->y_blocks_number; j_block++ ) 
            {
                params2d->block_values[i_block][j_block][R_DISP_2D] /= density_scale;     /* ��������� ���������� ���� */
                params2d->block_values[i_block][j_block][V_DISP_2D] /= velocity_scale1;    /* �������� ���������� ���� */
                params2d->block_values[i_block][j_block][U_DISP_2D] /= velocity_scale2;    /* �������� ���������� ���� */
                params2d->block_values[i_block][j_block][P_DISP_2D] /= pressure_scale;    /* �������� ���������� ���� */
                params2d->block_values[i_block][j_block][R_GAS_2D] /= density_scale;     /* ��������� ������� ���� */
                params2d->block_values[i_block][j_block][V_GAS_2D] /= velocity_scale1;    /* �������� ������� ���� */
                params2d->block_values[i_block][j_block][U_GAS_2D] /= velocity_scale2;    /* �������� ������� ���� */
                params2d->block_values[i_block][j_block][P_GAS_2D] /= pressure_scale;    /* �������� ������� ���� */
            }
        }

        if ( params2d->left_bc == INFLOW ) 
        {
            params2d->inflow_values_left_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_left_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->right_bc == INFLOW ) 
        {
            params2d->inflow_values_right_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_right_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->down_bc == INFLOW ) 
        {
            params2d->inflow_values_down_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_down_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_down_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_down_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_down_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_down_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_down_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_down_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->up_bc == INFLOW ) 
        {
            params2d->inflow_values_up_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_up_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_up_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_up_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_up_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_up_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_up_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_up_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->left_bc == COMPLEX_INFLOW ) 
        {
            params2d->inflow_values_left_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_left_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->right_bc == COMPLEX_INFLOW ) 
        {
            params2d->inflow_values_right_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_right_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_GAS_2D] /= pressure_scale;    
        }

    }
}

void current_block_number(struct ParametersCommon *params, struct Parameters1d *params1d, struct Parameters2d *params2d, int step_number_x, int step_number_y, int *number_of_block){

    number_of_block[0] = 0;

    if (params->program_name == ONED2PHC || params->program_name == ONED3PHC){
        for (int i = 0; i < params1d->ic_blocks_number; i++){
	    if (step_number_x <= params1d->cell_end[i] && step_number_x >= params1d->cell_begin[i])
	        number_of_block[0] = i;
        }
    }
    else if (params->program_name == TWOD2PHC){
        number_of_block[1] = 0;
        for (int i = 0; i < params2d->x_blocks_number; i++){
            for (int j = 0; j < params2d->y_blocks_number; j++){
                if (step_number_x <= params2d->x_end[i][j] && step_number_x >= params2d->x_begin[i][j] && step_number_y <= params2d->y_end[i][j] && step_number_y >= params2d->y_begin[i][j]){
                    number_of_block[0] = i;
                    number_of_block[1] = j;
                }
            }
        }
    }
    else{
        printf("Wrong program name is used in function current_block_number_2d\n");
        exit(EXIT_FAILURE);
    }
   
}

