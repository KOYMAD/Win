// main.cc
// ���������� (1D) ���������� (2PH) ��������� (C) ���
// (c) ����� �����, 2013 - 2018
// ������: 16 ������� 2013 �.

#include "main_1d2phc.h"

// ����� ��� ���� ������� �������
#include "grid.h"
#include "utils.h"
#include "utils_1d2phc.h"
#include "io_1d.h"
#include "exact_solution_1d2phc.h"
#include "memory.h"

// "������"
#include "physics_solver.h"
#include "source_terms.h"

// ��������� ������� �������������
#include "minmod.h"

int main( int argc, char *argv[] ) {

    double initial_total_mass; // �������� ��������� ����� �������� � �������
    double mass_diff; // ������������� ����������� ��������� ����� �������� � ������� ������������ ��������
    double time_gap; // ����� ����� �������� ������ � �������������� ������������, ����� ����� ��� params.is_output_on_time = 1
    int steps_gap; // ���������� ����� �� ������� ����� �������� ������ � �������������� ������������, ����� ����� ��� params.is_output_on_time = 0
    int files_counter = 0; // ������� ������ � �������������� ������������
    char output_filename[MAX_STRING_SIZE]; // ��� �������� ����� � �������������� ������������
    char tmp_str[MAX_STRING_SIZE]; // ��������� ���������� ��� �������� �������� ������

    printf( "\n1D two-phase compressible solver\n(c) Pavel Utkin, ICAD RAS, MIPT, 2013-2018\ne-mail: pavel_utk@mail.ru\n" );

    printf( "\nPreparing:\n" );

    // ��������� ���������� ��������� ������
    if ( argc != 3 ) {
        // argv[1] - ���� � ����������, ��� ��������� ���� � ����������� ������
        // argv[2] - ���� � ����������, ��� ����� �������� �������� �����
        printf( "\nmain -> wrong command line arguments number\n" );
        exit( EXIT_FAILURE );
    }

    // ���������� ����� � ����������� ������ � ���������� ��������� params
    int file_num = 0;
    struct ParametersCommon paramsc; // ��������� � ��������� ����������� ��������������� ������������
    struct Parameters1d params1d; // ��������� � ����������� ���������� ������
    struct Parameters2d empty_structure; // ��������� � ����������� 2d ������, �� ����������� � ���� ������
    fill_parameters_struct_1d( &paramsc, &params1d, argv[1], argv[2], file_num );
    printf( "\n> File Parameters1d.dat is processed successfully. Correspondent structure is filled.\n" );

    // ��������� ������ ��� �������
    double *xc; // ������ ��������� ������� ����� �����
    get_memory_for_1D_double_array( params1d.cells_number, &xc );
    double *x; // ������ ��������� ����� �����
    get_memory_for_1D_double_array( params1d.cells_number + 1, &x );
    double **u_prev; // ������� ����������� ���������� �� n-�� ����
    get_memory_for_2D_double_array( params1d.cells_number, M, &u_prev );
    double **u_next; // ������� ����������� ���������� �� (n+1)-�� ����
    get_memory_for_2D_double_array( params1d.cells_number, M, &u_next );
    double **slopes; // ������� ���� ��������� �������� �������������� ���������� �� ���� ������� ��������� ������� ��� ��������� ������� �������������
    get_memory_for_2D_double_array( params1d.cells_number, M, &slopes );
    double *beta; // configuration pressure
    get_memory_for_1D_double_array( params1d.cells_number, &beta );
    double *B; // ������� ��������������
    get_memory_for_1D_double_array( params1d.cells_number, &B );
    double *C; // ������� �������
    get_memory_for_1D_double_array( params1d.cells_number, &C );
    printf( "\n> All the necessary memory is allocated successfully.\n" );
    int *number_of_block;
    get_memory_for_1D_int_array( 1, &number_of_block );
    int *ignition_flag; // ��������� �������� ������� � ������
    get_memory_for_1D_int_array( params1d.cells_number, &ignition_flag );
    int *initial_ignition; // ��������� ���������� �������
    get_memory_for_1D_int_array( params1d.cells_number, &initial_ignition );
    double *S; // ������� ������
    get_memory_for_1D_double_array( params1d.cells_number, &S );
    // ���������� ��������� ���������� ������ � ������������� ����
    if ( paramsc.use_dimensions ) {
        dimensionalization( &paramsc, &params1d, &empty_structure );
        printf( "\n> Dimensionalization is done.\n" );
    }

    // ����������� ��������� ������� ����� �����
    build_grid( params1d.left_boundary_x, params1d.right_boundary_x, params1d.cells_number, xc, x );
    printf( "\n> Computational grid is prepared successfully.\n" );

    // ������������� �������-�������
    struct TimeMoment time_mom = { 0, 0.0 }; // ���������, ������������ ������� ������ �������
    init_solution_1d( &paramsc, &params1d, argv[1], &time_mom, u_prev );
    printf( "\n> The solution is initiated successfully.\n" );

    for (int i = 0; i < params1d.cells_number; i++){
	// ����������� �����, � ������� ������ ���������
	number_of_block[0] = 0;
	for (int j = 0; j < params1d.ic_blocks_number; j++){
	    if (i <= params1d.cell_end[j] && i >= params1d.cell_begin[j])
                number_of_block[0] = j;
	}
	ignition_flag[i] = 0;
	initial_ignition[i] = params1d.initial_burning[number_of_block[0]];
        beta[i] = calc_configuration_pressure(&paramsc, &params1d, &empty_structure, X_DIRECTION, u_prev[i], number_of_block);
	B[i] = compaction_energy(&paramsc, &params1d, &empty_structure, X_DIRECTION, u_prev[i], number_of_block);
        C[i] = calc_chemical_reaction(&paramsc, u_prev[i],ignition_flag[i]);
        S[i] = params1d.S[number_of_block[0]];
    }
    
    double curr_inflow_params[M]; // ��������� �������� � ������ ������ � ������ ������

    // ������ ����� � ���������� ���������������
    write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, 0, 0, output_filename, beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );

    // ���������������� ��������� ��� ������ � ��������� � �������� ������
    int sensors_cells[MAX_SENSORS_NUM]; // ��� ������� ������� �������� ����� ������, � ������� �� ���������
    FILE *sensors_files[MAX_SENSORS_NUM]; // ������ �������� ������������ ��� ������ ������ � ��������
    if ( paramsc.are_sensors ) {
        prepare_sensors( &paramsc, &params1d, x, sensors_cells, sensors_files );
    }

    // ������������� ����� ��������� debug_info
    struct DebugInfo debug_info; // ��������� � ���������� �����������
    init_debug_info( argv[2], &debug_info );

    // ���������� ������� �������
    if ( params1d.build_exact_solution ) {   
        build_exact_sol( argv[2], &paramsc, &params1d, &debug_info, M1D );
        printf( "\n> The exact solution is constructed successfully.\n" );
    }

    printf( "\nCalculation:\n" );

    // ������ ��������� �������� ��������� ����� ��������
    initial_total_mass = get_total_mass( &params1d, x, u_prev );

    // ���������� ������� ��� ���������� ����� ����� �������� ������ � �������������� ������������
    if ( paramsc.is_output_on_time )
        time_gap = paramsc.stop_time / paramsc.output_number;
    else
        steps_gap = paramsc.stop_steps_number / paramsc.output_number;

    // �������� �������  ����� ������� �������
    for ( int i = 0; i < params1d.cells_number; i++ ){
        for ( int j = 0; j < M1D + params1d.number_of_scalars; j++ )
            slopes[i][j] = 0.0;
	}
    double zero_slopes[M]; // ������ ������� �������� ��� ������� �� ��������
    for ( int i = 0; i < M; i++ )
        zero_slopes[i] = 0.0;

    double boun_v[M]; // ������ ����������� ���������� ��� ���������� ���������� �������
    double dt; // ��� �� �������

    // �������� ���� �� �������
    while ( ( ( paramsc.stop_steps_number - time_mom.steps_num > 0 ) && ( !paramsc.is_output_on_time ) ) ||
          ( ( paramsc.stop_time - time_mom.curr_t > 0 ) && ( paramsc.is_output_on_time ) ) ) {
	
        // ���������� ����� ���������� ��������� ��� ������� ��������� ����������
        for ( int i_relax_case = 0; i_relax_case < RELAX_CASES_NUM; i_relax_case++ )
            debug_info.relaxation_cases[i_relax_case] = 0;
        // ���������� ����� ���������� ��������� ��� ������� ������ ������ ��������
        for ( int i_godunov_case = 0; i_godunov_case < GODUNOV_CASES_NUM; i_godunov_case++ )
            debug_info.godunov_cases[i_godunov_case] = 0;

        // ������ �������� ���� �� �������
        if ( time_mom.steps_num >= 188 ){
            dt = get_time_step(&debug_info, &paramsc, &params1d, x, u_prev, 1 );
        }
        else {
            dt = get_time_step(&debug_info, &paramsc, &params1d, x, u_prev, 0 );
        }

        // ������ �������� �������� �������
        // ���������� ������������� �������� ������� ���������� � �������������, ��� ����� �����������
        if ( paramsc.approximation_order == 2 )
            reconstruction( &paramsc, &params1d, u_prev, xc, slopes, params1d.number_of_scalars );

        // ���� �� �������
        for ( int i = 0; i < params1d.cells_number; i++ ) {
            //printf("\n%lf", u_next[i][P_GAS]);

            // ���������� ����� ����� ��������� ��� ������� ����, ���������� ������ � ������� ������
           /* debug_info.current_cell = i;
            debug_info.current_cell_x = xc[i];
            for ( int j = 0; j < M1D + params1d.number_of_scalars; j++ )
                debug_info.current_cell_vncons[j] = u_prev[i][j];*/

            // ������ ���������� � i-�� ������ �� ��������� ��������� ����
            if ( i == 0 ) {
                // ��������� ������ ���������� �������
                current_inflow_parameters(&paramsc, &params1d, &empty_structure, LEFT_BOUNDARY, curr_inflow_params); 
                boundary( &params1d, &empty_structure, u_prev[0], params1d.left_bc, boun_v, time_mom.curr_t, curr_inflow_params );
                calc_step_in_cell( &paramsc, &params1d, &debug_info, boun_v, u_prev[0], u_prev[1], zero_slopes, zero_slopes, slopes[1], dt, x[1]-x[0], u_next[0], i, M1D + params1d.number_of_scalars, true, params1d.number_of_scalars, time_mom.curr_t, &beta[i], S[i], S[i], S[i+1], i );
            }
            else if ( i == params1d.cells_number - 1 ) {
                // ��������� ������� ���������� �������
                current_inflow_parameters(&paramsc, &params1d, &empty_structure, RIGHT_BOUNDARY, curr_inflow_params); 
                boundary( &params1d, &empty_structure, u_prev[params1d.cells_number-1], params1d.right_bc, boun_v, time_mom.curr_t, curr_inflow_params );

                calc_step_in_cell( &paramsc, &params1d,&debug_info, u_prev[params1d.cells_number-2], u_prev[params1d.cells_number-1], boun_v,
                    slopes[params1d.cells_number-2], zero_slopes, zero_slopes, dt, x[params1d.cells_number] - x[params1d.cells_number-1], u_next[params1d.cells_number-1], i, M1D + params1d.number_of_scalars, true, params1d.number_of_scalars, time_mom.curr_t, &beta[i], S[i], S[i-1], S[i], i );

            }
            else {
                // ������ ���������� �� ���������� ������
                calc_step_in_cell( &paramsc, &params1d, &debug_info, u_prev[i-1], u_prev[i], u_prev[i+1], slopes[i-1], slopes[i], slopes[i+1], dt, x[i+1] - x[i], u_next[i], i, M1D + params1d.number_of_scalars, true, params1d.number_of_scalars, time_mom.curr_t, &beta[i], S[i] , S[i-1], S[i+1], i);
            }


            if ( paramsc.is_physics ) {
                // ���� "������" ������	

		if( initial_ignition[i] == 1 && time_mom.curr_t < params1d.ignition_time){ // ������� �������
		    initial_ignition[i] = 1;
		}
		else{
			if (initial_ignition[i] == 2)
			{
			}
			else
				initial_ignition[i] = 0;
		}


                physics_solver( &paramsc, &params1d, &empty_structure, X_DIRECTION, dt, u_next[i], u_next[i+1], u_next[i-1], number_of_block, &beta[i], &B[i], &C[i], M1D + params1d.number_of_scalars, ignition_flag[i], initial_ignition[i], S[i], time_mom.curr_t );
            }

        } // ����� ����� �� �������

        // ������ ���������� � ����� ���������� �� ������ ����
        write_statistics( &paramsc, &params1d, &time_mom, &debug_info, u_prev );

        // ������ � ���� ���������� � ���������� � ������ 
        write_time_dependent_information_in_a_cell(&paramsc, &params1d, time_mom.curr_t, argv[2], "last_cell_infromation.dat", u_prev[0], M1D + params1d.number_of_scalars);

        for ( int i = 0; i < params1d.cells_number; i++ )
            for ( int j = 0; j < M1D + params1d.number_of_scalars; j++ )
                u_prev[i][j] = u_next[i][j];

        if (time_mom.steps_num == 188){
            printf("time_mom.steps_num == 188");
            write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, files_counter + 1 ,time_mom.curr_t, output_filename,  beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );
        }

        
        time_mom.curr_t += dt;
        time_mom.steps_num++;

        // �������� ���������� ����� � �������
        check_mass( &params1d, x, u_next, initial_total_mass, &mass_diff );

        // ������ ��������������� ��������� � ����� ����
        printf( "\nStep %d, time %f: dt = %.2e\n", time_mom.steps_num - 1, time_mom.curr_t - dt, dt );
        printf( "Pressure relaxation info - case 1: %d  case 2: %d  case 3: %d  case 4: %d  case 5: %d\n", debug_info.relaxation_cases[0], debug_info.relaxation_cases[1],
        debug_info.relaxation_cases[2], debug_info.relaxation_cases[3], debug_info.relaxation_cases[4] );
        printf( "Godunov algorithm info - case 1: %d  case 2: %d  case 3: %d\n", debug_info.godunov_cases[0], debug_info.godunov_cases[1], debug_info.godunov_cases[2] );

        // ������ ������������� �����������
        if ( paramsc.is_output_on_time  ) {
            // ����� �� �������� �������
            if ( time_mom.curr_t >= ( files_counter + 1 ) * time_gap ) {
                write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, files_counter + 1,time_mom.curr_t, output_filename, beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );
		// ������ ��������� ��������
		if ( paramsc.are_sensors ) {
		    write_sensors( &paramsc, u_prev, &time_mom, sensors_files, sensors_cells, M1D + params1d.number_of_scalars );
		}

                files_counter++;
            }
            sprintf_s( tmp_str, "%f", ( files_counter + 1 ) * time_gap );
            write_restart_info( argv[2], tmp_str, output_filename );
        }
        else {
            // ����� �� ���������� ����� �� �������
            if ( time_mom.steps_num >= ( files_counter + 1 ) * steps_gap ) {
                write_solution_1d2phc( &paramsc, &params1d, argv[2], params1d.cells_number, xc, u_prev, files_counter + 1 ,time_mom.curr_t, output_filename,  beta, B, C, M1D + params1d.number_of_scalars, ignition_flag, initial_ignition );
                files_counter++;
            }
            sprintf_s( tmp_str, "%d", ( files_counter + 1 ) * steps_gap );
            write_restart_info( argv[2], tmp_str, output_filename );
        }


    } // ����� ����� �� �������

    // ������������ ������
    free( xc );
    free( x );
    for ( int i = 0; i < params1d.cells_number; i++ ) {
        free( u_prev[i] );
        free( u_next[i] );
    }
    for ( int i = 0; i < params1d.cells_number; i++ ) {
        free( slopes[i] );
    }
	free( B );
	free( beta );
	free( ignition_flag );
	free( initial_ignition);
        free( S );
    // �������� ������ � ����������� ��������
    if ( paramsc.are_sensors ) {
        for ( int i_sensor = 0; i_sensor < paramsc.sensors_num; i_sensor++ )
            fclose( sensors_files[i_sensor] );
    }
    // �������� ������ ��� ������ ����������
    fclose( debug_info.relaxation_out );
    fclose( debug_info.godunov_out );

    printf( "\nNormal finish: total time steps = %d, total time = %f\n\n", time_mom.steps_num, time_mom.curr_t );
    
    return 0;

}