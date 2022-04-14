// exact_solution_1d2phc.cc
// ���������� ������� ������� ������ � ������� ������� ��� ������� ��������� ����� - ��������.
// (c) ����� �����, 2013
// ������: 24 ������� 2013 �.

#include "exact_solution_1d2phc.h"

// ���������� ������� ������� ������ � ������� ������� ��� ������� ��������� �����-��������
// output_directory - ����������, ���� ����� ������� ���� � ������ �������� (in) 
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ (in)
// debug_info - ��������� � ���������� ����������� (in)
void build_exact_sol( char *output_directory, struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, int n ) {

    char output_filename[MAX_STRING_SIZE]; // ��� ��������� ����� �����
    double *xc; // ������ ��������� ������� ����� ����� ��� ������� �������
    double *x; // ������ ��������� ������ ����� ����� ��� ������� �������
    double s; // ������� �������� ������������� ����������
    FILE *ex_sol_out;
    double flux[M]; // ������ ������ - ���������� �� ������������, �� ����� ��� ������������� ������� godunov_cons_flux
    
    double p_solid_cont_l, p_solid_cont_r; // �������� � ���������� ���� ����� � ������ �� ����������� ������� � ���������� ����
    double v_solid_cont; // �������� ����������� ������� � ���������� ����
   
    double v_ncons_res[M]; // ������-������� ������ � ������� �������

    double time_moment; // ������ �������, ��� �������� �������� ������ �������

    Disp_phase_cases solver_part; // ���������, ������� ����������, ����� �� ������� ������ ������� ������������ ��� ������� "������"
    
    // ����������� ������� �������, ��� �������� �������� ������ �������, ������ �� ����� ������
    if ( paramsc->is_output_on_time )
        time_moment = paramsc->stop_time;
    else
        time_moment = paramsc->stop_steps_number * paramsc->dt;

    // ������������ ����� ��������� ����� � ������ ��������
    strcpy_s( output_filename, output_directory );
    strcat_s( output_filename, "\\exact_solution.dat" );
    fopen_s( &ex_sol_out, output_filename, "wt" );
    if ( NULL == ex_sol_out ) {
        printf( "build_exact_sol -> Can't open file %s for writing.\n", output_filename );
    }

    // ��������� ������ ��� ����� ��� ������� �������
    get_memory_for_1D_double_array( params1d->cells_number_for_exact_solution, &xc );
    get_memory_for_1D_double_array( params1d->cells_number_for_exact_solution + 1, &x );

    // ���������� ��������� �����
    build_grid( LEFT_BOUN_EX_SOL, RIGHT_BOUN_EX_SOL, params1d->cells_number_for_exact_solution, xc, x );

    // ������ ����������� ����� ���������� �������� ���� ���������� ���� ����� � ������ �� �������
    solver_part = what_case( paramsc, (params1d->block_values[0])[B_DISP], (params1d->block_values[1])[B_DISP] );
    double gas_left_ncons_reduced[M_REDUCTION]; // ����������� ������ ��� ���������� ������� ���� ����� �� �������
    double gas_right_ncons_reduced[M_REDUCTION]; // ����������� ������ ��� ���������� ������� ���� ������ �� �������
    double disp_left_ncons_reduced[M_REDUCTION]; // ����������� ������ ��� ���������� ���������� ���� ����� �� �������
    double disp_right_ncons_reduced[M_REDUCTION]; // ����������� ������ ��� ���������� ���������� ���� ������ �� �������
    if ( solver_part == NO_GRAD ) {
        // ��������� ����������� �������
        convert_full_to_reduced( params1d->block_values[0], GAS_PHASE, gas_left_ncons_reduced );
        convert_full_to_reduced( params1d->block_values[1], GAS_PHASE, gas_right_ncons_reduced );
        convert_full_to_reduced( params1d->block_values[0], DISPERSED_PHASE, disp_left_ncons_reduced );
        convert_full_to_reduced( params1d->block_values[1], DISPERSED_PHASE, disp_right_ncons_reduced );
    }

    // ���� �� �������
    for ( int i_cell = 0; i_cell < params1d->cells_number_for_exact_solution; i_cell++ ) {
        s = xc[i_cell] / time_moment; // ������� �������� ������������� ����������
        // ������ ���������� � ����
        fprintf( ex_sol_out, "%e ", xc[i_cell] + ( params1d->left_boundary_x - LEFT_BOUN_EX_SOL ) ); // ����������� � ������������ � ������� ��������� ��������
        if ( solver_part == NO_GRAD ) { // ������ ������� ����������� ���
                        
            // ������ ������� ��� ���������� ����
            double cl = calc_sound_velocity_reduced( paramsc, disp_left_ncons_reduced, DISPERSED_PHASE );
            double cr = calc_sound_velocity_reduced( paramsc, disp_right_ncons_reduced, DISPERSED_PHASE );
            // ������������ ��������� ������� �������� � �������� ���������� ���� �� ���������� �������
            double p_cont; // �������� �� ��������� �������
            double v_cont; // �������� �� ���������� �������
            calc_contact_pressure_velocity( paramsc, debug_info, disp_left_ncons_reduced, disp_right_ncons_reduced, M_REDUCTION,
                cl, cr, DISPERSED_PHASE, &p_cont, &v_cont );
            // ����� �������
            double v_ncons_disp_reduced[M_REDUCTION]; // ����������� ������
            sample_reduced( paramsc, disp_left_ncons_reduced, disp_right_ncons_reduced, M_REDUCTION, cl, cr, DISPERSED_PHASE,
                p_cont, v_cont, s, v_ncons_disp_reduced );
            
            // ������ ������� ��� ������� ����
            cl = calc_sound_velocity_reduced( paramsc, gas_left_ncons_reduced, GAS_PHASE );
            cr = calc_sound_velocity_reduced( paramsc, gas_right_ncons_reduced, GAS_PHASE );
            // ������������ ��������� ������� �������� � �������� ���� �� ���������� �������
            calc_contact_pressure_velocity( paramsc, debug_info, gas_left_ncons_reduced, gas_right_ncons_reduced, M_REDUCTION,
                cl, cr, GAS_PHASE, &p_cont, &v_cont );
            // ����� �������
            double v_ncons_gas_reduced[M_REDUCTION]; // ����������� ������
            sample_reduced( paramsc, gas_left_ncons_reduced, gas_right_ncons_reduced, M_REDUCTION, cl, cr, GAS_PHASE,
                p_cont, v_cont, s, v_ncons_gas_reduced );

            // ��������� ������ ������-�������
            v_ncons_res[B_DISP] = (params1d->block_values[0])[B_DISP];
            convert_reduced_to_full( v_ncons_disp_reduced, DISPERSED_PHASE, v_ncons_res );
            convert_reduced_to_full( v_ncons_gas_reduced, GAS_PHASE, v_ncons_res );

        }
        else { // ���� ������ �������� ���� ���������� ����
            ReturnCodes code = godunov_cons_flux( paramsc, debug_info, params1d->block_values[0], params1d->block_values[1], s, solver_part,
                flux, v_ncons_res, &v_solid_cont, &p_solid_cont_l, &p_solid_cont_r, n );
        }
        // ���������� ������ ������-�������
        for ( int i_component = 0; i_component < n; i_component++ )
            fprintf( ex_sol_out, "%e ", v_ncons_res[i_component] );
        fprintf( ex_sol_out, "\n" );
    }

    // ����������� ������
    free( xc );
    free( x );

    fclose( ex_sol_out );

}