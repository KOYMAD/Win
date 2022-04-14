// io_1d.cc
// ������� �����/������
// (c) ����� �����, 2013 - 2018
// ������: 17 ��� 2012 �.

#include "io.h"
#include "io_1d.h"
#include "eos.h"

// ���������� ��������� � ����������� ���������� ������
// params - ���������� ����������������� ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
void fill_parameters_1d( FILE *params, struct ParametersCommon* paramsc, struct Parameters1d* params1d ) {

    int dtmp; // ��� ���������� ������������� ����������
    char string[MAX_STRING_SIZE]; // ��� ���������� ��������� ���������� �� �����

    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ���������, ������������, ��� ����������� ����� ����������, ����������� ���������� ������

    // ��������� �����
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->left_boundary_x) ); // ���������� ������ ����� ��������� �������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->right_boundary_x) ); // ���������� ������� ����� ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->cells_number) ); // ���������� �����
    if ( params1d->cells_number < 0 ) {
        printf( "\nfill_parameters_struct -> cells_number should be a positive value\n\n" );
        exit( EXIT_FAILURE );
    }

    // ��������� �������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->number_of_scalars) ); // ���������� ����������� ��������
    if ( params1d->number_of_scalars < 0 ) {
        printf( "\nfill_parameters_struct -> number_of_scalars should be a non-negative value\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->ignition_time) ); // ����� �������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->additional_heat) ); // ������� �������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->mass_input) ); // �������� ��� ������� �����
    // ������������� ����� � ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->ic_blocks_number) ); // ���������� ������
    if ( params1d->ic_blocks_number <= 0 ) {
        printf( "\fill_parameters_1d -> ic_blocks_number should be positive.\n\n" );
        exit( EXIT_FAILURE );
    }
    for ( int i_block = 0; i_block < params1d->ic_blocks_number; i_block++ ) { // ���� �� ������ ��������� �������
        fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� - ����� �����
        fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->cell_begin[i_block]) ); // ����� ������-������
        if ( params1d->cell_begin[i_block] < 0 ) {
            printf( "\nread_Parameters1d -> cell_begin[%d] should be non-negative.\n\n", i_block );
            exit( EXIT_FAILURE );
        }
        fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->cell_end[i_block]) ); // ����� ������-�����
        if ( params1d->cell_end[i_block] < 0 ) {
            printf( "\nread_Parameters1d -> cell_end[%d] should be non-negative.\n\n", i_block );
            exit( EXIT_FAILURE );
        }

        // ������ ����������� ���������� � �����
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->block_values[i_block][B_DISP]) );
        if ( params1d->block_values[i_block][B_DISP] < 0.0 ) {
            printf( "\fill_parameters_1d -> block_values[%d][%d] should be non-negative.\n\n", i_block, B_DISP );
            exit( EXIT_FAILURE );
        }
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->block_values[i_block][R_DISP]) );
        if ( params1d->block_values[i_block][R_DISP] < 0.0 ) {
            printf( "\fill_parameters_1d -> block_values[%d][%d] should be non-negative.\n\n", i_block, R_DISP );
            exit( EXIT_FAILURE );
        }
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->block_values[i_block][V_DISP]) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->block_values[i_block][P_DISP]) );
        if ( params1d->block_values[i_block][P_DISP] < 0.0 ) {
            printf( "\fill_parameters_1d -> block_values[%d][%d] should be non-negative.\n\n", i_block, P_DISP );
            exit( EXIT_FAILURE );
        }
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->block_values[i_block][R_GAS]) );
        if ( params1d->block_values[i_block][R_GAS] <= 0.0 ) {
            printf( "\fill_parameters_1d -> block_values[%d][%d] should be positive.\n\n", i_block, R_GAS );
            exit( EXIT_FAILURE );
        }
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->block_values[i_block][V_GAS]) );
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->initial_burning[i_block]) );// ����� �� ���� ����� ������
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->S[i_block]) );// ������� ������� �����
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->block_values[i_block][P_GAS]) );
        if ( params1d->block_values[i_block][P_GAS] <= 0.0 ) {
            printf( "\fill_parameters_1d -> block_values[%d][%d] should be positive.\n\n", i_block, P_GAS );
            exit( EXIT_FAILURE );
        }
        fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
        for (int i_scalar = 0; i_scalar < params1d->number_of_scalars; i_scalar++){
            fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->scalar_values[i_block][i_scalar]) ); // �������������� �������
            if (params1d->number_of_scalars > 0 && params1d->scalar_values[i_block][i_scalar] < 0.0){
                printf( "\nfill_parameters_struct -> >scalar_values[%d][%d] should be a non-negative value\n\n", i_block, i_scalar );
                exit( EXIT_FAILURE );
            }
        }
    }

    // ��������� �������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->left_bc) ); // ����� ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->right_bc) ); // ������ ��������� �������
    // ������ ����������� ���������� ��� ���������� ������� ��������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[B_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[R_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[V_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[P_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[R_GAS]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[V_GAS]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[P_GAS]) );
    for (int i_scalar = 0; i_scalar < params1d->number_of_scalars; i_scalar++)
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->inflow_values[Z0 + i_scalar]) ); 

    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->change_inflow_time) ); // ������ ������� ������������ �� ������ ��������� ��������
    if ( ( params1d->left_bc == COMPLEX_INFLOW || params1d->left_bc == COMPLEX_INFLOW ) && params1d->change_inflow_time <= 0.0 ) {
        printf( "\nfill_parameters_struct -> change_inflow_time should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
    
    // ������ ������ ����������� ���������� ��� ��������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->changed_inflow_values[B_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->changed_inflow_values[R_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->changed_inflow_values[V_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->changed_inflow_values[P_DISP]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->changed_inflow_values[R_GAS]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->changed_inflow_values[V_GAS]) );
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->changed_inflow_values[P_GAS]) );

    // ������ ������� ��������� �����-��������, ����� ����� ������ ��� ��������� ������� �� ���� ������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ��������� �� ������� ������ �������
    switch ( dtmp ) {
        case 1:
            params1d->build_exact_solution = true;
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\nfill_parameters_struct -> Warning. The exact solution will be not constructed because it is constructed for Baer-Nunziato model only.\n\n" );
            }
            if ( ( params1d->ic_blocks_number != 2 ) ||
                 ( ( params1d->cell_end[0] - params1d->cell_end[0] ) != ( params1d->cell_end[1] - params1d->cell_end[1] ) ) ) {
                printf( "\nfill_parameters_struct -> Warning. The exact solution will be not constructed because the blocks number is not equal to 2 of the blocks are of different size.\n\n" );
            }
            break;
        case 0:
            params1d->build_exact_solution = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> build_exact_solution should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    // ���������� ����� ��� ���������� ������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(params1d->cells_number_for_exact_solution) );
    if ( params1d->build_exact_solution ) {
        if ( params1d->cells_number_for_exact_solution <= 0 ) {
            printf( "\nfill_parameters_struct -> cells_number_for_exact_solution should be a positive value\n\n" );
            exit( EXIT_FAILURE );
        }
    }

    // �������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    if ( paramsc->are_sensors ) {
        for ( int i_sensor = 0; i_sensor < paramsc->sensors_num; i_sensor++ ) {
            fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(params1d->sensors_location[i_sensor]) ); // ���������� �������
        }
    }

}

// ���������� ��������� � ����������� ���������� ������ � ���������� ���������� �������� �����, � ����� ��������� ���������� ��������� ������.
// �������� ������������ ������� ����������.
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
// input_file_directory - ����������, � ������� ��������� ���� � ����������� ������
// output_file_directory - ����������, � ������� ������ ���� �������� ����� � ������������
// file_num - ������� ����� ����� Parameters1d.dat ��� ����������
void fill_parameters_struct_1d( struct ParametersCommon* paramsc, struct Parameters1d* params1d, const char* input_file_directory, const char* output_file_directory, const int file_num ) {

    FILE *params; // ���������� ����������������� ����� ������
    char string[MAX_STRING_SIZE]; // ��� ���������� ��������� ���������� �� �����
            
    // ��� ��������� ������ ��������� � ����� parameters.dat
    strcpy_s( string, input_file_directory );
    strcat_s( string, "\\parameters" );
    if ( file_num != 0 ) { // ���� 0, �� ������ �� ���������� ��� ���������� ��������������� � ���������� ���������
        char tmp_str[MAX_STRING_SIZE]; // ��������� ���������� ��� �������� �������� ������
        sprintf_s( tmp_str, "%d", file_num );
        strcat_s( string, tmp_str );
    }
    strcat_s( string, ".dat" );
    if ( ( fopen_s( &params, string, "rt" ) ) != 0 ) {
        printf( "\nfill_parameters_struct_1d -> can't open file %s for reading\n\n", string );
        exit( EXIT_FAILURE );
    }

    // ��������� �������� ����������, � ������� ��������� ���� � ����������� ������, � params
    strcpy_s( paramsc->input_file_directory, input_file_directory );
    
    // ��������� �������� ����������, ���� ����� �������� ���������� ��������, � params
    strcpy_s( paramsc->output_file_directory, output_file_directory );
    
    // ���������� ��������� �����
    fscanf_s( params, "%s", string, MAX_STRING_SIZE );

    // ���������� ����� ���������� �������
    fill_parameters_common( params, paramsc );

    // ���������� ���������� ���������� ������
    fill_parameters_1d( params, paramsc, params1d );

    fclose( params );

}

// ������� ���������� ������ �����, � ������� ����������� �������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ ��������������� ������������ (in)
// nodes_coord - ������ ��������� ����� (in)
// sensors_cells - ������ ������� �����, � ������� ����������� ������� (out)
void get_sensors_cells( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *nodes_coord, int *sensors_cells ) {

    // �������������
    for ( int i_sensor = 0; i_sensor < paramsc->sensors_num; i_sensor++ )
        sensors_cells[i_sensor] = -1;

    // ��������� ������������
    for ( int i_sensor = 0; i_sensor < paramsc->sensors_num; i_sensor++ ) {
        for ( int i_node = 0; i_node < params1d->cells_number; i_node++ ) {
            if ( params1d->sensors_location[i_sensor] >= nodes_coord[i_node] &&
                 params1d->sensors_location[i_sensor] <= nodes_coord[i_node+1])
                 sensors_cells[i_sensor] = i_node;
        }
    }

    // ��������
    for ( int i_sensor = 0; i_sensor < paramsc->sensors_num; i_sensor++ ) {
        if ( sensors_cells[i_sensor] < 0 ) {
            printf( "\nget_sensors_cells -> wrong location of sensor %d\n", i_sensor );
            exit( EXIT_FAILURE );
        }
    }

}

// ������������ ������������ ����� ����������� �������� � �������� �����,
// � ����� �������� ������ ��� �������� �� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ (in)
// nodes_coord - ������ ��������� ����� (in)
// sensors_cells - ������ ������� �����, � ������� ����������� ������� (out)
// sensors_files[MAX_SENSORS_NUM] - ������ �������� ������������ ��� ������ ������ � �������� (out)
void prepare_sensors( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *nodes_coord, int *sensors_cells, FILE *sensors_files[MAX_SENSORS_NUM] ) {

    char curr_filename[MAX_STRING_SIZE]; // ��� �������� �����
    char tmp_str[MAX_STRING_SIZE]; // ��� ������������� ������ ������� � ������

    // ������������ ������������ ����� ����������� �������� � �������� �����
    get_sensors_cells( paramsc, params1d, nodes_coord, sensors_cells );

    // �������� ������ ��� ������ ������ � ��������
    for ( int i_sensor = 0; i_sensor < paramsc->sensors_num; i_sensor++ ) {
        // ������������ ����� �����
        strcpy_s( curr_filename, "sensor_" );
        sprintf_s( tmp_str, "%d", i_sensor );
        strcat_s( curr_filename, tmp_str );
        strcat_s( curr_filename, ".dat" );
        // �������� ����� ��� ������ � ��������� �������
        fopen_s( &(sensors_files[i_sensor]), curr_filename, "wt" );
        if ( NULL == sensors_files[i_sensor] ) {
            printf( "\nwrite_sensors -> can't open file %s for writing\n", curr_filename );
        }
    }

}

// ������ ������� � ����
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ (in)
// output_directory - ����������, ���� ����� ������� ���� � ������ �������� (in)
// cells_number - ����� ����� (in) 
// *xc - ������ ��������� ������� ����� (in)
// **v_ncons - ��������� ������ ����������� ���������� � ������� ���� ����� ����� (in)
// file_number - ����� ����� � �������������� ������������ (in)
// output_filename - ��� �������� ����� � �������������� ������������ (out)
// n - �������� ������ �������� (in)
void write_solution_1d2phc( struct ParametersCommon *paramsc, struct Parameters1d* params1d, char *output_directory, int cells_number, double *xc, double **v_ncons, int file_number, double time_mom, char *output_filename, double *beta, double *B, double *C, int n, int *ignition_flag, int *initial_ignition ) {

    FILE *f_out; // ���������� ����� � ������������
    char tmp_str[MAX_STRING_SIZE]; // ��������� ���������� ��� �������� �������� ������
        
    // ������������ ����� ��������� ����� � ��������� ��������
    strcpy_s( output_filename, MAX_STRING_SIZE, output_directory );
    strcat_s( output_filename, MAX_STRING_SIZE, "\\numerical_solution_" );
    // ���������� ����������� ������
    if ( file_number < TEN_TO_FIFTH )
        strcat_s( output_filename, MAX_STRING_SIZE, "0" );
    if ( file_number < TEN_TO_FOURTH )
        strcat_s( output_filename, MAX_STRING_SIZE, "0" );
    if ( file_number < TEN_TO_THIRD )
        strcat_s( output_filename, MAX_STRING_SIZE, "0" );
    if ( file_number < TEN_TO_SECOND )
        strcat_s( output_filename, MAX_STRING_SIZE, "0" );
    if ( file_number < TEN )
        strcat_s( output_filename, MAX_STRING_SIZE, "0" );
    sprintf_s( tmp_str, "%d", file_number );
    
    // ������ ������ �����
    strcat_s( output_filename, MAX_STRING_SIZE, tmp_str );
    strcat_s(output_filename, MAX_STRING_SIZE, "time");
    sprintf_s( tmp_str, "%lf", time_mom );
    strcat_s( output_filename, MAX_STRING_SIZE, tmp_str );
    strcat_s( output_filename, MAX_STRING_SIZE, ".dat" );
    
    // �������� ����� ��� ������ � ��������� �������
    fopen_s( &f_out, output_filename, "wt" );
    if ( NULL == f_out ) {
        printf( "\nwrite_solution -> can't open file %s for writing\n", output_filename );
    }

    /* ������ ��������� ��� Tecplot � ������ ������������� */
    if ( paramsc->is_tecplot_header )
        write_solution_tecplot_header( paramsc, params1d, params1d->cells_number, f_out, file_number );

    // ������ �����������
    for ( int i = 0; i < cells_number; i++ ) {


        
        fprintf( f_out, "%e", xc[i] );
        for ( int j = 0; j < n; j++ ) {
            if ( j == R_DISP && v_ncons[i][B_DISP] < paramsc->eps_cut_out && paramsc->nice_output ) {
                /* ��� ������ ������ ���������� ���������� ���� � ������ ���������� �������� ������� �������� */
                fprintf( f_out, " %e", paramsc->background_density );
                continue;
            }
            if ( j == P_DISP && v_ncons[i][B_DISP] < paramsc->eps_cut_out && paramsc->nice_output ) {
                /* ��� ������ ������ ���������� ���������� ���� � ������ ���������� �������� ������� �������� */
                fprintf( f_out, " %e", paramsc->background_pressure );
                continue;
            }
            if ( j == V_DISP && v_ncons[i][B_DISP] < paramsc->eps_cut_out && paramsc->nice_output ) {
                /* ��� ������ ������ ���������� ���������� ���� � ������ ���������� �������� ������� �������� */
                fprintf( f_out, " %e", paramsc->background_velocity );
                continue;
            }
	    if (j != P_DISP)
		fprintf( f_out, " %e", v_ncons[i][j] );
	    else
	//	fprintf( f_out, " %e", v_ncons[i][P_DISP] - beta[i] );
		fprintf( f_out, " %e", v_ncons[i][P_DISP] );
			
        }
        
	fprintf( f_out, " %e", beta[i] );
	fprintf( f_out, " %e", T_disp(paramsc, v_ncons[i][P_DISP], v_ncons[i][R_DISP]));
	fprintf( f_out, " %e", B[i] );
        fprintf( f_out, " %e", C[i] );
	fprintf(f_out," %e",T_gas(paramsc, v_ncons[i][P_GAS], v_ncons[i][R_GAS]));
        fprintf( f_out, "\n" );
    }

    fclose( f_out );

}

// ������ � ���� � ������������ ������������ ���������� ��� Tecplot
// params - ��������� � ��������� ����������� ��������������� ������������ (in)
// cells_num - ����� ����� (in)
// file_to_write - ���������� �����, ��������� ��� ������ (in)
// file_number - ����� ����� � �������������� ������������ (in)
void write_solution_tecplot_header( struct ParametersCommon *paramsc, struct Parameters1d *params1d, int cells_num, FILE *file_to_write, int file_number ) {

    char tmp_str[MAX_STRING_SIZE];  /* ��������� ���������� ��� ������ ������ */
    double time_gap;                /* ����� ����� �������� ������ � �������������� ������������,
                                       ����� ����� ��� params->is_output_on_time = 1 */
    int steps_gap;                  /* ���������� ����� �� ������� ����� �������� ������ � �������������� ������������,
                                       ����� ����� ��� params->is_output_on_time = 0 */

    /* �������� ������ � �������� ������ � ����������� */
    fprintf( file_to_write, "TITLE = \"1D Baer-Nunziato solver\"\nVARIABLES = " );
    
    /* �������� ���������� */
    /* ���������� */
    fprintf( file_to_write, "\"X\"\n" );
    /* �������� ���� ���������� ���� */
    fprintf( file_to_write, "\"Solid phase volume fraction\"\n" );
    /* ��������� ���������� ���� */
    fprintf( file_to_write, "\"Solid phase density\"\n" );
    /* �������� ���������� ���� */
    fprintf( file_to_write, "\"Solid phase velocity\"\n" );
    /* �������� ���������� ���� */
    fprintf( file_to_write, "\"Solid phase pressure\"\n" );
    /* ��������� ������� ���� */
    fprintf( file_to_write, "\"Gas phase density\"\n" );
    /* �������� ������� ���� */
    fprintf( file_to_write, "\"Gas phase velocity\"\n" );
    /* �������� ������� ���� */
    fprintf( file_to_write, "\"Gas phase pressure\"\n" );
	
    if(paramsc->program_name == ONED2PHC){
        /* ���������� ������� */
        for (int i = 0; i < params1d->number_of_scalars; i++)
        fprintf( file_to_write, "\"Scalar %d\"\n", i );
	/* �������� ��������������� */
	fprintf( file_to_write, "\"Configuration pressure\"\n" );
	/* ����������� ���������� ���� */
	fprintf( file_to_write, "\"Solid phase temperature\"\n" );
	/* ������� ��������������� */
	fprintf( file_to_write, "\"Compaction energy\"\n" );
        /* ���������� ������� */
	fprintf( file_to_write, "\"Chemical reaction\"\n" );
	/* ����������� ������� ���� */
	fprintf( file_to_write, "\"Gas phase temperature\"\n" );
    }

    if (paramsc->program_name == ONED3PHC){
	/* ���������� �������� ���������� ���� g */
	fprintf( file_to_write, "\"Solid phase g\"\n" );
	/* �������� P0 ����������� ��� ���������� ���� g */
	fprintf( file_to_write, "\"Solid phase P0\"\n" );
	/* ���������� �������� ������� ���� g */
	fprintf( file_to_write, "\"Gas phase g\"\n" );
	/* �������� P0 ����������� ��� ������� ���� g */
	fprintf( file_to_write, "\"Gas phase P0\"\n" );
	/* �������� �� �������� ���������� ���� � ������ ������ */
	fprintf( file_to_write, "\"Solid phase physical state\"\n" );
	/* �������� �� �������� ������� ���� � ������ ������ */
	fprintf( file_to_write, "\"Gas phase physical state\"\n" );
    }
    
    /* ������������ �������� ���� */
    if ( paramsc->is_output_on_time ) {
        /* ������ ������ �� ���������� ��������� ������� ������� */
        time_gap = paramsc->stop_time / paramsc->output_number;
        sprintf_s( tmp_str, "%f", time_gap * file_number );
        strcat_s( tmp_str, " time units" );
    }
    else {
        /* ������ ������ �� ���������� ��������� ���������� ����� */
        steps_gap = paramsc->stop_steps_number / paramsc->output_number;
        sprintf_s( tmp_str, "%d", steps_gap * file_number );
        strcat_s( tmp_str, " time steps" );
    }
    /* ������ �������� ���� */
    fprintf( file_to_write, "ZONE T=\"" );
    fprintf( file_to_write, "%s", tmp_str );
    fprintf( file_to_write, "\"\n" );

    /* ��������� ������� ���������� */
    fprintf( file_to_write, "STRANDID=0, SOLUTIONTIME=0\n" );

    /* ���������� �� ��������� ����� */
    fprintf( file_to_write, "I=%d, J=1, K=1, ZONETYPE=Ordered\n", cells_num );

    /* ������ �������������� ������ */
    fprintf( file_to_write, "DATAPACKING=POINT\n" );

    /* �������� ������������� ������ */
    if (paramsc->program_name == ONED2PHC){
	fprintf( file_to_write, "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE" );
        for (int i = 0; i < params1d->number_of_scalars; i++)
            fprintf( file_to_write, " SINGLE" );
        fprintf( file_to_write, ")\n" );
    }
    else if (paramsc->program_name == ONED3PHC){
	fprintf( file_to_write, "DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )\n" );
    }
}

// ������ ������ �� �������������� ����������� � ������ ��������� ���������� �������� � ��������� ������ ��������
// paramcs - ��������� � ��������� ����������� ��������������� ������������
// paramc1d - ��������� � ����������� ���������� ������
// time_mom - ���������� � ������� ������� �������
// debug_info - ���������� ����������
// u_prev - ������ ����������� ���������� �� ������� ����
void write_statistics( const struct ParametersCommon *paramsc, struct Parameters1d *params1d, const struct TimeMoment *time_mom,const struct DebugInfo *debug_info, double **u_prev ) {

    // ����������� ���������� �����, � ������� ������������ ���������� ����
    int disp_exist = 0;
    for ( int i_cell = 0; i_cell < params1d->cells_number; i_cell++ )
        if ( u_prev[i_cell][B_DISP] > paramsc->eps_disp_abs )
            disp_exist++;

    // ������ ���������� � ����� ���������� �� ������ ����
    fprintf( debug_info->relaxation_out, "%d %e %d %d %d %d %d\n", time_mom->steps_num, time_mom->curr_t, debug_info->relaxation_cases[0], debug_info->relaxation_cases[1],
        debug_info->relaxation_cases[2], debug_info->relaxation_cases[3], debug_info->relaxation_cases[4] );
    fprintf( debug_info->godunov_out, "%d %e %d %d %d %2.2f\n", time_mom->steps_num, time_mom->curr_t, debug_info->godunov_cases[0], debug_info->godunov_cases[1],
        debug_info->godunov_cases[2], (double)disp_exist / (double)params1d->cells_number * 100.0 );

}

// ������ � ���� �������� ������� ���������� � ���������� ������ � ������ ������ �������
// paramcs - ��������� � ��������� ����������� ��������������� ������������
// paramc1d - ��������� � ����������� ���������� ������
// curr_t - ������� ������ �������
// file_to_write - ���������� ����� ��� ������
// file_name - ��� ������������ �����
// v_ncons - ������ ����������� ���������� �� ������� ����
// number_of_variables - ���������� ���������� � ������� v_ncons
void write_time_dependent_information_in_a_cell (const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, double curr_t, char *output_directory, char *file_name, double *v_ncons, int number_of_variables) {
    
    FILE *file_to_write;
    
    char final_file_name[MAX_STRING_SIZE];
    // ������������ ����� ��������� ����� � ��������� ��������
    strcpy_s( final_file_name, MAX_STRING_SIZE, output_directory );
    strcat_s( final_file_name, MAX_STRING_SIZE, "//" );
    strcat_s( final_file_name, MAX_STRING_SIZE, file_name );

    fopen_s(&file_to_write,"final_file_name", "a");

    fprintf_s(file_to_write, "%.10lf \t", curr_t);
    for ( int i = 0; i < number_of_variables; i++ ){
        fprintf_s(file_to_write, "%lf \t", v_ncons[i]);
    }
    fprintf_s(file_to_write, "\n", curr_t);

    fclose(file_to_write);
}

// ������ ��������� ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// v_ncons - ������� ���������������� ���������� �� ���� ������� ��������� ������� (in)
// time_mom - ���������, ������������ ������� ������ ������� (in)
// sensors_cells - ������ ������� �����, � ������� ����������� ������� (in)
void write_sensors( struct ParametersCommon *paramsc, double **v_ncons, struct TimeMoment *time_mom, FILE *sensors_files[MAX_SENSORS_NUM],int sensors_cells[MAX_SENSORS_NUM], int n ) {

    for ( int i_sensor = 0; i_sensor < paramsc->sensors_num; i_sensor++ ) {
        // ������ ������� ������� � ����
        fprintf_s( sensors_files[i_sensor], "%e\t%d\t", time_mom->curr_t, time_mom->steps_num );
        // ������ ������
        for ( int i_component = 0; i_component < n; i_component++ )
            fprintf_s( sensors_files[i_sensor], "%e\t", v_ncons[sensors_cells[i_sensor]][i_component] );
        fprintf_s( sensors_files[i_sensor], "%e\t", v_ncons[sensors_cells[i_sensor]][B_DISP] * v_ncons[sensors_cells[i_sensor]][P_DISP] + (1.0 - v_ncons[sensors_cells[i_sensor]][B_DISP]) * v_ncons[sensors_cells[i_sensor]][P_GAS] );
        fprintf_s( sensors_files[i_sensor], "\n" );
    }

}