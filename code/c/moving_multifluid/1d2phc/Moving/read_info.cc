#include "read_info.h"

// ���������� ����� � ����������� ������, ���������� ����� ��������������� ��������� � �������� ������������ ������� ����������
// params - ��������� � ����������� ��������������� ������������
void read_parameters( struct Parameters *params )
{
    FILE *parameters_file;
    char string[MAX_STRING_SIZE]; // ��� ���������� ��������� ���������� �� �����
    int i_block; // ��� ���������� ������ � ���������� ���������
        
    // ��� ��������� ������ ��������� � ����� parameters2D.txt
    if ( ( fopen_s( &parameters_file, "parameters.dat", "rt" ) ) != 0 )
    {
        printf( "\nread_parameters -> can't open file parameters.dat for reading\n\n" );
        system("Pause");
    }

    // ���������� ��������� �����
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );

    // ��������� �����
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->number_of_cells) ); // ���������� �����
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->coordinate_of_left_boundary) ); // ���������� ����� ������� ��������� �������
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->coordinate_of_right_boundary) ); // ���������� ������ ������� ��������� �������

    // ��������� ����
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->initial_coordinate_of_left_boundary_of_body) ); // ���������� ����� ������� ����
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->initial_coordinate_of_right_boundary_of_body) ); // ���������� ������ ������� ����
    if ( params->initial_coordinate_of_right_boundary_of_body - params->initial_coordinate_of_left_boundary_of_body < 
        ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary ) / params->number_of_cells )
    {
        printf( "\nread_parameters -> too small body\n" );
        system("Pause");
    }

    // �������� ����
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->body_velosity) ); // �������� ����
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->body_cross_section_devided_to_mass) );

    // ��� �� ������� 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );              // ��������� ������� 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->exit_time_cycle) );    // ������� ������ �� ����� �� �������
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->time_cycle_iterations) );    // ����� ������� � ����
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->t_fin) );    // �������� ������ �������
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->time_step_method) );   // ����� ������� ���� �� ������
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->cfl) );   // ����� cfl ���������� ��� �� �������
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->dt) );  // �������� ����

    // ��������� ����������������
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );              // ��������� ������� 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->length_diml) ); 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->time_diml) ); 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->mass_diml) ); 

    // ������ �������� ������
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );              // ��������� ������� 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->number_of_output_files) );
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->step_of_output_files) );

    // ��������� ����� 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );                                  // ��������� ������� 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->flux_calculation_method) );  // ����� ���������� ����� 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->max_iter_num) );      // ������������ ���������� �������� 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->p_max_ratio) );  // ������������ ������� �� �������� �����
                                                                                           //� ������ �� �������, ��� ������� ���
                                                                                           //������������ ��������� ����������� ��
                                                                                           //������� ��������������� ������
                                                                                           //� ������ ��������
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->approximation_order) ); // ������� ���������� ����� �� ������������
    if ( params->approximation_order != 1 && params->approximation_order != 2 )
    {
        printf( "\nread_parameters -> approximation_order should be 1 or 2.%d\n\n", params->approximation_order );
        system("Pause");
    }
	
    // ��������� ������� 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );                                  // ��������� ������� 
    for ( i_block = 0; i_block < 2; i_block++ ) // ���� �� ������ ��������� ������� - ������ � �������
    {
        fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // ��������� - ����� �����
        fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).number_of_left_cell ) );
        fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).number_of_right_cell ) );
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).initial_primitive_vector.density ) );
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).initial_primitive_vector.velosity ) );
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).initial_primitive_vector.pressure ) );
    }

    // ����� ���������
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );  // ��������� ������� 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->eps_general) );  // ���������� ����� ����� ��� ��������� ������������ �����,
                                                                                           //������������� � �������������� ��������, ��������
                                                                                           //���������� ������������ ��������� (���� ���� �� ��������� �����) 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->eps_spatial) ); // ����� ����� ��� �������� ������� ���������� ���� �����
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->big) ); // ������� �����

    // ���������� �������� 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );  // ��������� ������� 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->g) );  // ��������� ������� 

    // ��� ������� �� 4 ���� ���������� ������������� �������
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );  // ��������� ������� - ��� �������
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->type_of_left_border ) );
    if ( params->type_of_left_border == PARAMETERS_BEHIND_SHOCK_WAVE )
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->Mach_number ) );
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->type_of_right_border ) );
    if ( params->type_of_right_border == PARAMETERS_BEHIND_SHOCK_WAVE )
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->Mach_number ) );
    

    // ���� ����������� ������
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    int dtmp; // ��������������� ���������� ��� ���������� ����� �����
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ���������� ����� ����������� ������
    switch ( dtmp )
    {
        case 0:
            params->is_debug_regime = false;
            break;
        case 1:
            params->is_debug_regime = true;
            break;
        default:
            printf( "\nread_parameters -> debug should be 0 or 1.\n\n" );
            system("Pause");
    }
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(dtmp) ); // ���������� ����� ������� ������� � ������ ��� � ���������� �������
    switch ( dtmp )
    {
        case 0:
            params->is_�ontinue = false;
            break;
        case 1:
            params->is_�ontinue = true;
            break;
        default:
            printf( "\nread_parameters -> is_�ontinue should be 0 or 1.\n\n" );
            system("Pause");
    }
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ���������� ����� ����������� ������
    switch ( dtmp )
    {
        case 0:
            params->is_calculate_drag_force = false;
            break;
        case 1:
            params->is_calculate_drag_force = true;
            break;
        default:
            printf( "\nread_parameters -> is_calculate_drag_force should be 0 or 1.\n\n" );
            system("Pause");
    }
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ���������� ����� ������������� �� ������� �������
    if ( 0 == dtmp ) params->is_initiate_parameters_in_outer_cells = false;
    else if ( 1 == dtmp ) params->is_initiate_parameters_in_outer_cells = true;
    else
    {
        printf( "\nread_parameters -> is_initiate_parameters_in_outer_cells should be 0 or 1.\n\n" );
        system("Pause");
    }

    fclose( parameters_file );
}

void dimless_parameters( struct Parameters *params )
{
    double length_diml = params->length_diml;
    double time_diml = params->time_diml;
    double dens_diml = params->mass_diml / pow(length_diml, 3);
    double speed_diml = length_diml / time_diml;
    double press_diml = params->mass_diml / ( length_diml * pow(time_diml, 2));

    params->coordinate_of_left_boundary /= length_diml;
    params->coordinate_of_right_boundary /= length_diml;

    params->body_velosity /= speed_diml;
    params->body_cross_section_devided_to_mass /= pow( length_diml, 2 ) / params->mass_diml;

    params->initial_coordinate_of_left_boundary_of_body /= length_diml;
    params->initial_coordinate_of_right_boundary_of_body /= length_diml;

    params->dt /= time_diml;
    params->t_fin /= time_diml;

    for ( int i_block = 0; i_block < 2; i_block++ ) // ���� �� ������ ��������� ������� 
    {
        (params->blocks[i_block]).initial_primitive_vector.density /= dens_diml;
        (params->blocks[i_block]).initial_primitive_vector.velosity /= speed_diml;
        (params->blocks[i_block]).initial_primitive_vector.pressure /= press_diml;
    }
}

// ���� � ��������� �������� �������, �� ������� ��������� ������ �������, �������� ��������� �����
void read_last_time_moment( struct Parameters *params, struct TimeMoment *time_mom )
{
    FILE *in;
    if ( ( fopen_s( &in, "last_time_moment.dat", "rt" ) ) != 0 )
    {
        printf( "\nread_last_time_moment -> can't open file read_last_time_moment.dat for reading%d\n\n", params->is_�ontinue );
        system ( "Pause" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( in, "%06d", &(time_mom->steps_num) );
    fscanf_s( in, "%lf", &(time_mom->curr_t) );
    time_mom->curr_t /= params->time_diml; // ���������������� ���������� ������� �������
    fclose( in );
}

void read_solution ( struct Parameters *params, struct Conservative_vector *conservative, struct TimeMoment *time_mom )
{
}