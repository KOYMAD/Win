// io.cc
// ����� ������� �����/������, �� ��������� �� ����������� ������
// (c) ����� �����, 2018
// ������: 26 ������� 2018 �.

#include "io.h"

// ���������� ��������� � ��������� ����������� ��������������� ������������
// params - ���������� ����������������� ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
void fill_parameters_common( FILE *params, struct ParametersCommon* paramsc ) {

    int dtmp; // ��� ���������� ������������� ����������
    char string[MAX_STRING_SIZE]; // ��� ���������� ��������� ���������� �� �����

    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->program_name) );
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ���������, ������������, ��� ����������� ����� �����, �� ��������� �� ����������� ������

    // ������ �����
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->media_model) ); // ������ ���������� �����
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ���������� ���� - ����������� �������� ������ ��� ���?
    switch ( dtmp ) {
        case 0:
            paramsc->porous_body = false;
            break;
        case 1:
            paramsc->porous_body = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> porous_body should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }

    // ��� �� �������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ��� - ���������� ��� ���
    switch ( dtmp ) {
        case 0:
            paramsc->constant_time_step = false;
            break;
        case 1:
            paramsc->constant_time_step = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> constant_time_step should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->cfl) ); // ����� cfl
    if ( !paramsc->constant_time_step ) {
        if ( paramsc->cfl <= 0.0 ) {
            printf( "\nfill_parameters_struct -> cfl number should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->dt) ); // ���������� ��� �� �������
    if ( paramsc->constant_time_step ) {
        if ( paramsc->dt <= 0.0 ) {
            printf( "\nfill_parameters_struct -> dt should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }

    // ��������� �����
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->numerical_method) ); // ����� ���������� �����
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->max_iter_num) ); // ������������ ���������� �������� ��� ������ ��������
    if ( paramsc->numerical_method == GODUNOV ) {
        if ( paramsc->max_iter_num <= 0 ) {
            printf( "\nfill_parameters_struct -> max_iter_num should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p_max_ratio) );  // ������������ ������� �� �������� ����� � ������ �� �������, ��� ������� ���
                                                                                     // ������������ ��������� ����������� �� ������� ��������������� ������
                                                                                     // � ������ ��������
    if ( paramsc->numerical_method == GODUNOV ) {
        if ( paramsc->p_max_ratio <= 0.0 ) {
            printf( "\nfill_parameters_struct -> p_max_ratio should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }

    // ���������� ��������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ������������ �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ���������� �������� ��� ������ ��� ���? ������ ��� SA � BN � ����������������
    switch ( dtmp ) {
        case 0:
            paramsc->pressure_relaxation = false;
            break;
        case 1:
            paramsc->pressure_relaxation = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> pressure_relaxation should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ��������� �� ��������������� � ��������� ���������� ��������?
    switch ( dtmp ) {
        case 0:
            paramsc->pressure_relaxation_compaction = false;
            break;
        case 1:
            paramsc->pressure_relaxation_compaction = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> pressure_relaxation_compaction should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� ������� ��� ������� ����������������� �������� � ���������� �������� � ����������������
    if ( paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction ) {
        switch ( dtmp ) {
            case 1:
                paramsc->pressure_relaxation_compaction_formula = SCHWENDEMAN_COMPACTION;
                break;
            case 2:
                paramsc->pressure_relaxation_compaction_formula = SAUREL_COMPACTION;
                break;
            default:
                printf( "\nfill_parameters_struct -> pressure_relaxation_compaction_formula should be 1 or 2.\n\n" );
                exit( EXIT_FAILURE );
        }
    } 
    if (paramsc->program_name == TWOD2PHC && paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction && paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION){
        printf( "\nfill_parameters_struct -> pressure_relaxation_compaction_formula should be 2. 1 has not been implemented for 2d2phc yet.\n\n" );
        exit( EXIT_FAILURE );
    }

    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->tau_parameter) ); 
    printf( "\nfill_parameters_struct -> tau_parameter = %lf.\n\n", paramsc->tau_parameter );
    if (paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction && paramsc->pressure_relaxation_compaction_formula == SAUREL_COMPACTION){
        if (paramsc->tau_parameter <= 0.0){
            printf( "\nfill_parameters_struct -> tau_parameter should be positive\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->n_parameter) ); 
    if (paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction && paramsc->pressure_relaxation_compaction_formula == SAUREL_COMPACTION){
        if (paramsc->n_parameter <= 1.0){
            printf( "\nfill_parameters_struct -> n_parameter should be greater than 1.0\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->volume_fraction_compaction) ); 
    if (paramsc->pressure_relaxation_compaction){
        if (paramsc->volume_fraction_compaction < 0.0 || paramsc->volume_fraction_compaction > 1.0){
            printf( "\nfill_parameters_struct -> volume_fraction_compaction should be equal to [0;1]\n\n" );
            exit( EXIT_FAILURE );
        }
    }

    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ���������� ��������� ��� ������ ��� ���?
    switch ( dtmp ) {
        case 0:
            paramsc->velocity_relaxation = false;
            break;
        case 1:
            paramsc->velocity_relaxation = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> velocity_relaxation should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->approximation_order) ); // ������� ������������� ������
    if ( paramsc->approximation_order != 1 && paramsc->approximation_order != 2 ) {
        printf( "\nfill_parameters_struct -> approximation_order should be 1 or 2\n\n" );
        exit( EXIT_FAILURE );
    }
            
    // ��������� ��������� ���
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ��������� �������� ��� ��� ���������� ���������� ����������� ��� ��� ���?
    switch ( dtmp ) {
        case 0:
            paramsc->use_real_eos = false;
            break;
        case 1:
            paramsc->use_real_eos = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> use_real_eos should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }

    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� ����������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->g1) ); // ���������� �������� � ���������� ����
    if ( paramsc->g1 <= 1.0 ) {
        printf( "\nfill_parameters_struct -> g1 should be grater than 1.0.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p01) ); // �������� p01 � ��� ��� ���������� ����
    if ( paramsc->p01 < 0.0 ) {
        printf( "\nfill_parameters_struct -> p01 should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->material_number1) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->specific_heat_cv1) );
	
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� ����������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->g2) ); // ���������� �������� � ������� ����
    if ( paramsc->g2 <= 1.0 ) {
        printf( "\nfill_parameters_struct -> g2 should be grater than 1.0.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p02) ); // �������� p02 � ��� ��� ������� ����
    if ( paramsc->p02 < 0.0 ) {
        printf( "\nfill_parameters_struct -> p02 should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->b_virial) ); // �������� b � ���������� ��� ��� ������� ����
    if ( paramsc->b_virial < 0.0 ) {
        printf( "\nfill_parameters_struct -> b should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mu2) ); // ����������� ������������ �������� � ����
    if ( paramsc->mu2 < 0.0 ) {
        printf( "\nfill_parameters_struct -> mu_g should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->material_number2) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->specific_heat_cv2) );

    // ����� �������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp );
    switch ( dtmp ) {
        case 1:
            paramsc->is_debug = true; // ������ � ���������� ��������� �������� � ������������ ������� ������ ���������
            break;
        case 0:
            paramsc->is_debug = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_debug should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }

    // �����
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp );
    switch ( dtmp ) {
        case 1:
            paramsc->is_output_on_time = true; // �������� ������ ����������� - ���������� ��������� ������� �������
            break;
        case 0:
            paramsc->is_output_on_time = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_output_on_time should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->stop_time) ); // ������ ������� ��������� �������
    if ( paramsc->is_output_on_time ) {
        if ( paramsc->stop_time <= 0.0 ) {
            printf( "\nfill_parameters_struct -> stop_time should be a positive value.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->stop_steps_number) ); // ��������� ���������� ����� �� �������
    if ( !paramsc->is_output_on_time ) {
        if ( paramsc->stop_steps_number <= 0.0 ) {
            printf( "\nfill_parameters_struct -> stop_steps_number should be a positive value.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->output_number) ); // ���������� ������ � �������������� ������������
    if ( paramsc->output_number <= 0 ) {
        printf( "\nfill_parameters_struct -> output_number should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // �������� �� ��������� ��� ������������ � Tecplot?
    switch ( dtmp ) {
        case 1:
            paramsc->is_tecplot_header = true;
            break;
        case 0:
            paramsc->is_tecplot_header = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_tecplot_header should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    // �������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� ����������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ��������� �� ������������ ������� � �������� ������?
    switch ( dtmp ) {
        case 1:
            paramsc->are_sensors = true;
            break;
        case 0:
            paramsc->are_sensors = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> are_sensors should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    if ( paramsc->are_sensors ) { 
        fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->sensors_num) ); // ����������� ���������� ��������
    }
    if ( paramsc->are_sensors && paramsc->sensors_num <= 0 ) {
        printf( "\nfill_parameters_struct -> sensors_num should be positive.\n\n" );
        exit( EXIT_FAILURE );
    }
   
    // ����� ���������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_general) ); // ���������� ����� ����� ��� ��������� ������������ �����,
                                                                                    // ������������� � �������������� ��������, ��������
                                                                                    // ���������� ������������ ��������� (���� ���� �� ��������� �����)
    if ( paramsc->eps_general <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_general should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_ludcmp) ); // ����� ����� � ������� ludcmp (math_utilc.cc),
                                                                                   // ��������� ��� ���� � ������������ ������
    if ( paramsc->eps_ludcmp <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_ludcmp should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_decouple) ); // ����������� ���������� ������� �� �������� ���� ����� � ������,
                                                                                     // ��� ������� ����������� ���������������� �����, ����� - �����������
                                                                                     // ������������ ����� ��������� ��� ���
    if ( paramsc->eps_decouple <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_decouple should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_thin_layer) ); // �������� ������� ������� ���������� �������������� ���������
                                                                                       // "������� ����" ��� ��������
    if ( paramsc->eps_thin_layer <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_thin_layer should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_contact) ); // ����� �������� ��� ������� �������� � ������ �������������
                                                                                    // ����������� �������
    if ( paramsc->eps_contact <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_contact should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_disp_abs) ); // ����� �������� �������� ���� ����������
                                                                                     // ����, ������� � ������� ���������, ��� ���
                                                                                     // �����������
    if ( paramsc->eps_disp_abs <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_disp_abs should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_cut_out) ); // ����� �������� �������� ���� ����������
                                                                                    // ����, ������� � ������� ��� ������
                                                                                    // ����������� ��������� ���������� ����
                                                                                    // �������� ������� ������� ����������
    if ( paramsc->eps_cut_out <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_cut_out should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }

    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_relax_compaction) ); // ����� ����� ��� ��������� ���������� 
                                                                                             // ��� ��������� � ���������� �������� 
                                                                                             // � ����������������, ������ ��������� ����
    if ( paramsc->eps_relax_compaction <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_relax_compaction should be a positive value.\n\n" ); 
        exit( EXIT_FAILURE );
    }

    // ��������� ����������������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ������ ����� ������������ ���������� ����������� ��� ���?
    switch ( dtmp ) {
        case 1:
            paramsc->use_dimensions = true;
            break;
        case 0:
            paramsc->use_dimensions = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> use_dimensions should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mass_scale) ); // ����������� ������� �����
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->time_scale) ); // ����������� ������� �������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->length_scale) ); // ����������� ������� �����
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->temperature_scale) ); // ����������� ������� �����������

    // ������� �������� ���������� ���������� ���� ��� ����� � ����� ��������� �������� ���� ���������� ����
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->background_density) ); // ������� �������� ��������� ���������� ����
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->background_velocity) ); // ������� �������� �������� ���������� ����
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->background_pressure) ); // ������� �������� �������� ���������� ����
   
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ���������� � ������ � ����� ����������� ���������� ����
                                                                 // ������� ��������� ��� ������ �����������?
    switch ( dtmp ) {
        case 1:
            paramsc->nice_output = true;
            break;
        case 0:
            paramsc->nice_output = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> nice_output should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }

    // "������"
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� �� ����������� �����, ����������� ������?
    switch ( dtmp ) {
        case 1:
            paramsc->is_physics = true;
            break;
        case 0:
            paramsc->is_physics = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_physics should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->substeps_num) ); // ���������� ��������, �� �������
                                                                                    // ����������� ���������������� ���, ���
                                                                                    // �������������� ����, �����������
                                                                                    // ��������� ��������������
    if ( paramsc->is_physics && paramsc->substeps_num <= 0.0 ) {
        printf( "\nfill_parameters_struct -> substeps_num should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->particle_diameter) ); // ������� ������ ���������� ����
    if ( paramsc->is_physics && paramsc->particle_diameter <= 0.0 ) {
        printf( "\nfill_parameters_struct -> particle_diameter should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->numerical_method_physics) );
    if ( paramsc->is_physics && ( paramsc->numerical_method_physics != EXPLICIT_EULER && paramsc->numerical_method_physics != NEWTON ) ){
        printf( "\nfill_parameters_struct -> numerical_method_physics should be 1 or 2.\n\n" );
        exit( EXIT_FAILURE );
    }

    // ���� ���������� ������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� �� ����������� ���� ���������� ������?
    if ( paramsc->is_physics ) {
        switch ( dtmp ) {
            case 1:
                paramsc->is_friction_force = true;
                break;
            case 0:
                paramsc->is_friction_force = false;
                break;
            default:
                printf( "\nfill_parameters_struct -> is_friction_force should be 0 or 1.\n\n" );
                exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� ������� ��� ������� ���� ���������� ������
    if ( paramsc->is_physics && paramsc->is_friction_force ) {
        switch ( dtmp ) {
            case 1:
                paramsc->fric_force_formula_num = ROGUE;
                break;
            case 2:
                paramsc->fric_force_formula_num = HOUIM;
                break;
            case 3:
                paramsc->fric_force_formula_num = TANINO;
                break;
	    case 4:
		paramsc->fric_force_formula_num = SCHWENDEMAN;
                break;
            case 5:
                paramsc->fric_force_formula_num = HOMENKO;
                break;

            default:
                printf( "\nfill_parameters_struct -> fric_force_formula_num should be 1, 2, 3, 4 or 5.\n\n" );
                exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->interface_drag_coef) );
    if ( paramsc->program_name == ONED3PHC || TWOD2PHC){
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // �������� �� ����������� Cd ����������?
	if ( paramsc->is_physics ) {
	    switch ( dtmp ) {
		case 1:
		    paramsc->is_Cd_const = true;
		    break;
		case 0:
		    paramsc->is_Cd_const = false;
		    break;
		default:
		    printf( "\nfill_parameters_struct -> is_Cd_const should be 0 or 1.\n\n" );
		    exit( EXIT_FAILURE );
	    }
	}
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->Cd) );
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->Cd_formula) );
    }

    // ���������������
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� �� ����������� ���������������?
    if ( paramsc->is_physics ) {
        switch ( dtmp ) {
            case 1:
		paramsc->is_compaction = true;
                break;
            case 0:
                paramsc->is_compaction = false;
                break;
            default:
                printf( "\nfill_parameters_struct -> is_compaction should be 0 or 1.\n\n" );
                exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->compaction_viscosity) );

    // ������� � ������������ � ���������� ���������
    if (paramsc->program_name == ONED2PHC || paramsc->program_name == TWOD2PHC){
	// �������
	fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� �� ����������� �������?

	if ( paramsc->is_physics ) {
	    switch ( dtmp ) {
		case 1:
		    paramsc->is_burning = true;
		    break;
		case 0:
		    paramsc->is_burning = false;
		    break;
		default:
		    printf( "\nfill_parameters_struct -> is_burning should be 0 or 1.\n\n" );
		    exit( EXIT_FAILURE );
	    }
	}
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->ignition_condition)); // ����� ������� �������������

        fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� ������� ��� ������� �������
        if ( paramsc->is_physics && paramsc->is_burning ) {
            switch ( dtmp ) {
                case 1:
                    paramsc->burning_formula_num = SEREBRYAKOV_BURNING;
                    break;
                case 2:
                    paramsc->burning_formula_num = SCHWENDEMAN_BURNING;
                    break;
                default:
                    printf( "\nfill_parameters_struct -> burning_formula_num should be 1 or 2.\n\n" );
                    exit( EXIT_FAILURE );
            }

        }
        // ��������� ��� ������� SEREBRYAKOV_BURNING
        fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->e_coef) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->U_coef) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->nu_coef) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->X_coef1) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->lambda_coef1) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mu_coef1) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->X_coef2) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->lambda_coef2) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mu_coef2) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->zk) );
        // ��������� ��� ������� SCHWENDEMAN_BURNING
        fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->reaction_rate_prefactor) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p_ignition) );
        // ����������� ����� ���������� �������
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->heat_release) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->T_ignition) );
	// ������������
	fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // ��������� �������
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // ����� �� ����������� ������������?
	if ( paramsc->is_physics ) {
	    switch ( dtmp ) {
		case 1:
		    paramsc->is_heat_transfer = true;
		    break;
		case 0:
		    paramsc->is_heat_transfer = false;
		    break;
		default:
		    printf( "\nfill_parameters_struct -> is_heat_transfer should be 0 or 1.\n\n" );
		    exit( EXIT_FAILURE );
	    }
	}
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->heat_transfer_coefficient) );
    }
}

// ������ ����� � ����������� ��� ����������� ����������� �������
// files_directory - ����������, � ������� ���������� ������ (in)
// time_moment - ������ � ����������� � ������� ������� ��� ����� �����, ��� �������
// � ��������� ��� ��� ������� ���� � �������������� ������������ (in)
// filename - ��� ���������� ����������� ����� (in)
void write_restart_info( char *files_directory, char *time_moment, char *filename ) {

    FILE *restart_info_file; // ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������
    char string[MAX_STRING_SIZE]; // ��� ������������ ����� �����

    // ����������, ����������� ��� ��������, ��������� � ����� restart_info.dat
    strcpy_s( string, files_directory );
    strcat_s( string, "\\restart_info.dat" );
    if ( ( fopen_s( &restart_info_file, string, "wt" ) ) != 0 ) {
        printf( "\nwrite_restart_info -> can't open file %s for writing\n\n", string );
        exit( EXIT_FAILURE );
    }

    // ������ ���������� � ������� ������� � �� ����� ���������� ����� � �������������� ������������
    fprintf( restart_info_file, "%s %s", time_moment, filename );

    fclose( restart_info_file );

}	



/* ���������� ����� � �������� ��� ���������� ����������� ����������� ����������� ������� 

   params - ��������� � ����������� ��������������� ������������ (in)
   restart_info_file - ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������ (in)

   **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out) */
void read_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, FILE *restart_info_file, double **initial_solution ) {

    // char string[MAX_STRING_SIZE];   /* ��� ���������� ��������� ���������� �� ����� */

}

// ������� ������� ���� � ����������� � �������� �, � ������ ������, ��������� �����������
// ����� � ���������� ����������� ���������������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in) 
// files_directory - ����������, � ������� ���������� ������ (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// file_to_read - ���������� ����� � ��������������� ��� �������� (out)
void try_to_restart( struct ParametersCommon *paramsc, char *files_directory, struct TimeMoment *time_mom, FILE *file_to_read ) {

    FILE *restart_info_file; // ���������� ����� � ����������� � ��������� ��������� ����� � �������������� ������������
    char string[MAX_STRING_SIZE]; // ��� ������������ ����� ����� � ���������� ��������� ����������
    char mom_string[MAX_STRING_SIZE]; // ��� ���������� ���������� � ������� �������, �������� ������������� ��������� ���������� �������������

    // ����������, ����������� ��� ��������, ��������� � ����� restart_info.dat */
    strcpy_s( string, files_directory );
    strcat_s( string, "\\restart_info.dat" );
    if ( ( fopen_s( &restart_info_file, string, "rt" ) ) != 0 ) {
        // ��� ����� � ����������� � ��������, �������� ������ � ��������� �������
        paramsc->isContinue = false;
        return;
    }
    else {
        // ���� ���� � ����������� � ��������, ���������� ������
        paramsc->isContinue = true;
        printf( "\n> File with restart data restart_info.dat is detected. " );
    }

    // ��������� ����������� ����� restart_info.dat
    fscanf_s( restart_info_file, "%s %s", mom_string, MAX_STRING_SIZE, string, MAX_STRING_SIZE );

    // ������ �������� ������� �������
    if ( paramsc->is_output_on_time )
        time_mom->curr_t = atof( mom_string );
    else
        time_mom->steps_num = atoi( mom_string );

    // �������� ����� � ����������� ��������������� �� ����������
    if ( ( fopen_s( &file_to_read, string, "rt" ) ) != 0 ) {
        printf( "\ntry_to_restart -> can't open file %s for writing\n\n", string );
        exit( EXIT_FAILURE );    
    }
    printf( "Data from %s is loading.\n", string );

}
