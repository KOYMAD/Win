// relaxation.cc
// ���������� ��������� � �������� ��� �� ��������� �������

// ���������� ��������� ����������� �� ������:
// Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467.

// ���������� �������� ����������� �� ������:
// ������ �.�. ��������� ������������� ����������� ������� � ������� ����������� ���������� ���� // ������� ���. � 2009. � �. 16, � 2. � �. 62 � 70.
// �������� �������� ������ � \science\utkin\docs\���������� ��������.docx
// ������ �����, ��� �������� � ������:
// Saurel R., Lemetayer O. A multiphase model for compressible flows with interfaces, shocks,
// detonation waves and cavitation // Journal of Fluid Mechanics. � 2001. � V. 431. � P. 239 � 271.

// ��������� ����� ���������: - � ������ (Saurel R., Lemetayer O.) ������������ ������ ������� ��� pi;
//                            - �� ������������� �������� ������ �������� � ���� �������������;
//                            - �������� �� ��������, ����� ��� ����� �������������, � ��� ��������?
//                            - ���������� ��������� ���������� ��������, ��������� ����� ��� ������ ���������� �� ������ ��������, � �������� ����;
//                            - ���������� �������� ��� �������� ���������� �����. 

// (c) ����� �����, 2017
// ������: 29 ����� 2017 �.

#include "relaxation.h"

// �������������� �������� ��� ������� ������� ��������� ���� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� �����������
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ��������������� ���������
// n - �������� ������ �������
// curr_time - ������� ������ �������
// configuration_pressure - ���������������� ��������
void Lr( const struct ParametersCommon* params�, const struct Parameters1d* params1d, struct DebugInfo *debug_info, const double left_ncons[M], const double center_ncons[M],
         const double right_ncons[M], const double dt, const double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure) {
    if ( params�->velocity_relaxation == true ) {
        double solution_ncons_Lrv[M]; // ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ��������� � ��������� ���������� ��������
        Lrv( params�, center_ncons, dt, h, solution_ncons_Lrv, n ); // �������� ��������� ���������� ��������
        if (!params�->pressure_relaxation_compaction)
            Lrp_Ivanov( params�, debug_info, solution_ncons_Lrv, solution_ncons, n ); // �������� ��������� ���������� ��������
        else
            Lrp_compaction(params�, params1d, debug_info, solution_ncons_Lrv, solution_ncons, step_number, n, configuration_pressure);
    }
    else{
        if (!params�->pressure_relaxation_compaction)
            Lrp_Ivanov( params�, debug_info, center_ncons, solution_ncons, n ); // �������� ��������� ���������� ��������
        else{
            if (step_number == 11200)
                write_time_dependent_information_in_a_cell(params�, params1d, curr_time, "relaxation_information", "before_pressure_relaxation_center.dat", solution_ncons, n);
            Lrp_compaction(params�, params1d, debug_info, center_ncons, solution_ncons, step_number, n, configuration_pressure);
            if (step_number == 11200)
                write_time_dependent_information_in_a_cell(params�, params1d, curr_time, "relaxation_information", "after_pressure_relaxation_center.dat", solution_ncons, n);
        }
    }
}

// �������� ���������� �������� ��� ������� ������� ��������� ���� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ��������������� ��������� ��������
// n - �������� ������ ��������
void Lrv( const struct ParametersCommon* paramsc, const double center_ncons[M], const double dt, const double h, double solution_ncons[M], int n ) {

    // �������������� ������� ������� �������� � �������������� ������
    for ( int i = 0; i < n; i++ )
        solution_ncons[i] = center_ncons[i];

    // ������������ ��������� ���
    double t1 = center_ncons[B_DISP] * center_ncons[R_DISP];
    double t2 = ( 1.0 - center_ncons[B_DISP] ) * center_ncons[R_GAS];
    solution_ncons[V_DISP] = ( t1 * center_ncons[V_DISP] + t2 * center_ncons[V_GAS] ) / ( t1 + t2 );
    solution_ncons[V_GAS] = solution_ncons[V_DISP];

    // �������� ���������� ������� ���
    double eg_0 = e_gas( paramsc, center_ncons[P_GAS], center_ncons[R_GAS] ); // ���������� ������� ���� �� ���������
    double ed_0 = e_disp( paramsc, center_ncons[P_DISP], center_ncons[R_DISP] ); // ���������� ������� ���������� ���� �� ���������
    double u_i = calc_u_i( paramsc, center_ncons ); // ����������� ��������� �������� �� ���������
    // ���������� ������� ���� ����� ���������
    double eg_1 = eg_0 + 0.5 * ( solution_ncons[V_GAS] - center_ncons[V_GAS] ) * ( u_i - center_ncons[V_GAS] );
    // ���������� ������� ���������� ���� ����� ���������
    double ed_1 = ed_0 - 0.5 * ( solution_ncons[V_DISP] - center_ncons[V_DISP] ) * ( u_i - center_ncons[V_DISP] );

    // �������� �������� �� ���������� ���������� �������
    solution_ncons[P_GAS] = p_gas( paramsc, eg_1, center_ncons[R_GAS] );
    solution_ncons[P_DISP] = p_disp( paramsc, ed_1, center_ncons[R_DISP] );

}

// �������� ���������� �������� ��� ������� ������� ��������� ���� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� �����������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ��������� � ��������� ���������� ��������
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� � ������� ��������������� ����������
// n - �������� ������ ��������
void Lrp_Ivanov( const struct ParametersCommon* paramsc, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int n ) {

    // �������������� ������� ������� �������� � �������������� ������
    for ( int i = 0; i < n; i++ )
        solution_ncons[i] = center_ncons[i];

    double pi = center_ncons[B_DISP] * center_ncons[P_DISP] + ( 1.0 - center_ncons[B_DISP] ) * center_ncons[P_GAS]; // �������� �� ��������� ������� �� ������� ����

    double g1 = paramsc->g1;
    double g2 = paramsc->g2;
    double p01 = paramsc->p01;
    double p02 = paramsc->p02;

    double g2m = g2 - 1.0; // ���������� �������� ���� ��� �������
    double g1m = g1 - 1.0; // ���������� �������� ���������� ���� ��� �������
    
    // ������ ������������� ����������� ��������� ��� ��������� ������� ���� �� ��������� ���� �� �������

    double a1 = 0.5 * center_ncons[R_GAS] * ( g2 + 1.0 );
    double b1 = - ( center_ncons[P_GAS] + p02 * paramsc->g2 ) - 0.5 * g2m * pi;
    double c1 = - 0.5 * g2m;
    double d1 = - p02 * g2 * center_ncons[R_GAS] - 0.5 * g2m * pi * center_ncons[R_GAS];

    double a2 = 0.5 * center_ncons[R_DISP] * ( g1 + 1.0 );
    double b2 = - ( center_ncons[P_DISP] + p01 * paramsc->g1 ) - 0.5 * g1m * pi;
    double c2 = - 0.5 * g1m;
    double d2 = - p01 * g1 * center_ncons[R_DISP] - 0.5 * g1m * pi * center_ncons[R_DISP];

    double a3 = center_ncons[B_DISP] * center_ncons[R_DISP];
    double b3 = ( 1.0 - center_ncons[B_DISP] ) * center_ncons[R_GAS];
    double c3 = - 1.0;

    // �� ������ ������, ����� �������� �������� �� ��, ��� ������������ ���������

    double G =  - a2 * b1 * c3 - a3 * b2 * c1 + a3 * b1 *c2 - c1 * c3 * d2;
    double H =  - a2 * b1 * b3 + a2 * c3 * d1 - a1 * a3 * b2 - a3 * c2 * d1 - b3 * c1 * d2 - a1 * c3 * d2;
    double I = a2 * b3 * d1 - a1 * b3 * d2;
    
    // ������� ��������� G * y^2 + H * y + I = 0 � ������ ������

    if ( fabs( G ) < paramsc->eps_general ) {
        // ��������� �� ����� ���� �� �������� ����������
        if ( fabs( H ) < paramsc->eps_general ) {
            printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
            vector_debug_print( debug_info );
            exit( EXIT_FAILURE );
        }
        else {
            solution_ncons[R_GAS] = - I / H;
            solution_ncons[P_DISP] = ( d1 - b1 * solution_ncons[R_GAS] ) / ( a1 + c1 * solution_ncons[R_GAS] );
            if ( solution_ncons[P_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[P_GAS] = solution_ncons[P_DISP];
            solution_ncons[R_DISP] = - a3 * solution_ncons[R_GAS] / ( b3 + c3 * solution_ncons[R_GAS] );
            if ( solution_ncons[R_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
            if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 1\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
        }
        debug_info->relaxation_cases[0]++;
    }
    else {
        double D = H * H - 4.0 * G * I;
        if ( D > 0 ) { // ������� - D > epsilon
            // ����� ����������� ���������
            double y1 = 0.5 * ( - H + sqrt( D ) ) / G;
            double y2 = 0.5 * ( - H - sqrt( D ) ) / G;
            bool first_root = false; // ������ ������ �� ������ ������ ������� ������?
            bool second_root = false; // ������ ������ �� ������ ������ ������� ������?
            if ( y1 < 0 ) {
                // ������ ������ �������� �� ��������
                if ( y2 < 0 ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
                else {
                    // ������������ ��� �� ������� �����
                    second_root = true;
                    solution_ncons[R_GAS] = y2;
                    solution_ncons[P_DISP] = ( d1 - b1 * y2 ) / ( a1 + c1 * y2 );
                    if ( solution_ncons[P_DISP] < 0.0 ) {
                        printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                        vector_debug_print( debug_info );
                        exit( EXIT_FAILURE );
                    }
                    solution_ncons[P_GAS] = solution_ncons[P_DISP];
                    solution_ncons[R_DISP] = - a3 * y2 / ( b3 + c3 * y2 );
                    if ( solution_ncons[R_DISP] < 0.0 ) {
                        printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                        vector_debug_print( debug_info );
                        exit( EXIT_FAILURE );
                    }
                    solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
                    if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                        printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 2\n", debug_info->current_cell );
                        vector_debug_print( debug_info );
                        exit( EXIT_FAILURE );
                    }
                    debug_info->relaxation_cases[1]++;
                }
            }
            else {
                // �������� ���������� ��� �� ������� ������
                first_root = true;
                solution_ncons[R_GAS] = y1;
                solution_ncons[P_DISP] = ( d1 - b1 * y1 ) / ( a1 + c1 * y1 );
                if ( solution_ncons[P_DISP] < 0.0 ) {
                    // ������ ������ ��������� �� ���������� �������� �� ��������� �������, ����� ��������� ������
                    first_root = false;
                }
                if ( first_root == true )
                    solution_ncons[P_GAS] = solution_ncons[P_DISP];
                if ( first_root == true ) {
                    solution_ncons[R_DISP] = - a3 * y1 / ( b3 + c3 * y1 );
                    if ( solution_ncons[R_DISP] < 0.0 ) {
                        // ������ ������ ��������� �� ��������� ��������� ���������� ����, ����� ��������� ������
                        first_root = false;
                    }
                }
                if ( first_root == true ) {
                    solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
                    if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                        // ������ ������ ��������� �� ��������� �������� ���� ���������� ����, ����� ��������� ������
                        first_root = false;
                    }
                    debug_info->relaxation_cases[2]++;
                }
            }
            if ( first_root == false && second_root == false ) {
                // ������ ������ �� ������� � ���������� �������, ������� ������
                // ����� �������� y2 < 0 - ���� ���, �� ������
                solution_ncons[R_GAS] = y2;
                solution_ncons[P_DISP] = ( d1 - b1 * y2 ) / ( a1 + c1 * y2 );
                if ( solution_ncons[P_DISP] < 0.0 ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 4\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
                solution_ncons[P_GAS] = solution_ncons[P_DISP];
                solution_ncons[R_DISP] = - a3 * y2 / ( b3 + c3 * y2 );
                if ( solution_ncons[R_DISP] < 0.0 ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 4\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
                solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
                if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                    printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 4\n", debug_info->current_cell );
                    vector_debug_print( debug_info );
                    exit( EXIT_FAILURE );
                }
            }
            debug_info->relaxation_cases[3]++;
        }
        else if ( fabs( D ) < paramsc->eps_general ) {
            // ���� ������
            solution_ncons[R_GAS] = - 0.5 * H / G;
            solution_ncons[P_DISP] = ( d1 - b1 * solution_ncons[R_GAS] ) / ( a1 + c1 * solution_ncons[R_GAS] );
            if ( solution_ncons[P_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 5\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[P_GAS] = solution_ncons[P_DISP];
            solution_ncons[R_DISP] = - a3 * solution_ncons[R_GAS] / ( b3 + c3 * solution_ncons[R_GAS] );
            if ( solution_ncons[R_DISP] < 0.0 ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 5\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            solution_ncons[B_DISP] = center_ncons[B_DISP] * center_ncons[R_DISP] / solution_ncons[R_DISP];
            if ( ( solution_ncons[B_DISP] < 0.0 ) || ( solution_ncons[B_DISP] > 1.0 ) ) {
                printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, case 6\n", debug_info->current_cell );
                vector_debug_print( debug_info );
                exit( EXIT_FAILURE );
            }
            debug_info->relaxation_cases[4]++;
        }
        else {
            printf( "\nLrp_Ivanov -> pressure relaxation fails in cell %d, no roots\n", debug_info->current_cell );
            vector_debug_print( debug_info );
            exit( EXIT_FAILURE );
        }
    }

}

// �������� ���������� �������� ��� ������� ������� ��������� ���� Baer-Nunziato � ������ ��������������� Schwendeman � Saurel
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� c ����������� ��������������� ������������, ��������� 1d ������
// debug_info - ��������� � ���������� �����������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ��������� � ��������� ���������� ��������
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� � ������� ��������������� ����������
// step_number - ����� ������� ������
void Lrp_compaction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int step_number, int n, double *configuration_pressure){

    // �������������� ������� ������� �������� � �������������� ������
    for ( int i = 0; i < n; i++ )
        solution_ncons[i] = center_ncons[i];    

    double g1 = paramsc->g1;
    double g2 = paramsc->g2;
    double p01 = paramsc->p01;
    double p02 = paramsc->p02;

    double z; // ��������� ���������� ���� - ������� ����������
    double x, y; // �������� ������� ���� � ��������� ������� ����

    double dx, dy; // ����������� x � y �� z

    double beta; // ���������������� ��������
    double f, df; // ���������������� �������� ���� � ����������� �� z. ����������� �� ������ ����

    double A, a, b; // ��������������� ������������

    a = solution_ncons[B_DISP] * solution_ncons[R_DISP];
    b = (1.0 - solution_ncons[B_DISP]) * solution_ncons[R_GAS];

    double common_part; // ����� ����� (B / tau)^(1/n) � ������� SAUREL_COMPACTION

    double p0, p0_disp, b0, b0_disp, r0_disp;

    // ����������� �����, � ������� ������ ���������
    int number_of_block = 0;

    if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION){
        for (int i = 0; i < params1d->ic_blocks_number; i++){
	    if (step_number <= params1d->cell_end[i] && step_number >= params1d->cell_begin[i])
	        number_of_block = i;
        }

        p0 = params1d->block_values[number_of_block][P_GAS];
        p0_disp = params1d->block_values[number_of_block][P_DISP];
        b0 = 1.0 - params1d->block_values[number_of_block][B_DISP];
        b0_disp = params1d->block_values[number_of_block][B_DISP];
        r0_disp = params1d->block_values[number_of_block][R_DISP];
  
        if (solution_ncons[B_DISP] < paramsc->volume_fraction_compaction)
            A = 0.0;
        else
            A = -(p0 - p0_disp) / b0_disp / r0_disp * (2.0 - b0_disp) * (2.0 - b0_disp) / log(1.0 - b0_disp); // ��������� � ��������� ��� ����������������� ��������

        beta = A * a * log(1.0 - solution_ncons[B_DISP]) / pow(2.0 - solution_ncons[B_DISP], 2.0);
    }
    else{
        
        common_part = (1.0 - solution_ncons[B_DISP]) * log( (1.0 - solution_ncons[B_DISP]) / (1.0 - paramsc->volume_fraction_compaction) ) + solution_ncons[B_DISP] - paramsc->volume_fraction_compaction;
        if (common_part < 0.0 && fabs(common_part) < paramsc->eps_relax_compaction){
            common_part = 0.0;
        }
        if (solution_ncons[B_DISP] > paramsc->volume_fraction_compaction)
            beta = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - solution_ncons[B_DISP]) / (1.0 - paramsc->volume_fraction_compaction) ) * pow (common_part , paramsc->n_parameter - 1.0);
        else
            beta = 0.0;
    }

    double fi, dfi;

    z = solution_ncons[R_DISP];
 //   printf("Before relaxation procedure\n");
//    for (int i = 0; i < n; i++)
//        printf("solution_ncons[%d] = %lf \n", i, solution_ncons[i]);
//    printf("\n");
 //   if (solution_ncons[P_DISP] - beta > solution_ncons[P_GAS]){

    do{
		p0 = params1d->block_values[number_of_block][P_GAS];
		p0_disp = params1d->block_values[number_of_block][P_DISP];
		b0 = 1.0 - params1d->block_values[number_of_block][B_DISP];
		b0_disp = params1d->block_values[number_of_block][B_DISP];
		r0_disp = params1d->block_values[number_of_block][R_DISP];
        if (a / z < paramsc->volume_fraction_compaction)
            A = 0.0;
        else
            A = -(p0 - p0_disp) / b0_disp / r0_disp * (2.0 - b0_disp) * (2.0 - b0_disp) / log(1.0 - b0_disp); // ��������� � ��������� ��� ����������������� ��������

        y = b * z / (z - a);
        x = ( 1.0 / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) + (solution_ncons[P_GAS] + g2 * p02) / (g2 - 1.0) / solution_ncons[R_GAS] + solution_ncons[P_GAS] / 2.0 / solution_ncons[R_GAS] ) / 
            ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 0.5 / solution_ncons[R_GAS] );
        
        if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION){
            f = A * a * log(1.0 - a / z) / pow ( 2.0 - a / z , 2.0 ); 
            df = A * a * a / z / z / (2 - a / z) / (2 - a / z) * ( 1.0 / (1.0 - a / z) - 2.0 * log (1.0 - a / z) / (2.0 - a / z) );// (2 - a/z)
        }
        else{
            if (a / z > paramsc->volume_fraction_compaction){
                common_part = ((1.0 - a / z) * log ( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ) + a / z - paramsc->volume_fraction_compaction);
                if (common_part < 0.0 && fabs(common_part) < paramsc->eps_relax_compaction){
                    common_part = 0.0;
                }
                f = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ) * pow ( common_part , paramsc->n_parameter - 1.0);
                df = - a * a * paramsc->tau_parameter * paramsc->n_parameter / z / z * ( z / (z - a) * pow(common_part, paramsc->n_parameter - 1.0) + 
                    (paramsc->n_parameter - 1.0) * pow(log( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ), 2.0) * pow(common_part, paramsc->n_parameter - 2.0) );

                if (f == 0.0)
                    df = 0.0;

            }
            else{
                f = 0;
                df = 0;
            }
        }
        fi = (x + f + g1 * p01) / (g1 - 1.0) / z - (solution_ncons[P_DISP] + g1 * p01) / (g1 - 1.0) / solution_ncons[R_DISP] + 0.5 * (x + f + solution_ncons[P_GAS] + beta) * (1.0 / z - 1.0 / solution_ncons[R_DISP]);

        dy = - b * a / (z - a) / (z - a);
        dx = (- dy / y / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) * ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 1.0 / 2.0 / solution_ncons[R_GAS] ) + 
            (1.0 / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) + (solution_ncons[P_GAS] + g2 * p02) / (g2 - 1.0) / solution_ncons[R_GAS] + solution_ncons[P_GAS] / 2.0 / solution_ncons[R_GAS] ) * (- dy / y / y * (g2 + 1.0) / 2.0 / (g2 - 1.0) ) ) /
            pow ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 1.0 / 2.0 / solution_ncons[R_GAS] , 2.0 );
        
        dfi = (dx + df) / (g1 - 1.0) / z - 1.0 / z / z * (x + f + g1 * p01) / (g1 - 1.0) - 1.0 / 2.0 / z / z * (x + f + solution_ncons[P_GAS] + beta) + 0.5 * (dx + df) * (1.0 / z - 1.0 / solution_ncons[R_DISP]);

        z = z - fi / dfi;

    }while(fabs(fi/dfi / solution_ncons[R_DISP]) > 10 * paramsc->eps_general);

    y = b * z / (z - a);
    x = ( 1.0 / y * ( -g2 * p02 / (g2 - 1.0) - solution_ncons[P_GAS] / 2.0 ) + (solution_ncons[P_GAS] + g2 * p02) / (g2 - 1.0) / solution_ncons[R_GAS] + solution_ncons[P_GAS] / 2.0 / solution_ncons[R_GAS] ) / 
            ( (g2 + 1.0) / 2.0 / (g2 - 1.0) / y - 0.5 / solution_ncons[R_GAS] );

    double v_ncons_new[M];
    double beta_new;
    v_ncons_new[R_DISP] = z;
    v_ncons_new[R_GAS] = y;
    v_ncons_new[P_GAS] = x;
    v_ncons_new[B_DISP] = solution_ncons[B_DISP] * solution_ncons[R_DISP] / v_ncons_new[R_DISP];
    if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION)
        beta_new = A * a * log(1.0 - v_ncons_new[B_DISP]) / pow(2.0 - v_ncons_new[B_DISP], 2.0);
    else{
        if (v_ncons_new[B_DISP] > paramsc->volume_fraction_compaction)
            beta_new = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - a / v_ncons_new[R_DISP]) / (1.0 - paramsc->volume_fraction_compaction) ) * pow ( common_part , paramsc->n_parameter - 1.0);
        else
            beta_new = 0.0;
    }
    v_ncons_new[P_DISP] = beta_new + v_ncons_new[P_GAS];
    v_ncons_new[V_DISP] = solution_ncons[V_DISP];
    v_ncons_new[V_GAS] = solution_ncons[V_GAS];
    
   /* double check_f;
    if (paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION)
        check_f = A * a * log(1.0 - a / z) / pow ( 2.0 - a / z , 2.0 );
    else
        if (solution_ncons[B_DISP] > paramsc->volume_fraction_compaction )
            check_f = - a * paramsc->tau_parameter * paramsc->n_parameter * log( (1.0 - a / z) / (1.0 - paramsc->volume_fraction_compaction) ) * pow ( common_part , paramsc->n_parameter - 1.0);
        else
            check_f = 0.0;
    double check_fi = (x + check_f + g1 * p01) / (g1 - 1.0) / z - (solution_ncons[P_DISP] + g1 * p01) / (g1 - 1.0) / solution_ncons[R_DISP] + 0.5 * (x + check_f + solution_ncons[P_GAS] + beta) * (1.0 / z - 1.0 / solution_ncons[R_DISP]);
*/
   // if (v_ncons_new[P_DISP] <= 0.0 )
        //printf("current step is %d, check_fi = %lf, v_ncons_new[P_DISP] = %lf\n", step_number, check_fi, v_ncons_new[P_DISP]);

   // if (v_ncons_new[P_GAS] <= 0.0 )
        //printf("current step is %d, check_fi = %lf, v_ncons_new[P_DISP] = %lf\n", step_number, check_fi, v_ncons_new[P_GAS]);

    //if (check_fi > 1.e-5)
        //printf("current step is %d, check_fi = %lf\n", step_number, check_fi);

    solution_ncons[B_DISP] = v_ncons_new[B_DISP];          
    solution_ncons[R_DISP] = v_ncons_new[R_DISP];  
    solution_ncons[V_DISP] = v_ncons_new[V_DISP];  
    solution_ncons[P_DISP] = v_ncons_new[P_DISP];  
    solution_ncons[R_GAS] = v_ncons_new[R_GAS];  
    solution_ncons[V_GAS] = v_ncons_new[V_GAS];  
    solution_ncons[P_GAS] = v_ncons_new[P_GAS];  
    *configuration_pressure = beta_new;

  //  printf("After relaxation procedure\n");
  //  for (int i = 0; i < n; i++)
 //       printf("solution_ncons[%d] = %lf \n", i, solution_ncons[i]);
 //   printf("\n");
}