// godunov.cc
// ����� �������� ������� ������� ��� ���������� ������� ��������� �����-��������.
// Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
// two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526.
// (c) ����� �����, 2013
// ������: 4 ���� 2013 �.

#include "godunov_bn.h"

// ����� �������� �������������� ���������� ������� ��������� �����-��������,
// ������ ������ ���� �� ������� � ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// left_params[M] - ������ ����������� ���������� � ������ ����� �� �������������� (in)
// center_params[M] - ������ ����������� ���������� � �������������� ������ (in)
// right_params[M] - ������ ����������� ���������� � ������ ������ �� �������������� (in)
// slopes_left - ������ �������� � ������ ����� �� �������������� (in)
// slopes_center - ������ �������� � �������������� ������ (in)
// slopes_right - ������ �������� � ������ ������ �� �������������� (in)
// dt - ��������� ��� (in)
// h - ���������������� ��� (in)
// solution[M] - ������ ����������� ���������� � �������������� ������ �� ��������� ���� (out)
// n - �������� ������ ��������
// configuration_pressure - ���������������� ��������
void godunov( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_params[M], double center_params[M],
              double right_params[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
              double dt, double h, double solution[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, double curr_time, double *configuration_pressure ) {

    double left_godunov_flux[M];            /* ������ "�����" ����� ����� ����� ������ */
    double right_godunov_flux[M];           /* ������ "�����" ����� ������ ����� ������ */
    double center_cons_params[M];           /* ������ �������������� ���������� � �������������� ������ �� ������� ���� */
    double solution_cons[M];                /* ������ �������������� ���������� � �������������� ������ �� ����� ���� */
    ReturnCodes left_flux_code, right_flux_code, pressure_code; /* ��� ��������� ������� ������� ������� � ��������
                                                                   "�������������� ����������" � ����������� */
    Disp_phase_cases solver_part_left;      /* ���������, ������� ����������, ����� �� ������� ������ ������� ������������
                                               ��� ������� "������" ����� ����� ����� ������ */
    Disp_phase_cases solver_part_right;     /* ���������, ������� ����������, ����� �� ������� ������ ������� ������������
                                               ��� ������� "������" ����� ������ ����� ������ */

    // ������������� �������� �������
    double left_minus_ncons[M];
    double left_minus_cons[M];
    double left_plus_ncons[M];
    double left_plus_cons[M];
    double right_minus_ncons[M];
    double right_minus_cons[M];
    double right_plus_ncons[M];
    double right_plus_cons[M];
    for ( int j = 0; j < n; j++ ) {
        // ������������� �������� �� ����� ������� �����
        convert_noncons_to_cons( paramsc, left_params, left_minus_cons, 0 );
        left_minus_cons[j] += 0.5 * h * slopes_left[j];
        convert_cons_to_noncons( paramsc, left_minus_cons, left_minus_ncons, 0 );
        // ������������� �������� �� ����� ������� ������
        convert_noncons_to_cons( paramsc, center_params, left_plus_cons, 0 );
        left_plus_cons[j] -= 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, left_plus_cons, left_plus_ncons, 0 );
        // ������������� �������� �� ������ ������� �����
        convert_noncons_to_cons( paramsc, center_params, right_minus_cons, 0 );
        right_minus_cons[j] += 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, right_minus_cons, right_minus_ncons, 0 );
        // ������������� �������� �� ������ ������� ������
        convert_noncons_to_cons( paramsc, right_params, right_plus_cons, 0 );
        right_plus_cons[j] -= 0.5 * h * slopes_right[j];
        convert_cons_to_noncons( paramsc, right_plus_cons, right_plus_ncons, 0 );
    }

    // ������ ����������� ����� ���������� �������� ���� ���������� ���� ����� �� ������� ������ � � �������
    solver_part_left = what_case( paramsc, left_minus_ncons[B_DISP], left_plus_ncons[B_DISP] );
    
    // ������ ����������� ����� ���������� �������� ���� ���������� ���� � ������� ������ � ������ �� ���
    solver_part_right = what_case( paramsc, right_minus_ncons[B_DISP], right_plus_ncons[B_DISP] );

    if ( solver_part_left == NO_GRAD && solver_part_right == NO_GRAD ) {
        // ������ ����������� ��� - ������� ������� � �������������� ������������ �������
        full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution, n );
        debug_info->godunov_cases[0]++;
        debug_info->godunov_cases[0]++;
    }
    else { 
        // ���������� ������ ������, ���� ���� ���������� �������� ���� � ������ ����
        
        // ������ ������� "������" ����� ����� �����
        
        // ���������� ����� ����� ��������� ��� ������� ����
        debug_info->neighbour_cell = debug_info->current_cell - 1;
        for ( int i = 0; i < n; i++ )
            debug_info->current_cell_vncons[i] = left_plus_ncons[i];
        for ( int i = 0; i < n; i++ )
            debug_info->neighbour_cell_vncons[i] = left_minus_ncons[i];
        
        left_flux_code = godunov_flux( paramsc, debug_info, left_minus_ncons, left_plus_ncons,
            solver_part_left, LEFT, center_params[B_DISP], left_godunov_flux, n );
        if ( left_flux_code == GODUNOV_FAILS ) {
            // ����� �������� �� �������� - ������� ��������� ��������� � ����� �������
            full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution, n );
            return;
        }

        // ������ ������� "������" ����� ������ �����
        
        // ���������� ����� ����� ��������� ��� ������� ����, ���������� ������ � �������� ������
        debug_info->neighbour_cell = debug_info->current_cell + 1;
        for ( int i = 0; i < n; i++ )
            debug_info->current_cell_vncons[i] = right_minus_ncons[i];
        for ( int i = 0; i < n; i++ )
            debug_info->neighbour_cell_vncons[i] = right_plus_ncons[i];
        
        right_flux_code = godunov_flux( paramsc, debug_info, right_minus_ncons, right_plus_ncons,
            solver_part_right, RIGHT, center_params[B_DISP], right_godunov_flux, n );
        if ( right_flux_code == GODUNOV_FAILS ) {
            // ����� �������� �� �������� - ������� ��������� ��������� � ����� �������
            full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution , n);
            return;
        }

        // ���������� ���������� � ������
        convert_noncons_to_cons( paramsc, center_params, center_cons_params, 0 );
        for ( int i = 0; i < n; i++ ) {
            solution_cons[i] = center_cons_params[i] - dt * ( right_godunov_flux[i] - left_godunov_flux[i] ) / h;
        }
        pressure_code = convert_cons_to_noncons( paramsc, solution_cons, solution, 0 );
        if ( pressure_code == NEGATIVE_PRESSURE ) {
            /* � ���������� ������� ���������� ������������� �������� � ���������� ���� - ������� ��������� ��������� */
            full_decouple_case_sol( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, dt, h, solution, n );
        }
    }

    // ��� ������������� ������ ���������� ��������� � ��������
    double solution_ncons_Lh[M]; // ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
    for ( int i = 0; i < n; i++ )
        solution_ncons_Lh[i] = solution[i];
    if ( paramsc->pressure_relaxation == true && is_pressure_relaxation_after_this_step == true )
        Lr( paramsc, params1d, debug_info, left_params, solution_ncons_Lh, right_params, dt, h, solution, step_number, n, curr_time, configuration_pressure ); // �������� ��������������� ���������

}

// ������ �������� �������� �������� ���� ���������� ���� ����� � ������ �� �������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// left_beta - �������� �������� ���� ���������� ���� ����� �� ������� (in)
// right_beta - �������� �������� ���� ���������� ���� ������ �� ������� (in)
// ���������� ���� �� �������� ������������ Disp_phase_cases
Disp_phase_cases what_case( struct ParametersCommon *paramsc, double left_beta, double right_beta ) {
    
    if ( fabs( left_beta - right_beta ) <= paramsc->eps_decouple ) {
        // ������� �������� �������� ���� ����� � ������ �� ������� �������������� - ������� ��������� ������������,
        // � ���������� ��������� �������� ���� ���������� ����; ��������������, ��� ������ ������ �������� ��������
        // ���������� ���������� ���� �� ��� ������� �� �������, �� ���� ��������� params->eps_decouple �����������
        // ������ ��������� params->eps_disp_abs
        return NO_GRAD;
    }
    else { 
        // ������� �������� �������� ���� ����� � ������ �� ������� ������������ - ������������� ��������� � �����������������
        // ������� � ������ �����
        if ( left_beta > paramsc->eps_disp_abs ) {
            // ���������� ���� ����� �� ������� ������������
            if ( right_beta > paramsc->eps_disp_abs ) {
                // ���������� ���� ������ �� ������� ������������
                return BOTH_GRAD;
            }
            else {
                // ���������� ���� ����� �� ������� ������������, � ������ �����������
                return LEFT_ONLY_GRAD;
            }
        }
        else {
            // ���������� ���� ����� �� ������� �����������
            if ( right_beta > paramsc->eps_disp_abs ) {
                // ���������� ���� ������ �� ������� ������������
                return RIGHT_ONLY_GRAD;
            }
            else {
                // ���������� ���� ����������� �� ��� ������� �� ������� - �������� ���� �� ������, ��� ������
                // ���������� �������� �������� �������� ���� ����� � ������ �� �������, �� �����������������
                return NO_GRAD;
            }
        }
    }

}

// ��������� ������� �� ��������� ���� �� ������� ��� ������ ���������� ����������������� ����� � ������ ����� - 
// ������ ����������� ���
// params� - ��������� � ��������� ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// center_params[M] - ������ ����������� ���������� � �������������� ������ (in)
// left_minus_ncons[M] - ������������������ ������ ����������� ���������� �� ������ ������ ������ ����� �� �������������� (in)
// left_plus_ncons[M] - ������������������ ������ ����������� ���������� �� ����� ����� �������������� ������ (in)
// right_minus_ncons[M] - ������������������ ������ ����������� ���������� �� ������ ����� �������������� ������ (in)
// right_plus_ncons[M] - ������������������ ������ ����������� ���������� �� ����� ����� ������ ������ �� �������������� (in)
// dt - ��������� ��� (in)
// h - ���������������� ��� (in)
// solution[M] - ������ ����������� ���������� � �������������� ������ �� ��������� ���� (out)
// n (in)
void full_decouple_case_sol( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double center_params[M], double left_minus_ncons[M], 
                             double left_plus_ncons[M], double right_minus_ncons[M], double right_plus_ncons[M], double dt, double h, double solution[M], int n ) {

    double v_ncons_res_solid_left[M];
    double v_ncons_res_solid_right[M];

    // ��������� ���������� ������� solution, ��������������� ������� ����
    godunov_classical_one_phase( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, GAS_PHASE, dt, h, 
        v_ncons_res_solid_left, v_ncons_res_solid_right, solution, n );
    // ��������� ���������� ������� solution, ��������������� ���������� ����
    godunov_classical_one_phase( paramsc, debug_info, center_params, left_minus_ncons, left_plus_ncons, right_minus_ncons, right_plus_ncons, DISPERSED_PHASE, dt, h, 
        v_ncons_res_solid_left, v_ncons_res_solid_right, solution, n );
    // �������� ���� ��������� ��� ���������
    solution[B_DISP] = center_params[B_DISP];

}

// ���� � ������ ���������� ����������� ���������� ����, �� ���������� � ��� ������� ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// solution[M] - ������ ����������� ���������� � �������������� ������ (in/out)
void set_background_state( struct ParametersCommon *paramsc, double solution[M] ) {

    if ( solution[B_DISP] < paramsc->eps_disp_abs ) {
        solution[R_DISP] = paramsc->background_density;
        solution[V_DISP] = paramsc->background_velocity;
        solution[P_DISP] = paramsc->background_pressure;
    }

}

// ������ ������� "������" ����� ����� ������ ������� ��������
// params - ��������� � ��������� ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
// right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
// solver_part - ���������, ������� ����������, ����� �� ������� ������ ������� ������������ ��� ������� "������" (in)
// edge - ������������� �����, ����� ������� ��������� "�����" - LEFT ��� RIGHT (in)
// curr_cell_beta - �������� �������� ���� � �������������� ������ (in)
// godunov_flux[M] - ������ "�����" �������� (out)
// n - �������� ������ ������� left_params, right_params � godunov_dlux (in)
// ����������: SUCCESS          � ������ ������
//             GODUNOV_FAILS    � ������ ������������� ��������� ������� ������ ������
ReturnCodes godunov_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                          Disp_phase_cases solver_part, Direction edge, double curr_cell_beta, double godunov_flux[M], int n ) {

    double cons_flux[M]; // �������������� ������������ ������� "������"
    double v_ncons_res[M]; // ������-������� ������ � ������� �������, �� ������������, �� ����� � ��������
                           // ��������� ��� godunov_cons_flux
    double cont_ncons[M]; //����� ��� ���������
    double s = 0.0; // �������� ������������� ���������� x/t
    double ncons_flux[M]; // ���������������� ������������ ������� "������"
    double v_solid_cont; // �������� ����������� ������� � ���������� ����
    double p_solid_cont_l, p_solid_cont_r; // �������� ����� � ������ �� ����������� ������� � ���������� ����
    int i; // ������ ���������� �������
    ReturnCodes return_code; // ��� �������� ������� ���������� ���������� �� �������
    
    if ( solver_part == NO_GRAD ) {
        /* ��������� ��������� ������� ��������� � ���������� ����������� ������� �������, �� �����
           ������������ ������ "�����", ��������� � �������� �������� ���� �������� � ������� ������ */
        full_decouple_case_flux( paramsc, debug_info, left_params, right_params, curr_cell_beta, godunov_flux );
        debug_info->godunov_cases[0]++;
    }
    else {
        /* ������� ��������� �� ������������, ����� ������ */

        /* ������ �������������� ������������ ������� "������" */
        return_code = godunov_cons_flux( paramsc, debug_info, left_params, right_params, s, solver_part,
            cons_flux, v_ncons_res, &v_solid_cont, &p_solid_cont_l, &p_solid_cont_r, n, cont_ncons );
        if ( return_code == GODUNOV_FAILS ) {
            /* ������ ������ ������ �� ������� */
            return return_code;
        }

        /* ������ ���������������� ������������ ������� "������" */
        for ( i = 0; i < n; i++ )
            ncons_flux[i] = 0.0;   /* ������������� */
        if ( v_solid_cont > 0.0 && edge == LEFT )
            calc_ncons_term( v_solid_cont, left_params[B_DISP], right_params[B_DISP], p_solid_cont_l, p_solid_cont_r, ncons_flux );
        if ( v_solid_cont < 0.0 && edge == RIGHT )
            calc_ncons_term( v_solid_cont, left_params[B_DISP], right_params[B_DISP], p_solid_cont_l, p_solid_cont_r, ncons_flux );


        /* ������ ������� "������" */
        if ( edge == LEFT ) {
            for ( i = 0; i < n; i++ )
                godunov_flux[i] = cons_flux[i] + ncons_flux[i];
        }
        else {
            for ( i = 0; i < n; i++ )
                godunov_flux[i] = cons_flux[i] - ncons_flux[i];
        }
    }

    return SUCCESS;

}

/* �������� ��������� ������� ��������� � ���������� ����������� ������� �������, �� �����
   ������������ ������ "�����", ��������� � �������� �������� ���� �������� � ������� ������.

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   debug_info - ��������� � ���������� ����������� (in)
   left_params_full[M] - ������ ������ ����������� ���������� ����� �� ������� (in)
   right_params_full[M] - ������ ������ ����������� ���������� ������ �� ������� (in)
   beta - �������� �������� ���� � �������������� ������ (in)
   
   flux_full[M] - ������ "�����" (out) */
void full_decouple_case_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params_full[M], double right_params_full[M],
                              double beta, double flux_full[M] ) {

    double flux_reduced[M_REDUCTION];           /* ����������� ���������� ������ ������ */
    double left_params_reduced[M_REDUCTION];    /* ���������� ������ ����������� ���������� � ������ ����� �� ���������������� ������� */
    double right_params_reduced[M_REDUCTION];   /* ���������� ������ ����������� ���������� � ������ ������ �� ���������������� ������� */
    double v_ncons_res[M_REDUCTION];            /* ���������� ������ ����������� ����������, ���������� � ���������� ������� ������ � ������� ������� */
    double cont_red[M_REDUCTION];
    /* ������ ������� ������ � ������� ���� */
    convert_full_to_reduced( left_params_full, GAS_PHASE, left_params_reduced );
    convert_full_to_reduced( right_params_full, GAS_PHASE, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, GAS_PHASE, v_ncons_res, flux_reduced,cont_red );
    convert_reduced_to_full( flux_reduced, GAS_PHASE, flux_full );

    /* ������ ������� ������ � ���������� ���� */
    convert_full_to_reduced( left_params_full, DISPERSED_PHASE, left_params_reduced );
    convert_full_to_reduced( right_params_full, DISPERSED_PHASE, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, DISPERSED_PHASE, v_ncons_res, flux_reduced,cont_red );
    convert_reduced_to_full( flux_reduced, DISPERSED_PHASE, flux_full );

    double beta_new = full_decouple_case_flux_volume_fraction( 0, left_params_full[B_DISP], beta, right_params_full[B_DISP], v_ncons_res[V]);

    if (beta_new != beta){
        printf("here\n");
    }

    /* ���� �������� ���� */
    flux_full[B_DISP] = 0.0;
    flux_full[R_DISP] *= beta_new;
    flux_full[V_DISP] *= beta_new;
    flux_full[P_DISP] *= beta_new;
    flux_full[R_GAS] *= 1.0 - beta_new;
    flux_full[V_GAS] *= 1.0 - beta_new;
    flux_full[P_GAS] *= 1.0 - beta_new;
}

/*  ������� ������� �������������� ������������ "������" ������� ��������
    ��� ���������� ������� ��������� �����-��������

    paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
    debug_info - ��������� � ���������� ����������� (in)
    left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
    right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
    s - �������� ������������� ���������� x/t (in)
    solver_part - ���������, ������� ����������, ����� �� ������� ������ ������� ������������ ��� ������� "������" (in)

    flux[M] - �������������� ������������ ������� "������" (out)
    v_ncons_res[M] - ������-������� ������ � ������� �������, ����� � �������� ���������� ��� ���������� ������� ������� (out)
    v_solid_cont - �������� ����������� ������� � ���������� ���� (out)
    p_solid_cont_l - �������� � ���������� ���� ����� �� ����������� ������� � ���������� ���� (out)
    p_solid_cont_r - �������� � ���������� ���� ������ �� ����������� ������� � ���������� ���� (out)
    n - �������� ������ ������� left_params, right_params � godunov_dlux (in)
    ����������: SUCCEESS        �������� ���������� ��������
                GODUNOV_FAILS   �������� �� �������, ������� ������ ������ ��������� �� ������� */
ReturnCodes godunov_cons_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                               double s, Disp_phase_cases solver_part, double flux[M], double v_ncons_res[M], double *v_solid_cont,
                               double *p_solid_cont_l, double *p_solid_cont_r, int n, double cont_ncons[M] ) {

    double c1l, c2l;                    /* �������� ����� � ���������� � ������� ����� ����� �� ������� */
    double c1r, c2r;                    /* �������� ����� � ���������� � ������� ����� ������ �� ������� */
    double p_cont_disp, p_cont_gas;     /* �������� �� ���������� ������� � ������� � ���������� ����� ��� ����� �������
                                           �������� ���� ���������� ���� */
    double v_cont_disp, v_cont_gas;     /* �������� �� ���������� ������� � ������� � ���������� ����� ���
                                           ����� ������� �������� ���� ���������� ���� */
    double solid_discontinuity_pressures[K_GENERAL_CASE];   /* [GAS_LEFT] - �������� ����� � ����, [GAS_RIGHT] - �������� ������ � ����,
                                                               [DISP_LEFT] - �������� ����� � ���������� ����,
                                                               [DISP_RIGHT] - �������� ������ � ���������� ���� */
    double sound_velocities[K_GENERAL_CASE];    /* ������ �� ���������� ����� ��� �� ������ ������� �� ������� ���
                                                   �������� � ������� */
    double v1, v2;                      /* �������� ���������� �������� � ���������� � ������� ����� */
    double v21, v22;                    /* �������� ���� ����� � ������ �� ����������� ������� ���������� ���� */
    int i_comp;                         /* ������� ��������� �������-������� */
    int return_code;                    /* ��� �������� ������� ���������� ���������� �� ������� */
    
    /* 0. ������������ �������� ����� � ����� ����� ����� � ������ �� ������� ��� ������������ ������������� */
    calc_sound_velocity( paramsc, left_params, &c1l, &c2l );
    calc_sound_velocity( paramsc, right_params, &c1r, &c2r );
    fill_sound_velocities( c2l, c2r, c1l, c1r, sound_velocities );

    /* 1. ���������� ���������� ����������� - �������� � �������� �� ���������� ������� � ���������� ���� ��� �����
          ��� ��� ������������� �������� */
    return_code = calc_p_v_initial_guess( paramsc, debug_info, left_params, right_params, sound_velocities,
        solver_part, &p_cont_gas, &v_cont_gas, &p_cont_disp, &v_cont_disp );
    if ( return_code == GODUNOV_FAILS )
        /* ��������� ����������� ��� ������������� �������� ��������� �� ������� */
        return GODUNOV_FAILS;
    
    /* 2. ������ �������� � ��������� ����� ��� �� ���������� ������� � ���������� ���� */
    /* ������������� �������� � ��������� �� ������ ������� �� ������� � ���������� ���� */
    init_solid_discontinuity_pressures( p_cont_gas, p_cont_disp, solid_discontinuity_pressures );
    v1 = v_cont_disp;
    v2 = v_cont_gas;
    v21 = v_cont_gas;
    v22 = v_cont_gas;
    return_code = calc_p_v( paramsc, debug_info, left_params, right_params, sound_velocities,
        solver_part, solid_discontinuity_pressures, &v1, &v2, &v21, &v22 );
    if ( return_code == GODUNOV_FAILS )
        /* �������� � �������� � ������ ������ ����� �� ������� */
        return GODUNOV_FAILS;
    cont_ncons[P_GAS] = solid_discontinuity_pressures[GAS_RIGHT];
    cont_ncons[P_DISP] = solid_discontinuity_pressures[DISP_RIGHT];
    cont_ncons[V_GAS] = v_cont_gas;
    cont_ncons[V_DISP] = v_cont_disp;
    /* 3. ����� ������� � ���������� ���� */
    /* ������������� �������-������� ������ */
    for ( i_comp = 0; i_comp < n; i_comp++ )
        v_ncons_res[i_comp] = 0.0;
    sample_solid_solution( paramsc, left_params, right_params, c1l, c1r, solid_discontinuity_pressures[DISP_LEFT],
        solid_discontinuity_pressures[DISP_RIGHT], v1, s, v_ncons_res, cont_ncons );

    /* 4. ����� ������� � ������� ���� */
    sample_gas_solution( paramsc, left_params, right_params, c2l, c2r, solid_discontinuity_pressures[GAS_LEFT],
        solid_discontinuity_pressures[GAS_RIGHT], v1, v2, v21, v22, s, v_ncons_res, cont_ncons );

    // 5. ������ �������������� ������������ "������"
	array1D tmp_flux( M );
    diff_flux_ncons( paramsc, v_ncons_res, &tmp_flux, 0 );
	for ( int i = 0; i < n; i++ )
		flux[i] = tmp_flux[i];

    /* 6. ���������� ���������� ��� ������������� ���������������� ������������ ������� "������" */
    *v_solid_cont = v1;
    *p_solid_cont_l = solid_discontinuity_pressures[DISP_LEFT];
    *p_solid_cont_r = solid_discontinuity_pressures[DISP_RIGHT];

    return SUCCESS;

}

// ���������� ���������� ����������� - �������� � �������� �� ���������� ������� � ���������� ���� ��� ����� ���
// ��� ������������� ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
// right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
// c - ������ �� ���������� ����� ��� �� ������ ������� �� ������� (in)
// p_cont_gas - ��������� ����������� ��� �������� � ���� �� ���������� ������� � ���������� ���� (out)
// v_cont_gas - ��������� ����������� ��� �������� ���� �� ���������� ������� � ���������� ���� (out)
// p_cont_disp - ��������� ����������� ��� �������� � ���������� ���� �� ���������� ������� � ���������� ���� (out)
// v_cont_disp - ��������� ����������� ��� �������� ���������� ���� �� ���������� ������� � ���������� ���� (out)
// solver_part - ���������, ������� ����������, ����� �� ������������ �������� ���� ���������� ���� ����������� (in)
// ����������: SUCCEESS         ��������� ����������� ��������� �������
//             GODUNOV_FAILS    ��������� ����������� ��������� �� �������
ReturnCodes calc_p_v_initial_guess( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M],
                                    double right_params[M], double c[K_GENERAL_CASE], Disp_phase_cases solver_part,
                                    double *p_cont_gas, double *v_cont_gas, double *p_cont_disp, double *v_cont_disp ) {

    double F, DF;                                   // �������, ������������ �������� ����� �� ���������� �������, � �� �����������
    double c1l = c[DISP_LEFT], c1r = c[DISP_RIGHT]; // �������� ����� � ���������� ���� ����� � ������ �� ��������������� �������
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   // �������� ����� � ������� ���� ����� � ������ �� ��������������� �������
    int return_code;                                // ��� ��������

    // ������ ������� �������� � �������� �� ���������� ������� � ������� ���� �� ������������ ��������� ������ ��������
    return_code = calc_contact_pressure_velocity( paramsc, debug_info, left_params, right_params, M, c2l, c2r, GAS_PHASE,
        p_cont_gas, v_cont_gas );
    if ( return_code == GODUNOV_FAILS )
        // �� ������� ��������� ��������� ����������� ��� ������� ������ ������
        return GODUNOV_FAILS;
    
    // ��� ���������� ���� - � ����������� �� ������������ �������� ���� ��������� ����
    switch ( solver_part ) {
        case BOTH_GRAD:
            // ���������� ���� ������������ � ������, � ����� �� �������. ������� �������� � �������� �� ����������
            // ������� � ���������� ���� �� ������������ ��������� ������ �������� ��� ������ ����������� ��������� ���������
            // Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
            // two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 499, formula (24).
            return_code = calc_contact_pressure_velocity( paramsc, debug_info, left_params, right_params, M, c1l, c1r,
                DISPERSED_PHASE, p_cont_disp, v_cont_disp );
            if ( return_code == GODUNOV_FAILS )
                // �� ������� ��������� ��������� ����������� ��� ������� ������ ������
                return GODUNOV_FAILS;
            break;
        case LEFT_ONLY_GRAD:
            // ��������� ���� ������������ ����� �� ������� � ����������� ������ �� �������, ������ ������� ������������, ������ ������
            // Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
            // two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 501 + �����������
            *p_cont_disp = *p_cont_gas;
            calc_F_and_DF( paramsc, *p_cont_disp, left_params, M, c1l, DISPERSED_PHASE, &F, &DF ); 
            *v_cont_disp = left_params[V_DISP] - F;
            break;
        case RIGHT_ONLY_GRAD:
            // ��������� ���� ������������ ������ �� ������� � ����������� ����� �� �������, ������ ������� ������������, ������ ������
            // Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
            // two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 501 + �����������
            *p_cont_disp = *p_cont_gas;
            calc_F_and_DF( paramsc, *p_cont_disp, right_params, M, c1r, DISPERSED_PHASE, &F, &DF );
            *v_cont_disp = right_params[V_DISP] + F;
            break;
        default:
            printf( "\ncalc_p_v_initial_guess -> wrong solver_part value\n\n" );
            exit( EXIT_FAILURE );
    }

    return SUCCESS;

}

/* ������ �������� � ��������� ����� ��� �� ���������� ������� � ���������� ����
   
   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   debug_info - ��������� � ���������� ����������� (in)
   left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
   right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
   sound_velocities[K_GENERAL_CASE] - ������ �� ���������� ����� ��� �� ������ ������� �� ������� - 
                                      sound_velocities[GAS_LEFT] - � ���� �����, sound_velocities[GAS_RIGHT] - � ���� ������,
                                      sound_velocities[DISP_LEFT] - � ���������� ���� �����, sound_velocities[DISP_RIGHT] - � ���������� ���� ������ (in)
   solver_part - ���������, ������� ����������, ����� �� ������������ �������� ���� ���������� ���� ����������� (in)
   
   solid_discontinuity_pressures[K_GENERAL_CASE] - �������� ���� � ���������� ���� �� ������ ������� �� ����������� �������
                                                   � ���������� ����, �� ����� - ��������� �����������, �� ������ - ���������;
                                                   [GAS_LEFT] - �������� ����� � ����, [GAS_RIGHT] - �������� ������ � ����,
                                                   [DISP_LEFT] - �������� ����� � ���������� ����, [DISP_RIGHT] - �������� ������ � ���������� ���� (in/out)
   v1 - �������� ����������� ������� � ���������� ���� (out)
   v2 - �������� ����������� ������� � ������� ���� (out)
   v21 - �������� ���� ����� �� ����������� ������� ���������� ���� (out)
   v22 - �������� ���� ������ �� ����������� ������� ���������� ���� (out)
   
   ����������: SUCCEESS         �������� ���������� ��������
               GODUNOV_FAILS    �������� �� �������, ������� ������ ������ ��������� �� ������� */ 
ReturnCodes calc_p_v( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M], double right_params[M],
                      double sound_velocities[K_GENERAL_CASE], Disp_phase_cases solver_part,
                      double solid_discontinuity_pressures[K_GENERAL_CASE], double *v1, double *v2, double *v21, double *v22 ) {
    
    int return_code;    /* ��� �������� ������� */

    /* ��� ���������� ���� - � ����������� �� ������������ �������� ���� ��������� ���� */
    switch ( solver_part ) {
        case BOTH_GRAD:
            return_code = calc_solid_discontinuity_pressures( paramsc, debug_info, left_params, right_params, BOTH_GRAD,
                sound_velocities, solid_discontinuity_pressures, v1, v2, v21, v22 );
            debug_info->godunov_cases[1]++;
            break;
        case LEFT_ONLY_GRAD:
            /* ��������� ���� ������������ ����� �� ������� � ����������� ������ �� �������, ������ ������
               Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of
               compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 502 + ����������� */
            return_code = calc_solid_discontinuity_pressures( paramsc, debug_info, left_params, right_params, LEFT_ONLY_GRAD,
                sound_velocities, solid_discontinuity_pressures, v1, v2, v21, v22 );
            /* �������� � ���������� ���� ������ �� ������� �� ����� ������, ��� �������� ������ �������� �������� */
            solid_discontinuity_pressures[DISP_RIGHT] = paramsc->background_pressure;
            debug_info->godunov_cases[2]++;
            break;
        case RIGHT_ONLY_GRAD:
            /* ��������� ���� ������������ ������ �� ������� � ����������� ����� �� �������, ������ ������
               ���������� Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a
               model of compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 502 + ����������� */
            return_code = calc_solid_discontinuity_pressures( paramsc, debug_info, left_params, right_params, RIGHT_ONLY_GRAD,
                sound_velocities, solid_discontinuity_pressures, v1, v2, v21, v22 );
             /* �������� � ���������� ���� ����� �� ������� �� ����� ������, ��� �������� ������ �������� �������� */
             solid_discontinuity_pressures[DISP_LEFT] = paramsc->background_pressure;
             debug_info->godunov_cases[2]++;
            break;
        default:
            printf( "\ncalc_p_v -> wrong solver_part value\n\n" );
            exit( EXIT_FAILURE );

    }

    // � ������ ������������� ����������� ������� �������� ����������� ������� �������� ������ ������ �����,
    // ����� �������� ���������������� � ��������� ������ �������
    if ( fabs( *v1 ) < paramsc->eps_contact ) *v1 = paramsc->eps_contact;
    if ( fabs( *v2 ) < paramsc->eps_contact ) *v2 = paramsc->eps_contact;
    if ( fabs( *v21 ) < paramsc->eps_contact ) *v21 = paramsc->eps_contact;
    if ( fabs( *v22 ) < paramsc->eps_contact ) *v22 = paramsc->eps_contact;

    return SUCCESS;

}

/* ������� ������� ���������� ���������������� ������� �� ���������� ����� ��� ����� � ������ �� �������

   c2l - �������� ����� � ������� ���� ����� �� ������� (in)
   c2r - �������� ����� � ������� ���� ������ �� ������� (in)
   c1l - �������� ����� � ���������� ���� ����� �� ������� (in)
   c1r - �������� ����� � ���������� ���� ������ �� ������� (in)

   sound_velocities[K_GENERAL_CASE] - ������ ��������� ����� (out) */
void fill_sound_velocities( double c2l, double c2r, double c1l, double c1r, double sound_velocities[K_GENERAL_CASE] ) {

    sound_velocities[GAS_LEFT] = c2l;
    sound_velocities[GAS_RIGHT] = c2r;
    sound_velocities[DISP_LEFT] = c1l;
    sound_velocities[DISP_RIGHT] = c1r;

}

/* ������� ������� ������������� ������� �������� ����� � ������ �� ������� � ���������� ����

   p_cont_gas - �������� ���� �� ���������� ������� ��� ����� ������� �������� ���� ���������� ���� (in)
   p_cont_solid - �������� ���������� ���� �� ���������� ������� ��� ����� ������� �������� ���� ���������� ���� (in)

   solid_discontinuity_pressures[K_GENERAL_CASE] - ������ ��������� �������� ����� � ������ �� ������� �������� ����
   ���������� ���� (out) */
void init_solid_discontinuity_pressures( double p_cont_gas, double p_cont_solid,
                                         double solid_discontinuity_pressures[K_GENERAL_CASE] ) {

    solid_discontinuity_pressures[GAS_LEFT] = p_cont_gas;
    solid_discontinuity_pressures[GAS_RIGHT] = p_cont_gas;
    solid_discontinuity_pressures[DISP_LEFT] = p_cont_solid;
    solid_discontinuity_pressures[DISP_RIGHT] = p_cont_solid;

}

/* ������������ ��������� ������� �������� � �������� �� ���������� ������� � ����� ��� ������� �������� ���� ���������� ����

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 155. - Subroutine STARPU.

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   debug_info - ��������� � ���������� ����������� (in)
   v_ncons_l - ������ ����������� ���������� ����� �� ������� (in)
   v_ncons_r - ������ ����������� ���������� ������ �� ������� (in)
   vector_size - ������ ������� ���������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   phase - ������������� ����, ��� ������� ������ ��������� ����������� - ������� ��� ���������� (in)

   p_cont - �������� �� ���������� ������� (out)
   v_cont - �������� �� ���������� ������� (out)
   
   ����������: SUCCEESS         �������� ���������� ��������
               GODUNOV_FAILS    �������� �� �������, ������� ������ ������ ��������� �� ������� */
int calc_contact_pressure_velocity( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double *v_ncons_l, double *v_ncons_r,
                                    int vector_size, double cl, double cr, Phase phase, double *p_cont, double *v_cont ) {

    double vl, vr;      /* �������� ����� � ������ �� ������� */
    double p_old;       /* �������� �������� �� ���������� �������� */
    double fl, fr;      /* �������� ������� */
    double fld, frd;    /* �������� ����������� */
    int iter_num = 0;   /* ���������� ����������� �������� */
    double criteria;    /* ���������� ��� ����������� ���������� */
    double g;           /* ���������� �������� */
    double p_prev = 0;// �������� ��  ���������� �������� ����� ����;
    if ( vector_size == M_REDUCTION ) {
        /* ����������� ���������� ������ */
        vl = v_ncons_l[V];
        vr = v_ncons_r[V];
        printf("\n cl = %lf cr = %lf", cl, cr);
    }
    switch ( phase ) {
        case GAS_PHASE:
            if ( vector_size == M ) {
                /* ������ ���������� ������ */
                vl = v_ncons_l[V_GAS];
                vr = v_ncons_r[V_GAS];
            }
            g = paramsc->g2;
            break;
        case DISPERSED_PHASE:
            if ( vector_size == M ) {
                vl = v_ncons_l[V_DISP];
                vr = v_ncons_r[V_DISP];
            }
            g =paramsc->g1;
            break;
        default:
            printf( "\ncalc_contact_pressure_velocity -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    if ( 2.0 * ( cl + cr ) / ( g - 1.0 ) <= vr - vl ) {
        /* ������ ������������� ������� */
        printf( "\ncalc_contact_pressure_velocity -> vacuum is generated in Godunov flux calculation " );
        debug_print( debug_info, phase );
        exit( EXIT_FAILURE );
    }

    /* ������ ���������� ����������� ��� �������� */
    p_old = pressure_initial_guess( paramsc, v_ncons_l, v_ncons_r, vector_size, cl, cr, phase );
    if ( p_old < 0.0 ) {
        printf( "\ncalc_contact_pressure_velocity -> initial pressure guess is negative " );
        debug_print( debug_info, phase );
        printf(" \n pressure_old = %lf", p_old);
        exit( EXIT_FAILURE );
    }
    
    /* ������� ����������� ��������� ��� ���������� �������� �� ���������� ������� ������� �������-������� */
    do {
        calc_F_and_DF( paramsc, p_old, v_ncons_l, vector_size, cl, phase, &fl, &fld );
        calc_F_and_DF( paramsc, p_old, v_ncons_r, vector_size, cr, phase, &fr, &frd );
        *p_cont = p_old - ( fl + fr + vr - vl ) / ( fld + frd );
        printf("\n p_cont = %lf, iter_num = %d fl = %lf, fr = %lf, vl = %lf, vr = %lf, fld = %lf, frd = %lf", *p_cont, iter_num, fl, fr, vl, vr, fld, frd);
        criteria = 2.0 * fabs( ( *p_cont - p_old ) / ( *p_cont + p_old ) );
        if(iter_num == 0)
            p_prev = *p_cont;
        printf( "\n %lf %lf %lf", p_prev, p_old, *p_cont);
        if( iter_num>1){
            if(*p_cont - p_prev < paramsc->eps_general){
                *p_cont = max(*p_cont, p_old);
                criteria = paramsc->eps_general;
                printf("check");
            }
            p_prev = p_old;
        }
        
        iter_num++;
        if ( iter_num > paramsc->max_iter_num ) {
            printf( "\ncalc_contact_pressure_velocity -> number of iterations exceeds the maximum value " );
            debug_print( debug_info, phase );
            return GODUNOV_FAILS;
        }
        if ( *p_cont < 0.0 && p_old < 0) {
            printf( "\ncalc_contact_pressure_velocity -> pressure is negative \n" );
            debug_print( debug_info, phase );
            printf(" \n pressure_cont = %lf", *p_cont);
            /* ���������� ������� ��� �������� */
            /* draw_adiabatic_curve( params, v_ncons_l, vector_size, cl, phase, LEFT );
            draw_adiabatic_curve( params, v_ncons_r, vector_size, cr, phase, RIGHT );
             */exit( EXIT_FAILURE );
            return GODUNOV_FAILS;
        }
        p_old = *p_cont;
    } while ( criteria > paramsc->eps_general );

    /* �������� ����������� ������� */
    *v_cont = 0.5 * ( vl + vr + fr - fl );
    return SUCCESS;

}

/* ����������� ���������� ����������� ��� ������� �������� �� ���������� ������� � ����� ��� ������� �������� ���� ���������� ����

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons_l - ������ ����������� ���������� ����� �� ������� (in)
   v_ncons_r - ������ ����������� ���������� ������ �� ������� (in)
   vector_size - ������ ������� ���������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   phase - ������������� ����, ��� ������� ������ ��������� ����������� - ������� ��� ���������� (in)

   ���������� ������� ��������� ����������� */
double pressure_initial_guess( struct ParametersCommon *paramsc, double *v_ncons_l, double *v_ncons_r, int vector_size,
                               double cl, double cr, int phase ) {

    double rl, vl, pl;                  /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;                  /* ����������� ���������� ������ �� ������� */
    double g;                           /* ���������� �������� */
    double p01;                          /* �������� ��� ����������� ��������� ��������� */
    /* ��������� �����������, ������������ �� ����������� ������������ ��������������� �������
       � ����������� ���������� */
    double p_lin;
    double p_min, p_max;                /* ����������� � ������������ �������� ����� � ������ �� ������� */
    double p_ratio;                     /* ������� �� �������� ����� � ������ �� ������� */
    double p1, p2, g1, g2;              /* ��������������� ���������� ��� ������������� �������� */
    double p_cand;                      /* �������� �� ���� ���������� ����������� */
    
    if ( vector_size == M_REDUCTION ) {
        /* ����������� ���������� ������ */
        rl = v_ncons_l[R];
        vl = v_ncons_l[V];
        pl = v_ncons_l[P];
        rr = v_ncons_r[R];
        vr = v_ncons_r[V];
        pr = v_ncons_r[P];
    }
    switch ( phase ) {
        case GAS_PHASE:
            if ( vector_size == M ) {
                rl = v_ncons_l[R_GAS];
                vl = v_ncons_l[V_GAS];
                pl = v_ncons_l[P_GAS];
                rr = v_ncons_r[R_GAS];
                vr = v_ncons_r[V_GAS];
                pr = v_ncons_r[P_GAS];
            }
            g = paramsc->g2;
            p01 = 0.0;
            break;
        case DISPERSED_PHASE:
            if ( vector_size == M ) {
                rl = v_ncons_l[R_DISP];
                vl = v_ncons_l[V_DISP];
                pl = v_ncons_l[P_DISP];
                rr = v_ncons_r[R_DISP];
                vr = v_ncons_r[V_DISP];
                pr = v_ncons_r[P_DISP];
            }
            g = paramsc->g1;
            p01 =paramsc->p01;
            break;
        default:
            printf( "\npressure_initial_guess -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    /* ��������� ����������� �� �������� ������
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = min( pl, pr );
    p_max = max( pl, pr );
    p_ratio = p_max / p_min;
    printf("\n p_ratio = %lf", p_ratio);
    if ( ( p_ratio <= paramsc->p_max_ratio ) &&
        ( ( p_min < p_lin && p_lin < p_max ) || ( fabs( p_min - p_lin ) < paramsc->eps_general || fabs( p_max - p_lin ) < paramsc->eps_general ) ) ) {
        /* ��������� ����������� �� ��������������� ������ */
        p_cand = p_lin;
        printf("\n pressure initial guess linear p_lin = %lf", p_lin);
    } else {
        if ( p_lin < p_min ) {
            /* ��������� ����������� �� ���� ������ ����������
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 301. - Formula (9.32) + �������� �� ���������� ��������� ��������� */
            g1 = 0.5 * ( g - 1.0 ) / g;
            p_cand = pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / pow( pl + p01, g1 ) + cr / pow( pr + p01, g1 ) ) ),
                1 / g1 ) - p01;
            printf("\n pressure initial guess two fan p_cand = %lf  pr = %lf  pl = %lf cr = %lf cl = %lf vr = %lf, vl = %lf", p_cand, pr, pl, cr, cl, vr, vl);
        } else {
            /* ��������� ����������� �� ���� ������� ������
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48) + �������� �� ���������� ��������� ��������� */
            g1 = 2.0 / ( g + 1.0 );
            g2 = ( g - 1.0 ) / ( g + 1.0 );
            p1 = sqrt( g1 / rl / ( g2 * ( pl + p01 ) + p_lin + p01 ) );
            p2 = sqrt( g1 / rr / ( g2 * ( pr + p01 ) + p_lin + p01 ) );
            p_cand = ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
            printf("\n pressure initial guess two shock p_cand = %lf  pr = %lf  pl = %lf", p_cand, pr, pl);
        }
    }
    
    if ( p_cand < 0.0 ) {
        /* ���� �� ��������� �������� Toro, �� ������� �������� �� "��������� ������� �������"
           ������� �.�. � ��. ��������� ������� ����������� ����� ������� ��������. - 
           �.: �����, 1976. - �. 113. - ������� (13.26). */
        p_cand = ( pl * rr * cr + pr * rl * cl + ( vl - vr ) * rl * cl * rr * cr ) / ( rl * cl + rr * cr );
        printf("\n pressure initial guess sonic discontinue p_cand = %lf  pr = %lf  pl = %lf ", p_cand, pr, pl);
    }

    if ( p_cand < 0.0 ) {
        /* ���� ������ �� �������, �������� ��������� ���� �����-�� ��������� ����������� */
        p_cand = 0.5 * ( pl + pr );
        printf("\n pressure initial guess average p_cand = %lf  pr = %lf  pl = %lf ", p_cand, pr, pl);
    }
    printf("\n final guess p_cand = %lf", p_cand);
    return p_cand;

}

/* ������ ������� F, ������������ �������� ����� �� ���������� ������� � ����� ��� ������� �������� ���� ���������� ����, � �� ����������� �� �������� ����� DF

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + �������� �� ���������� ��������� ���������:
   ������� �.�. � ��. ��������� ������� ����������� ����� ������� ��������. - 
   �.: �����, 1976. - �. 110 - 111. - ������� (13.16), (13.17).

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   curr_press - �������� � ���������� �������� (in)
   v_ncons - ������ ����������� ���������� (in)
   vector_size - ������ ������� ���������� (in)
   c - �������� ����� � ������������� ����� (in)
   phase - ������������� ����, ��� ������� ������ ��������� ����������� - ������� ��� ���������� (in)

   F - �������� ������� (out)
   DF - �������� ����������� (out) */
void calc_F_and_DF( struct ParametersCommon *paramsc, double curr_press, double *v_ncons, int vector_size, double c, int phase, double *F, double *DF ) {

    double r, v, p; // ����������� ����������
    double g; // ���������� ��������
    double p01; // �������� ��� ����������� ��������� ���������
    double p_ratio, fg, q; // ��������������� ����������

    if ( vector_size == M_REDUCTION ) {
        // ����������� ���������� ������
        r = v_ncons[R];
        v = v_ncons[V];
        p = v_ncons[P];
     }
    switch ( phase ) {
        case GAS_PHASE:
            if ( vector_size == M ) {
                r = v_ncons[R_GAS];
                v = v_ncons[V_GAS];
                p = v_ncons[P_GAS];
            }
            g = paramsc->g2;
            p01 = 0.0;
            break;
        case DISPERSED_PHASE:
            if ( vector_size == M ) {
                r = v_ncons[R_DISP];
                v = v_ncons[V_DISP];
                p = v_ncons[P_DISP];
            }
            g = paramsc->g1;
            p01 =paramsc->p01;
            break;
        default:
            printf( "\ncalc_F_and_DF -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    p_ratio = ( curr_press + p01 ) / ( p + p01 );
    if ( curr_press <= p ) {
        // ����� ����������
        fg = 2.0 / ( g - 1.0 );
        *F = fg * c * ( pow( p_ratio, 1.0 / fg / g ) - 1.0 );
        *DF = ( 1.0 / r / c ) * pow( p_ratio, - 0.5 * ( g + 1.0 ) / g );
    }
    else {
        // ������� �����
        double a, b;
        a = 2.0 / r / ( g + 1.0 );
        b = ( g - 1.0 ) * p / ( g + 1.0 );
        q = sqrt( a / ( b + curr_press ) );
        *F = ( curr_press - p ) * q;
        *DF = ( 1.0 - 0.5 * ( curr_press - p ) / ( b + curr_press ) ) * q;
        /*q = sqrt( 0.5 * ( g + 1.0 ) / g * p_ratio + 0.5 * ( g - 1.0 ) / g );
        *F = ( curr_press - p ) / c / r / q;
        *DF = 0.25 * ( ( g + 1.0 ) * p_ratio + 3 * g - 1.0 ) / g / r / c / pow( q, 3.0 );*/
    }

}

/* ������ ������� G, ������������ ��������� ���� �� ������ ������� �� ����������� ������� � ����� ��� ������� �������� ���� ���������� ����,
   � �� ����������� �� �������� DG

   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - ������� G ������������ ��������� (7).

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   curr_press - �������� � ���������� �������� (in)
   v_ncons[M] - ������ ����������� ���������� (in)
   
   G - �������� ������� (out)
   DG - �������� ����������� (out) */
void calc_G_and_DG( struct ParametersCommon *paramsc, double curr_press, double v_ncons[M], double *G, double *DG ) {

    double r = v_ncons[R_GAS], p = v_ncons[P_GAS]; // ��������� � �������� ���� � ������������� �����
    double p_ratio; // ��������� �������� �������� �� ���������� ������� � �������� � ������������� �����
    double g1 = paramsc->g2 - 1.0, g2 = paramsc->g2 + 1.0; // ���������� ��� �������� ����������
    double A, B;
        
    p_ratio = curr_press / p;
    if ( curr_press < p ) {
        // ����� ����������
        *G = r * pow( p_ratio, 1.0 / paramsc->g2 );
        *DG = *G / paramsc->g2 / curr_press;
    }
    else {
        // ������� �����
        A = g1 + g2 * p_ratio;
        B = g1 * p_ratio + g2;
        *G = r * A / B;
        *DG = r * ( g2 * B - g1 * A ) / p / pow( B, 2.0 );
    }

}

/* ������������ ��������� ��� ������ �������� � ������� � ���������� ����� �� ���������� ������� � ���������� ����
   ������� �������-�������
   
   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   debug_info - ��������� � ���������� ����������� (in)
   left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
   right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
   dsolver_part - ���������, ������� ����������, ����� �� ������������ �������� ���� ���������� ���� ����������� (in)
   sound_velocities[K_GENERAL_CASE] - ������ �� ���������� ����� ��� �� ������ ������� �� ������� - 
                                      sound_velocities[GAS_LEFT] - � ���� �����, sound_velocities[GAS_RIGHT] - � ���� ������,
                                      sound_velocities[DISP_LEFT] - � ���������� ���� �����,
                                      sound_velocities[DISP_RIGHT] - � ���������� ���� ������ (in)

   solid_discontinuity_pressures[K_GENERAL_CASE] - �������� ���� � ���������� ���� �� ������ ������� �� ����������� �������
                                                   � ���������� ����, �� ����� - ��������� �����������, �� ������ - ���������;
                                                   [GAS_LEFT] - �������� ����� � ����, [GAS_RIGHT] - �������� ������ � ����,
                                                   [DISP_LEFT] - �������� ����� � ���������� ����, [DISP_RIGHT] - �������� ������ � ���������� ���� (in/out)
   v1 - �������� ����������� ������� � ���������� ���� (out)
   v2 - �������� ����������� ������� � ������� ���� (out)
   v21 - �������� ���� ����� �� ����������� ������� � ���������� ���� (out)
   v22 - �������� ���� ������ �� ����������� ������� � ���������� ���� (out)
   
   ����������: SUCCEESS         �������� ���������� ��������
               GODUNOV_FAILS    �������� �� �������, ������� ������ ������ ��������� �� ������� */
int calc_solid_discontinuity_pressures( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params[M],
                                        double right_params[M], Disp_phase_cases solver_part, double sound_velocities[K_GENERAL_CASE],
                                        double solid_discontinuity_pressures[K_GENERAL_CASE], double *v1, double *v2,
                                        double *v21, double *v22 ) {

    int actual_sol_cont_pres_size; /* ���������� �������������� �������� - K_GENERAL_CASE � ����� ������ � K_SPECIAL_CASE � ������ */
    // all arrays are of the maximum possible constant size
    double dp[M]; /* ���������� �������-������� � ������������ �������� */
    double P[M]; /* ������-�������, ������������ ����������� �� ������� �������� ���� � ���������� ���� */
    double DP[M][M]; 

    double criteria; /* ���������� ��� ����������� ���������� */
    double average_pressures[M]; /* ������ �� �������� ���������� ����� ���������� - ��� ������� �������� �������� */
    double check_system_sol[M]; /* ������ ��� �������� ������� ���� */

    // solution of system of non-linear algebraic equations to find both phases pressures on the solid discontinuity
    // by Newton-Raphson method

    /* ����������� ������� �������� ������� */
    switch ( solver_part ) {
        case BOTH_GRAD:
            actual_sol_cont_pres_size = K_GENERAL_CASE; /* ����� ������ ������� ���������� ���� �� ��� ������� �� ������� */
            break;
        case LEFT_ONLY_GRAD:
        case RIGHT_ONLY_GRAD:
            actual_sol_cont_pres_size = K_SPECIAL_CASE; /* ������ ������ ���������� ���������� ���� �� ���� �� ������ �� ������� */
            break;
        default:
            printf( "\ncalc_solid_discontinuity_pressures -> wrong value of solver_part variable.\n" );
            exit( EXIT_FAILURE );
    }

    do {
        
        /* ������ ������������ ������� � ������� ����� */
        switch ( solver_part ) {
            case BOTH_GRAD:
                /* ����� ������ ������� ���������� ���� �� ��� ������� �� ������� */
                calc_P_and_DP_general_case( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case LEFT_ONLY_GRAD:
                /* ������ ������ - ���������� ���� ������ ����� �� ������� */
                calc_P_and_DP_left_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case RIGHT_ONLY_GRAD:
                /* ������ ������ - ���������� ���� ������ ������ �� ������� */
                calc_P_and_DP_right_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            default:
                printf( "\ncalc_solid_discontinuity_pressures -> wrong value of solver_part variable.\n" );
                exit( EXIT_FAILURE );
        }
                
        // right-hand side vector calculation for Newton iterations
        for ( int i = 0; i < actual_sol_cont_pres_size; i++ )
            P[i] = - P[i];
        
        // solution of the system for residuals
        solve_linear_system( DP, P, actual_sol_cont_pres_size, paramsc->eps_general, dp );
        if ( paramsc->is_debug ) {
            // check the correctness of system solution
            mult_matrix_vector( actual_sol_cont_pres_size, DP, dp, check_system_sol );
            if ( !compare_vectors( actual_sol_cont_pres_size, check_system_sol, P, paramsc->eps_general ) ) {
                printf( "\ncalc_solid_discontinuity_pressures -> linear system of equations is solved incorrectly.\n" );
                exit( EXIT_FAILURE );
            }
        }
        
        // solution vector update
        if ( solver_part == BOTH_GRAD || solver_part == LEFT_ONLY_GRAD ) {
            for ( int i = 0; i < actual_sol_cont_pres_size; i++ )
                /* ����������� 4 ���������� � ����� ������ ��� 3 ������ � ������ ������� ��������� ���� ������ ����� �� ������� */
                solid_discontinuity_pressures[i] += dp[i];
        }
        else {
            solid_discontinuity_pressures[GAS_LEFT] += dp[0];   /* �������� ���� ����� �� ����������� ������� � ���������� ���� */
            solid_discontinuity_pressures[GAS_RIGHT] += dp[1];  /* �������� ���� ������ �� ����������� ������� � ���������� ���� */
            solid_discontinuity_pressures[DISP_RIGHT] += dp[2]; /* �������� ���������� ���� ������ �� ����������� ������� � ���������� ���� */
        }

        // stop criteria calculation
        if ( solver_part == BOTH_GRAD || solver_part == LEFT_ONLY_GRAD ) {
            for ( int i = 0; i < actual_sol_cont_pres_size; i++ )
                average_pressures[i] = solid_discontinuity_pressures[i] - 0.5 * dp[i];
        }
        else {
            average_pressures[0] = solid_discontinuity_pressures[GAS_LEFT] - 0.5 * dp[0];   /* ������������� �������� ���� ����� ��
                                                                                               ����������� ������� � ���������� ���� */
            average_pressures[1] = solid_discontinuity_pressures[GAS_RIGHT] - 0.5 * dp[1];  /* ������������� �������� ���� ������ ��
                                                                                               ����������� ������� � ���������� ���� */
            average_pressures[2] = solid_discontinuity_pressures[DISP_RIGHT] - 0.5 * dp[2]; /* ������������� �������� ���������� ���� ������ ��
                                                                                               ����������� ������� � ���������� ���� */
        }
        criteria = norm( dp, actual_sol_cont_pres_size ) / norm( average_pressures, actual_sol_cont_pres_size );

        // the necessary check of the correctness of the resulats
        for ( int i = 0; i < actual_sol_cont_pres_size; i++ ) {
            if ( solid_discontinuity_pressures[i] < 0.0 ) {
                
                /* ����������� ���������� ������������ �������� ���������� ���������� � ����� ������ ������ */
                /* printf( "\ncalc_solid_discontinuity_pressures -> pressure is negative " );
                debug_print( debug_info, i ); */
                
                return GODUNOV_FAILS;
            }
        }

    } while ( criteria > paramsc->eps_general );

    if ( paramsc->is_debug ) {
        /* ��������, ��� ��������� ������� ������������� �������� ���������� ������� */
        switch ( solver_part ) {
            case BOTH_GRAD:
                /* ����� ������ ������� ���������� ���� �� ��� ������� �� ������� */
                calc_P_and_DP_general_case( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case LEFT_ONLY_GRAD:
                /* ������ ������ - ���������� ���� ������ ����� �� ������� */
                calc_P_and_DP_left_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            case RIGHT_ONLY_GRAD:
                /* ������ ������ - ���������� ���� ������ ������ �� ������� */
                calc_P_and_DP_right_disp_phase( paramsc, left_params, right_params, sound_velocities, solid_discontinuity_pressures,
                    P, DP, v1, v2, v21, v22 );
                break;
            default:
                printf( "\ncalc_solid_discontinuity_pressures -> wrong value of disp_phase variable.\n" );
                exit( EXIT_FAILURE );
        }
        for ( int i = 0; i < actual_sol_cont_pres_size; i++ ) {
            if ( fabs( P[i] ) > paramsc->eps_thin_layer ) {
                printf( "\ncalc_solid_discontinuity_pressures -> the 'thin-layer' algebraic system of non-linear equations is solved wrong " );
                debug_print( debug_info, i );
                printf( "calc_solid_discontinuity_pressures -> P[%d] = %e instead of acceptable tolerance %e.\n", i, P[i], paramsc->eps_thin_layer );
                exit( EXIT_FAILURE );
            }
        }
    }
    return SUCCESS;

}

/* ������ ������-������� ��� ����������� ����������� �� ���������� ������� � ���������� ����, �� ������� �����, � ����� ���������
   ���������� �������� � ��������� ���� ����� � ������ �� ����������� ������� � ���������� ����.
   ����� ������ ������� ���������� ���� �� ��� ������� �� �������.
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - ������� ������� ������������ ���������
   (23), (25).

   ����� ������� ����������� ���������� � ����� ���������:
   - ������ ������ - ����� ����: 1 - ���������� ����, 2 - ������� ����
   - ������ ������ - ��������� ������������ ����������� �������: 1 - �����, 2 - ������
   ������� l � r ������������� ���������� � ������������� ����� ����� � ������ �� �������, ��������������.
   
   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
   right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
   c[K_GENERAL_CASE] - ������ ��������� ����� � ���������� � ������� ����� �� ������ ������� �� ������� (in)
   curr_p[K_GENERAL_CASE] - ������� ������ �������� ���������� � ������� ��� ����� � ������ �� ������� ���������� (in)

   P[K_GENERAL_CASE] - ������� ������-������� (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - ������� ������� ����� (out)
   v_cont_disp - �������� ����������� ������� � ���������� ���� (out)
   v_cont_gas - �������� ����������� ������� � ������� ���� (out)
   v_gas_left - �������� ���� ����� �� ����������� ������� � ���������� ���� (out)
   v_gas_right - �������� ���� ������ �� ����������� ������� � ���������� ���� (out) */
void calc_P_and_DP_general_case( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                 double c[M], double curr_p[M], double P[M],
                                 double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                 double *v_gas_left, double *v_gas_right ) {

    double p11 = curr_p[DISP_LEFT], p12 = curr_p[DISP_RIGHT]; /* ������� �������� � ���������� ���� ����� � ������ �� ����������� ������� � ���������� ���� */
    double p21 = curr_p[GAS_LEFT], p22 = curr_p[GAS_RIGHT];   /* ������� �������� � ������� ���� ����� � ������ �� ����������� ������� � ���������� ���� */
    
    double c1l = c[DISP_LEFT], c1r = c[DISP_RIGHT]; /* �������� ����� � ���������� ���� ����� � ������ �� ��������������� ������� */
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   /* �������� ����� � ������� ���� ����� � ������ �� ��������������� ������� */

    double v1l = left_params[V_DISP], v1r = right_params[V_DISP];  /* �������� ���������� ���� ����� � ������ �� ��������������� ������� */
    double v2l = left_params[V_GAS], v2r = right_params[V_GAS];  /* �������� ������� ���� ����� � ������ �� ��������������� ������� */

    double b1l = left_params[B_DISP], b1r = right_params[B_DISP];  /* �������� ���� ���������� ���� ����� � ������ �� ��������������� ������� */
    double b2l = 1.0 - b1l, b2r = 1.0 - b1r;               /* �������� ���� ������� ���� ����� � ������ �� ��������������� ������� */

    double f11, df11;   /* ������� � ����������� ��� ����������� �������� ���������� ���� ����� �� ����������� ������� � ���������� ���� */
    double f12, df12;   /* ������� � ����������� ��� ����������� �������� ���������� ���� ������ �� ����������� ������� � ���������� ���� */

    double f21, df21;   /* ������� � ����������� ��� ����������� �������� ������� ���� ����� �� ����������� ������� � ���������� ���� */
    double f22, df22;   /* ������� � ����������� ��� ����������� �������� ������� ���� ������ �� ����������� ������� � ���������� ���� */

    double g21, dg21;   /* ������� � ����������� ��� ����������� ��������� ������� ���� ����� �� ����������� ������� � ���� */
    double g22, dg22;   /* ������� � ����������� ��� ����������� ��������� ������� ���� ����� �� ����������� ������� � ���� */

    double v11, v12;    /* ������� �������� ���������� ���� ����� � ������ �� ����������� ������� � ���������� ���� */
    double v21, v22;    /* ������� �������� ������� ���� ����� � ������ �� ����������� ������� � ���������� ���� */

    double r1, r2;  /* ��������� ���� ����� � ������ �� ����������� ������� � ���� */

    double g, A, B, C, D, E, dv2, dv1, dv2c;    /* ��������������� ���������� ��� ���������� ������ ���������� */

    /* �� �������� ������� ��������� ���������� �������� � ���������� � ������� ����� ����������, ����� �� ���� ���������
       ������������ ����������� */
    /* ���������� ���� */
    calc_F_and_DF( paramsc, p11, left_params, M, c1l, DISPERSED_PHASE, &f11, &df11 );
    v11 = v1l - f11;
    calc_F_and_DF( paramsc, p12, right_params, M, c1r, DISPERSED_PHASE, &f12, &df12 );
    v12 = v1r + f12;
    /* ������� ���� */
    calc_F_and_DF( paramsc, p21, left_params, M, c2l, GAS_PHASE, &f21, &df21 );
    v21 = v2l - f21;
    calc_F_and_DF( paramsc, p22, right_params, M, c2r, GAS_PHASE, &f22, &df22 );
    v22 = v2r + f22;

    /* ������ ���������� ���� ����� � ������ �� ����������� ������� � ���� �� �������� �������� �� ���������� ������� */
    calc_G_and_DG( paramsc, p21, left_params, &g21, &dg21 );
    r1 = g21;
    calc_G_and_DG( paramsc, p22, right_params, &g22, &dg22 );
    r2 = g22;

    /* ������ ��������������� ���������� */
    g = paramsc->g2 / ( paramsc->g2 - 1.0 );
    A = pow( p22 / p21, 1.0 / paramsc->g2 );
    B = b1r * p12 + b2r * p22 - b1l * p11 - b2l * p21;
    dv2 = v22 - v12;
    dv1 = v21 - v11;
    dv2c = v22 - v21;
    C = 0.5 * ( pow( dv2, 2.0 ) - pow( dv1, 2.0 ) );
    D = b2r * A * dv2;
    E = g / r1 / A;

    /* ������ ����������� ������-������� P */
    P[0] = v12 - v11;
    if ( v21 > v11 ) {
        P[1] = D - b2l * dv1;
        P[2] = B + b2l * r1 * dv1 * dv2c;
        P[3] = E * p22 - g * p21 / r1 + C;
        *v_cont_gas = v22;  /* ������� �������� ����������� ������� � ������� ���� */
    }
    else {
        P[1] = - b2l / A * dv1 + b2r * dv2;
        P[2] = B + b2r * r2 * dv2 * dv2c;
        P[3] = - g * p21 * A / r2 + g * p22 / r2 + C;
        *v_cont_gas = v21;  /* ������� �������� ����������� ������� � ������� ���� */
    }

    /* ������ ����������� ������� ����� */
    /* ������ ������ */
    DP[0][0] = 0.0;
    DP[0][1] = 0.0;
    DP[0][2] = df11;
    DP[0][3] = df12;
    if ( v21 > v11 ) {
        /* ������ ������ */
        DP[1][0] = - D / p21 / paramsc->g2 + b2l * df21;
        DP[1][1] = D / p22 / paramsc->g2 + b2r * A * df22;
        DP[1][2] = - b2l * df11;
        DP[1][3] = - b2r * A * df12;
        /* ������ ������ */
        DP[2][0] = b2l * ( - 1.0 + dv1 * dv2c * dg21 + r1 * df21 * ( dv1 - dv2c ) );
        DP[2][1] = b2r + b2l * r1 * dv1 * df22;
        DP[2][2] = - b1l + b2l * r1 * dv2c * df11;
        DP[2][3] = b1r;
        /* ��������� ������ */
        DP[3][0] = g / r1 * ( p22 / A / paramsc->g2 / p21 - 1.0 ) + dv1 * df21 + dg21 * g / pow( r1, 2.0 ) * ( p21 - p22 / A );
        DP[3][1] = 1.0 / r1 / A + dv2 * df22;
        DP[3][2] = - dv1 * df11;
        DP[3][3] = - dv2 * df12;
    }
    else {
        /* ������ ������ */
        DP[1][0] = b2l / A * ( df21 - dv1 / paramsc->g2 / p21 );
        DP[1][1] = b2r * df22 + b2l * dv1 / A / paramsc->g2 / p22;
        DP[1][2] = - b2l / A * df11; DP[1][3] = - b2r * df12;
        /* ������ ������ */
        DP[2][0] = - b2l + b2r * r2 * dv2 * df21;
        DP[2][1] = b2r * ( 1.0 + dg22 * dv2 * dv2c + r2 * df22 * ( dv2c + dv2 ) );
        DP[2][2] = - b1l;
        DP[2][3] = b1r - b2r * r2 * dv2c * df12;
        /* ��������� ������ */
        DP[3][0] = - A / r2 + dv1 * df21;
        DP[3][1] = g / r2 + dv2 * df22 - A * p21 / ( paramsc->g2 - 1.0 ) / r2 / p22 + dg22 * g / pow( r2, 2.0 ) * ( p21 * A - p22 );
        DP[3][2] = - dv1 * df11;
        DP[3][3] = - dv2 * df12;
    }

    *v_cont_disp = v11; /* ������� �������� ����������� ������� � ���������� ���� */
    *v_gas_left = v21;  /* ������� �������� ���� ����� �� ����������� ������� � ���������� ���� */
    *v_gas_right = v22; /* ������� �������� ���� ������ �� ����������� ������� � ���������� ���� */

}

/* ������� ������ ������� � ���������� ����

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
   + �������� �� ���������� ��������� ���������: Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann
   problem and a high-resolution Godunov method for a model of compressible two-phase flow // Journal of Computational Physics.
   - 2006. - V. 212. - P. 490 - 526. - ������� (9) - (11).

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons_l[M] - ������ ���������������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ���������������� ���������� ������ �� ������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   p1 - �������� ����� �� ����������� ������� � ���������� ���� (in)
   p2 - �������� ������ �� ����������� ������� � ���������� ���� (in)
   v_cont - �������� �� ���������� ������� (in)
   s - �������� x/t, ��� �������� ���������� ������� (in)

   v_ncons_res[M] - ������ ���������������� ���������� � ����������� ������������ ��� ���������� ���� (out) */
void sample_solid_solution( struct ParametersCommon *paramsc, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr,
                            double p1, double p2, double v_cont, double s, double v_ncons_res[M], double cont_ncons[M] ) {

    double betal, rl, vl, pl;           /* ����������� ���������� ����� �� ������� */
    double betar, rr, vr, pr;           /* ����������� ���������� ������ �� ������� */
    double g1, g2, g3, g4, g5, g6, g7;  /* ��������������� ����������, ����������� �� ���������� ��������,
                                           � ������������ � Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics.
                                           - 2nd Edition. - Springer, 1999. - P. 153. */

    /* �������� ����� ���� */
    double shl, stl;    /* �������� "������" � "������" ����� ����� ���������� */
    double sl;          /* �������� ����� ������� ����� */

    /* �������� ������ ���� */
    double shr, str;    /* �������� "������" � "������" ������ ����� ���������� */
    double sr;          /* �������� ������ ������� ����� */

    double cml, cmr;        /* �������� ����� ����� � ������ �� ����������� ������� */
    double c;               /* ��������� �������� ����� ������ ����� ���������� */
    double p_ratio;
    double beta, r, v, p;   /* ���������� �������� �������� ����, ���������, �������� � �������� */

    /* ��������������� ���������� */
    /* ��������� ����� �� ������� */
    betal = v_ncons_l[B_DISP];
    rl = v_ncons_l[R_DISP];
    vl = v_ncons_l[V_DISP];
    pl = v_ncons_l[P_DISP];
    /* ��������� ������ �� ������� */
    betar = v_ncons_r[B_DISP];
    rr = v_ncons_r[R_DISP];
    vr = v_ncons_r[V_DISP];
    pr = v_ncons_r[P_DISP];
    /* ����������� �� ���������� �������� */
    g1 = 0.5 * ( paramsc->g1 - 1.0 ) / paramsc->g1;
    g2 = 0.5 * ( paramsc->g1 + 1.0 ) / paramsc->g1;
    g3 = 2.0 * paramsc->g1 / ( paramsc->g1 - 1.0 );
    g4 = 2.0 / ( paramsc->g1 - 1.0 );
    g5 = 2.0 / ( paramsc->g1 + 1.0 );
    g6 = ( paramsc->g1 - 1.0 ) / ( paramsc->g1 + 1.0 );
    g7 = 0.5 * ( paramsc->g1 - 1.0 );
    cont_ncons[B_DISP] = betal;
    if ( s <= v_cont ) {
        /* ��������������� ����� - ����� �� ����������� ������� */
        beta = betal;
        if ( p1 <= pl ) {
            /* ����� ����� ���������� */
            cont_ncons[R_DISP] = rl * pow( p1 / pl, 1.0 / paramsc->g1 );
            shl = vl - cl;
            if ( s <= shl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( ( p1 + paramsc->p01 ) / ( pl + paramsc->p01 ), g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* ��������� ����� �� ����������� ������� */
                    r = rl * pow( ( p1 + paramsc->p01 ) / ( pl + paramsc->p01 ), 1.0 / paramsc->g1 );
                    v = v_cont;
                    p = p1;
                }
                else {
                    /* ��������� ������ ����� ����� ���������� */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = ( pl + paramsc->p01 ) * pow( c / cl, g3 ) - paramsc->p01;
                }
            }
        }
        else {
            /* ����� ������� ����� */
            p_ratio = ( p1 + paramsc->p01 ) / ( pl + paramsc->p01 );
            cont_ncons[R_DISP] = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* ��������� �� ����� ������� ������ */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p1;
            }
        }
    }
    else {
        /* ��������������� ����� - ������ �� ����������� ������� */
        beta = betar;
        if ( p2 > pr ) {
            /* ������ ������� ����� */
            p_ratio = ( p2 + paramsc->p01 ) / ( pr + paramsc->p01 );
            cont_ncons[R_DISP] = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* ��������� �� ������ ������� ������ */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p2;
            }
        }
        else {
            /* ������ ����� ���������� */
            cont_ncons[R_DISP] = rr * pow( p1 / pr, 1.0 / paramsc->g1 );
            shr = vr + cr;
            if ( s >= shr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( ( p2 + paramsc->p01 ) / ( pr + paramsc->p01 ), g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* ��������� ������ �� ����������� ������� */
                   r = rr * pow( ( p2 + paramsc->p01 ) / ( pr + paramsc->p01 ), 1.0 / paramsc->g1 );
                   v = v_cont;
                   p = p2;
               }
               else {
                    /* ��������� ������ ������ ����� ���������� */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = ( pr + paramsc->p01 ) * pow( c / cr, g3 ) - paramsc->p01;
               }
            }
        }
    }
    
    /* ������������ ��������� ������� � ����������� */
    v_ncons_res[B_DISP] = beta;
    v_ncons_res[R_DISP] = r;
    v_ncons_res[V_DISP] = v;
    v_ncons_res[P_DISP] = p;
    
}

/* ������� ������ ������� � ������� ����

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons_l[M] - ������ ���������������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ���������������� ���������� ������ �� ������� (in)
   cl - �������� ����� � ���� ����� �� ������� (in)
   cr - �������� ����� � ���� ������ �� ������� (in)
   p1 - �������� ���� ����� �� ����������� ������� � ���������� ���� (in)
   p2 - �������� ���� ������ �� ����������� ������� � ���������� ���� (in)
   v_cont_solid - �������� ����������� ������� � ���������� ���� (in)
   v_cont_gas - �������� ����������� ������� � ������� ���� (in)
   v1 - �������� ���� ����� �� ����������� ������� � ���������� ���� (in)
   v2 - �������� ���� ������ �� ����������� ������� � ���������� ���� (in)
   s - �������� x/t, ��� �������� ���������� ������� (in)

   v_ncons_res[M] - ���������� ������ ���������������� ���������� (out) */
void sample_gas_solution( struct ParametersCommon *paramsc, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr, double p1,
                          double p2, double v_cont_solid, double v_cont_gas, double v1, double v2, double s, double v_ncons_res[M], double cont_ncons[M] ) {

    double rl, vl, pl;                  /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;                  /* ����������� ���������� ������ �� ������� */
    double g1, g2, g3, g4, g5, g6, g7;  /* ��������������� ����������, ����������� �� ���������� ��������,
                                           � ������������ � Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics.
                                           - 2nd Edition. - Springer, 1999. - P. 153. */

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

    /* ��������������� ���������� */
    /* ��������� ����� �� ������� */
    rl = v_ncons_l[R_GAS];
    vl = v_ncons_l[V_GAS];
    pl = v_ncons_l[P_GAS];
    /* ��������� ������ �� ������� */
    rr = v_ncons_r[R_GAS];
    vr = v_ncons_r[V_GAS];
    pr = v_ncons_r[P_GAS];
    /* ����������� �� ���������� �������� */    
    g1 = 0.5 * ( paramsc->g2 - 1.0 ) / paramsc->g2;
    g2 = 0.5 * ( paramsc->g2 + 1.0 ) / paramsc->g2;
    g3 = 2.0 * paramsc->g2 / ( paramsc->g2 - 1.0 );
    g4 = 2.0 / ( paramsc->g2 - 1.0 );
    g5 = 2.0 / ( paramsc->g2 + 1.0 );
    g6 = ( paramsc->g2 - 1.0 ) / ( paramsc->g2 + 1.0 );
    g7 = 0.5 * ( paramsc->g2 - 1.0 );

    if ( s <= v_cont_gas ) {
        /* ��������������� ����� - ����� �� ����������� ������� � ������� ���� */
        if ( p1 <= pl ) {
            /* ����� ����� ���������� */
            cont_ncons[R_GAS] = rl * pow( p1 / pl, 1.0 / paramsc->g2 );
            shl = vl - cl;
            if ( s <= shl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( p1 / pl, g1 );
                stl = v1 - cml;
                if ( s > stl ) {
                    /* ��������� ����� �� ����������� ������� */
                    if ( s <= v_cont_solid ) {
                        r = rl * pow( p1 / pl, 1.0 / paramsc->g2 );
                        v = v1;
                        p = p1;                        
                    }
                    else {
                        r = rl * pow( p2 / pl, 1.0 / paramsc->g2 );
                        v = v2;
                        p = p2;
                    }
                }
                else {
                    /* ��������� ������ ����� ����� ���������� */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = pl * pow( c / cl, g3 );
                }
            }
        }
        else {
            /* ����� ������� ����� */
            if ( s <= v_cont_solid )
                p_ratio = p1 / pl;
            else
                p_ratio = p2 / pl;
            cont_ncons[R_GAS] = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* ��������� �� ����� ������� ������ */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                if ( s <= v_cont_solid ) {
                    v = v1;
                    p = p1;
                }
                else {
                    v = v2;
                    p = p2;
                }
            }
        }
    }
    else {
        /* ��������������� ����� - ������ �� ����������� ������� */
        if ( p2 > pr ) {
            /* ������ ������� ����� */
            if ( s <= v_cont_solid )
                p_ratio = p1 / pr;
            else
                p_ratio = p2 / pr;
            cont_ncons[R_GAS] = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* ��������� �� ������ ������� ������ */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                if ( s <= v_cont_solid ) {
                    v = v1;
                    p = p1;
                }
                else {
                    v = v2;
                    p = p2;
                }
            }
        }
        else {
            /* ������ ����� ���������� */
            cont_ncons[R_GAS] = rr * pow( p1 / pr, 1.0 / paramsc->g2 );
            shr = vr + cr;
            if ( s >= shr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( p2 / pr, g1 );
               str = v2 + cmr;
               if ( s <= str ) {
                   /* ��������� ������ �� ����������� ������� */
                   if ( s <= v_cont_solid ) {
                       r = rr * pow( p1 / pr, 1.0 / paramsc->g2 );
                       v = v1;
                       p = p1;
                   }
                   else {
                       r = rr * pow( p2 / pr, 1.0 / paramsc->g2 );
                       v = v2;
                       p = p2;
                   }
               }
               else {
                    /* ��������� ������ ������ ����� ���������� */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = pr * pow( c / cr, g3 );
               }
            }
        }
    }
    
    /* ������������ ��������� ������� � ����������� */
    v_ncons_res[R_GAS] = r;
    v_ncons_res[V_GAS] = v;
    v_ncons_res[P_GAS] = p;
    
}

/* ������ ���������������� ������������ ������� "������"
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model
   of compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - ������� (30).
   
   v_solid_cont - �������� ����������� ������� � ���������� ���� (in)
   beta_l - �������� ���� ���������� ���� ����� �� ������� (in)
   beta_r - �������� ���� ���������� ���� ������ �� ������� (in)
   p_solid_l - �������� � ���������� ���� ����� �� ����������� ������� � ���������� ���� (in)
   p_solid_r - �������� � ���������� ���� ������ �� ����������� ������� � ���������� ���� (in)
   
   v_ncons_term[M] - ���������������� ������������ ������� "������" (out) */
void calc_ncons_term( double v_solid_cont, double beta_l, double beta_r, double p_solid_l, double p_solid_r,
                      double v_ncons_term[M] ) {

    v_ncons_term[B_DISP] = - v_solid_cont * ( beta_r - beta_l );    /* �������� � ��������� �������� �������� ���� ���������� ���� */
    v_ncons_term[R_DISP] = 0.0;                                     /* �������� � ��������� ������������� ��� ���������� ���� */
    v_ncons_term[V_DISP] = p_solid_r * beta_r - p_solid_l * beta_l; /* �������� � ��� ��� ���������� ���� */
    v_ncons_term[P_DISP] = v_solid_cont * v_ncons_term[V_DISP];         /* �������� � ��� ��� ���������� ���� */
    v_ncons_term[R_GAS] = 0.0;                                     /* �������� � ��������� ������������� ��� ������� ���� */
    v_ncons_term[V_GAS] = - v_ncons_term[V_DISP];                      /* �������� � ��� ��� ������� ���� */
    v_ncons_term[P_GAS] = - v_ncons_term[P_DISP];                      /* �������� � ��� ��� ������� ���� */

}

/* ������ ������-������� ��� ����������� ����������� �� ���������� ������� � ���������� ����, �� ������� �����, � ����� ���������
   ���������� �������� � ��������� ���� ����� � ������ �� ����������� ������� � ���������� ����.
   ������ ������ ���������� ���������� ���� ������ �� �������.
   
   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 501 - 502.

   ����� ������� ����������� ���������� � ����� ���������:
   - ������ ������ - ����� ����: 1 - ���������� ����, 2 - ������� ����
   - ������ ������ - ��������� ������������ ����������� �������: 1 - �����, 2 - ������
   ������� l � r ������������� ���������� � ������������� ����� ����� � ������ �� �������, ��������������.
   
   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
   right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
   c[K_GENERAL_CASE] - ������ ��������� ����� � ���������� � ������� ����� �� ������ ������� �� ������� (in)
   curr_p[K_GENERAL_CASE] - ������� ������ �������� ���������� � ������� ��� ����� � ������ �� ������� ���������� (in)

   curr_p[GAS_LEFT] - �������� ���� ����� �� ����������� ������� � ���������� ����
   curr_p[GAS_RIGHT] - �������� ���� ������ �� ����������� ������� � ���������� ����
   curr_p[DISP_LEFT] - �������� ���������� ���� ����� �� ����������� ������� � ���������� ����

   P[K_GENERAL_CASE] - ������� ������-������� (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - ������� ������� ����� (out)
   v_cont_disp - �������� ����������� ������� � ���������� ���� (out)
   v_cont_gas - �������� ����������� ������� � ������� ���� (out)
   v_gas_left - �������� ���� ����� �� ����������� ������� � ���������� ���� (out)
   v_gas_right - �������� ���� ������ �� ����������� ������� � ���������� ���� (out) */
void calc_P_and_DP_left_disp_phase( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                    double c[K_GENERAL_CASE], double curr_p[K_GENERAL_CASE], double P[M],
                                    double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                    double *v_gas_left, double *v_gas_right ) {

    /* ��������  */
                                        
    double p11 = curr_p[DISP_LEFT];                         /* ������� �������� � ���������� ���� ����� �� ����������� �������
                                                               � ���������� ���� */
    double p21 = curr_p[GAS_LEFT], p22 = curr_p[GAS_RIGHT]; /* ������� �������� � ������� ���� ����� � ������ �� ����������� �������
                                                               � ���������� ���� */
    
    /* �������� ����� */

    double c1l = c[DISP_LEFT];                      /* �������� ����� � ���������� ���� ����� �� ��������������� ������� */
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   /* �������� ����� � ������� ���� ����� � ������ �� ��������������� ������� */

    /* �������� ��� */

    double v1l = left_params[V_DISP];                           /* �������� ���������� ���� ����� �� ��������������� ������� */
    double v2l = left_params[V_GAS], v2r = right_params[V_GAS];   /* �������� ������� ���� ����� � ������ �� ��������������� ������� */

    double b1l = left_params[B_DISP];   /* �������� ���� ���������� ���� ����� �� ��������������� ������� */
    double b2l = 1.0 - b1l;         /* �������� ���� ������� ���� ����� �� ��������������� ������� */

    double f11, df11;  /* ������� � ����������� ��� ����������� �������� ���������� ���� �����
                           �� ����������� ������� � ���������� ���� */
    
    double f21, df21;  /* ������� � ����������� ��� ����������� �������� ������� ���� ����� �� ����������� ������� � ���������� ���� */
    double f22, df22;  /* ������� � ����������� ��� ����������� �������� ������� ���� ������ �� ����������� ������� � ���������� ���� */

    double g21, dg21;  /* ������� � ����������� ��� ����������� ��������� ������� ���� ����� �� ����������� ������� � ���� */
    double g22, dg22;  /* ������� � ����������� ��� ����������� ��������� ������� ���� ����� �� ����������� ������� � ���� */

    double v11;       /* ������� �������� ���������� ���� ����� �� ����������� ������� � ���������� ���� */
    double v21, v22;  /* ������� �������� ������� ���� ����� � ������ �� ����������� ������� � ���������� ���� */

    double r1, r2;  /* ��������� ���� ����� � ������ �� ����������� ������� � ���� */

    double g, A, B, C, D, E, dvL, dv1, dv2c;    /* ��������������� ���������� ��� ���������� ������ ���������� */

    /* �� �������� ������� ��������� ���������� �������� � ���������� � ������� ����� ����� ����������, ����� �� ���� ���������
       ������������ ����� ���������������. ������� - ���� �������� ���������� ��������. */

    /* ���������� ���� */
    calc_F_and_DF( paramsc, p11, left_params, M, c1l, DISPERSED_PHASE, &f11, &df11 );
    v11 = v1l - f11;
    /* ������� ���� */
    calc_F_and_DF( paramsc, p21, left_params, M, c2l, GAS_PHASE, &f21, &df21 );
    v21 = v2l - f21;
    calc_F_and_DF( paramsc, p22, right_params, M, c2r, GAS_PHASE, &f22, &df22 );
    v22 = v2r + f22;

    /* ������ ���������� ���� ����� � ������ �� ����������� ������� � ���� �� �������� �������� �� ���������� ������� */
    calc_G_and_DG( paramsc, p21, left_params, &g21, &dg21 );
    r1 = g21;
    calc_G_and_DG( paramsc, p22, right_params, &g22, &dg22 );
    r2 = g22;

    /* ������ ��������������� ���������� */
    g = paramsc->g2 / ( paramsc->g2 - 1.0 );
    A = pow( p22 / p21, 1.0 / paramsc->g2 );
    B = p22 - b1l * p11 - b2l * p21;
    dvL = v22 - v11;
    dv1 = v21 - v11;
    dv2c = v22 - v21;
    C = 0.5 * ( pow( dvL, 2.0 ) - pow( dv1, 2.0 ) );
    D = A * dvL;
    E = g / r1 / A;

    /* ������ ����������� ������-������� P */
    if ( v21 > v11 ) {
        P[0] = D - b2l * dv1;
        P[1] = B + b2l * r1 * dv1 * dv2c;
        P[2] = E * p22 - g * p21 / r1 + C;
        *v_cont_gas = v22;  /* ������� �������� ����������� ������� � ������� ���� */
    }
    else {
        /* ����� ��������, ���� 16 */
        P[0] = dvL - b2l * dv1 / A;
        P[1] = B + r2 * dvL * dv2c;
        P[2] = - g * p21 * A / r2 + g * p22 / r2 + C;
        *v_cont_gas = v21;  /* ������� �������� ����������� ������� � ������� ���� */
    }

    /* ������ ����������� ������� ����� */
    if ( v21 > v11 ) {
        /* ������ ������ */
        DP[0][0] = - D / p21 / paramsc->g2 + b2l * df21;
        DP[0][1] = D / p22 / paramsc->g2 + A * df22;
        DP[0][2] = ( A - b2l ) * df11;
        /* ������ ������ */
        DP[1][0] = b2l * ( - 1.0 + dv1 * dv2c * dg21 + r1 * df21 * ( dv1 - dv2c ) );
        DP[1][1] = 1.0 + b2l * r1 * dv1 * df22;
        DP[1][2] = - b1l + b2l * r1 * dv2c * df11;
        /* ������ ������ */
        DP[2][0] = g / r1 * ( p22 / A / paramsc->g2 / p21 - 1.0 ) + dv1 * df21 + dg21 * g / pow( r1, 2.0 ) * ( p21 - p22 / A );
        DP[2][1] = 1.0 / r1 / A + dvL * df22;
        DP[2][2] = dv2c * df11;
    }
    else {
        /* ����� ��������, ���� 16 */
        /* ������ ������ */
        DP[0][0] = b2l / A * ( df21 - dv1 / paramsc->g2 / p21 );
        DP[0][1] = df22 + b2l * dv1 / A / paramsc->g2 / p22;
        DP[0][2] = ( 1.0 - b2l / A ) * df11;
        /* ������ ������ */
        DP[1][0] = - b2l + r2 * dvL * df21;
        DP[1][1] = 1.0 + dg22 * dvL * dv2c + r2 * df22 * ( dv2c + dvL );
        DP[1][2] = - b1l + r2 * dv2c * df11;
        /* ������ ������ */
        DP[2][0] = - 1.0 / r2 * A + dv1 * df21;
        DP[2][1] = g / r2 + dvL * df22 - A * p21 / ( paramsc->g2 - 1.0 ) / r2 / p22 + dg22 * g / pow( r2, 2.0 ) * ( p21 * A - p22 );
        DP[2][2] = dv2c * df11;
    }

    *v_cont_disp = v11;     /* ������� �������� ����������� ������� � ���������� ���� */
    *v_gas_left = v21;      /* ������� �������� ���� ����� �� ����������� ������� � ���������� ���� */
    *v_gas_right = v22;     /* ������� �������� ���� ������ �� ����������� ������� � ���������� ���� */

}

/* ������ ������-������� ��� ����������� ����������� �� ���������� ������� � ���������� ����, �� ������� �����, � ����� ���������
   ���������� �������� � ��������� ���� ����� � ������ �� ����������� ������� � ���������� ����.
   ������ ������ ���������� ���������� ���� ����� �� �������.
   
   ������� ��������� � �����������.

   ����� ������� ����������� ���������� � ����� ���������:
   - ������ ������ - ����� ����: 1 - ���������� ����, 2 - ������� ����
   - ������ ������ - ��������� ������������ ����������� �������: 1 - �����, 2 - ������
   ������� l � r ������������� ���������� � ������������� ����� ����� � ������ �� �������, ��������������.
   
   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   left_params[M] - ������ ����������� ���������� ����� �� ������� (in)
   right_params[M] - ������ ����������� ���������� ������ �� ������� (in)
   c[K_GENERAL_CASE] - ������ ��������� ����� � ���������� � ������� ����� �� ������ ������� �� ������� (in)
   curr_p[K_GENERAL_CASE] - ������� ������ �������� ���������� � ������� ��� ����� � ������ �� ������� ���������� (in)
   
   curr_p[GAS_LEFT] - �������� ���� ����� �� ����������� ������� � ���������� ����
   curr_p[GAS_RIGHT] - �������� ���� ������ �� ����������� ������� � ���������� ����
   curr_p[DISP_RIGHT] - �������� ���������� ���� ������ �� ����������� ������� � ���������� ����

   P[K_GENERAL_CASE] - ������� ������-������� (out)
   DP[K_GENERAL_CASE][K_GENERAL_CASE] - ������� ������� ����� - ��������� ������ � ������� �� ����� ������ (out)
   v_cont_disp - �������� ����������� ������� � ���������� ���� (out)
   v_cont_gas - �������� ����������� ������� � ������� ���� (out)
   v_gas_left - �������� ���� ����� �� ����������� ������� � ���������� ���� (out)
   v_gas_right - �������� ���� ������ �� ����������� ������� � ���������� ���� (out) */
void calc_P_and_DP_right_disp_phase( struct ParametersCommon *paramsc, double left_params[M], double right_params[M],
                                     double c[M], double curr_p[M], double P[M],
                                     double DP[M][M], double *v_cont_disp, double *v_cont_gas,
                                     double *v_gas_left, double *v_gas_right ) {

    double p12 = curr_p[DISP_RIGHT];                        /* ������� �������� � ���������� ���� ������ �� ����������� �������
                                                               � ���������� ���� */
    double p21 = curr_p[GAS_LEFT], p22 = curr_p[GAS_RIGHT]; /* ������� �������� � ������� ���� ����� � ������ �� ����������� �������
                                                               � ���������� ���� */
    
    double c1r = c[DISP_RIGHT];                     /* �������� ����� � ���������� ���� ������ �� ��������������� ������� */
    double c2l = c[GAS_LEFT], c2r = c[GAS_RIGHT];   /* �������� ����� � ������� ���� ����� � ������ �� ��������������� ������� */

    double v1r = right_params[V_DISP];                          /* �������� ���������� ���� ������ �� ��������������� ������� */
    double v2l = left_params[V_GAS], v2r = right_params[V_GAS];   /* �������� ������� ���� ����� � ������ �� ��������������� ������� */

    double b1r = right_params[B_DISP];  /* �������� ���� ���������� ���� ������ �� ��������������� ������� */
    double b2r = 1.0 - b1r;         /* �������� ���� ������� ���� ����� �� ��������������� ������� */

    double f12, df12;   /* ������� � ����������� ��� ����������� �������� ���������� ���� ������ �� ����������� ������� � ���������� ���� */
    
    double f21, df21;   /* ������� � ����������� ��� ����������� �������� ������� ���� ����� �� ����������� ������� � ���������� ���� */
    double f22, df22;   /* ������� � ����������� ��� ����������� �������� ������� ���� ������ �� ����������� ������� � ���������� ���� */

    double g21, dg21;   /* ������� � ����������� ��� ����������� ��������� ������� ���� ����� �� ����������� ������� � ���� */
    double g22, dg22;   /* ������� � ����������� ��� ����������� ��������� ������� ���� ����� �� ����������� ������� � ���� */

    double v12;         /* ������� �������� ���������� ���� ������ �� ����������� ������� � ���������� ���� */
    double v21, v22;    /* ������� �������� ������� ���� ����� � ������ �� ����������� ������� � ���������� ���� */

    double r1, r2;  /* ��������� ���� ����� � ������ �� ����������� ������� � ���� */

    double g, A, B, C, D, E, dv2, dvR, dv2c;  /* ��������������� ���������� ��� ���������� ������ ���������� */

    // �� �������� ������� ��������� ���������� �������� � ���������� � ������� ����� ����� ����������, ����� �� ���� ���������
    // ������������ ����� ���������������. ������� - ���� �������� ���������� ��������.

    // ���������� ����
    calc_F_and_DF( paramsc, p12, right_params, M, c1r, DISPERSED_PHASE, &f12, &df12 );
    v12 = v1r + f12;
    // ������� ����
    calc_F_and_DF( paramsc, p21, left_params, M, c2l, GAS_PHASE, &f21, &df21 );
    v21 = v2l - f21;
    calc_F_and_DF( paramsc, p22, right_params, M, c2r, GAS_PHASE, &f22, &df22 );
    v22 = v2r + f22;

    // ������ ���������� ���� ����� � ������ �� ����������� ������� � ���� �� �������� �������� �� ���������� �������
    calc_G_and_DG( paramsc, p21, left_params, &g21, &dg21 );
    r1 = g21;
    calc_G_and_DG( paramsc, p22, right_params, &g22, &dg22 );
    r2 = g22;

    // ������ ��������������� ����������
    g = paramsc->g2 / ( paramsc->g2 - 1.0 );
    A = pow( p22 / p21, 1.0 / paramsc->g2 );
    B = b1r * p12 + b2r * p22 - p21;
    dv2 = v22 - v12;
    dvR = v21 - v12;
    dv2c = v22 - v21;
    C = 0.5 * ( pow( dv2, 2.0 ) - pow( dvR, 2.0 ) );
    D = A * dv2;
    E = g / r1 / A;

    // ������ ����������� ������-������� P
    if ( v21 > v12 ) {
        // ����� ��������, ���� 17
        P[0] = b2r * D - dvR;
        P[1] = B + r1 * dvR * dv2c;
        P[2] = E * p22 - g * p21 / r1 + C;
        *v_cont_gas = v22; // ������� �������� ����������� ������� � ������� ����
    }
    else {
        P[0] = b2r * dv2 - dvR / A;
        P[1] = B + b2r * r2 * dv2 * dv2c;
        P[2] = - g * p21 * A / r2 + g * p22 / r2 + C;
        *v_cont_gas = v21; // ������� �������� ����������� ������� � ������� ����
    }

    // ������ ����������� ������� �����
    if ( v21 > v12 ) {
        // ����� ��������, ���� 17
        // ������ ������
        DP[0][0] = - b2r * D / p21 / paramsc->g2 + df21;
        DP[0][1] = b2r * ( D / p22 / paramsc->g2 + A * df22 );
        DP[0][2] = ( 1.0 - b2r * A ) * df12;
        // ������ ������
        DP[1][0] = - 1.0 - r1 * dv2c * df21 + r1 * dvR * df21 + dv2c * dvR * df21;
        DP[1][1] = b2r + r1 * dvR * df22;
        DP[1][2] = b1r - r1 * dv2c * df12;
        // ������ ������
        DP[2][0] = g / r1 * ( p22 / A / paramsc->g2 / p21 - 1.0 ) + dg21 * g / pow( r1, 2.0 ) * ( p21 - p22 / A ) + dvR * df21;
        DP[2][1] = 1.0 / r1 / A + dv2 * df22;
        DP[2][2] = - dv2c * df12;
    }
    else {
        // ������ ������
        DP[0][0] = ( df21 - dvR / paramsc->g2 / p21 ) / A;
        DP[0][1] = b2r * df22 + dvR / A / paramsc->g2 / p22;
        DP[0][2] = ( 1.0 / A - b2r ) * df12;
        // ������ ������
        DP[1][0] = - 1.0 + b2r * r2 * dv2 * df21;
        DP[1][1] = b2r * ( 1.0 + dv2 * dv2c * dg22 + r2 * df22 * ( dv2c + dv2 ) );
        DP[1][2] = b1r - b2r * r2 * dv2c * df12;
        // ������ ������
        DP[2][0] = - ( paramsc->g2 + 1.0 ) / ( paramsc->g2 - 1.0 ) / r2 / A + dvR * df21;
        DP[2][1] = g / r2 * ( 1.0 + p21 / paramsc->g2 / p22 / A ) - g * dg22 / pow( r2, 2.0 ) * ( p22 - p21 / A ) + dv2 * df22;
        DP[2][2] = - dv2c * df12;
    }

    *v_cont_disp = v12; // ������� �������� ����������� ������� � ���������� ����
    *v_gas_left = v21; // ������� �������� ���� ����� �� ����������� ������� � ���������� ����
    *v_gas_right = v22; // ������� �������� ���� ������ �� ����������� ������� � ���������� ����

}

// ����� �������� �������������� ���������� ������� ��������� �����-�������� ��� ���������� ��������� �������� ����
// ���������� ����, ������ ������ ���� �� ������� � ����� ������ ��� ����� ����
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// center_params_full[M] - ������ ����������� ���������� � �������������� ������ ��� ������ ������� (in)
// left_minus_ncons[M] - ������������������ ������ ����������� ���������� �� ������ ������ ������ ����� �� �������������� ��� ������ ������� (in)
// left_plus_ncons[M] - ������������������ ������ ����������� ���������� �� ����� ����� �������������� ������ ��� ������ ������� (in)
// right_minus_ncons[M] - ������������������ ������ ����������� ���������� �� ������ ����� �������������� ������ ��� ������ ������� (in)
// right_plus_ncons[M] - ������������������ ������ ����������� ���������� �� ����� ����� ������ ������ �� �������������� ��� ������ ������� (in)
// phase - ������������� ����, ��� ������� ������������� ������� ��������� - ������� ��� ���������� (in)
// dt - ��������� ��� (in)
// h - ���������������� ��� (in)
// v_ncons_res (out)
// solution_full[M] - ������ ����������� ���������� � �������������� ������ ��� ������ ������� �� ��������� ���� (out)
// n - �������� ������ �������� center_params_full, left_minus_ncons_full, left_plus_ncons_full, right_minus_ncons_full, right_plus_ncons_full, solution_full
void godunov_classical_one_phase( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double center_params_full[M],
                                  double left_minus_ncons_full[M], double left_plus_ncons_full[M], double right_minus_ncons_full[M],
                                  double right_plus_ncons_full[M], Phase phase, double dt, double h, double v_ncons_res_left[M], 
                                  double v_ncons_res_right[M], double solution_full[M], int n ) {

    double left_flux_reduced[M_REDUCTION]; // ���������� ����� ����� ����� ����� ������
    double right_flux_reduced[M_REDUCTION]; // ���������� ����� ����� ������ ����� ������
    double left_params_reduced[M_REDUCTION]; // ���������� ������ ����������� ���������� � ������ ����� �� ���������������� �������
    double right_params_reduced[M_REDUCTION]; // ���������� ������ ����������� ���������� � ������ ������ �� ���������������� �������
    double center_ncons_params_reduced[M_REDUCTION]; // ���������� ������ ����������� ���������� � �������������� ������
    double center_cons_params_reduced[M_REDUCTION]; // ���������� ������ �������������� ���������� � �������������� ������
    double solution_ncons_reduced[M_REDUCTION]; // ����������� ���������� ������ ����������� ���������� � �������������� ������
    double solution_cons_reduced[M_REDUCTION]; // ����������� ���������� ������ �������������� ���������� � �������������� ������
    int i_comp; // ������ ���������� �������
    double cont_red[M_REDUCTION];
    // ������ ����������� ������ ����� ����� �����
    // ���������� ����� ��������� � ���������� �����������
    debug_info->neighbour_cell = debug_info->current_cell - 1;
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->current_cell_vncons[i_comp] = left_plus_ncons_full[i_comp];
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->neighbour_cell_vncons[i_comp] = left_minus_ncons_full[i_comp];
    // ������������ ����������� ���������� �������� �� ������ ����������
    convert_full_to_reduced( left_minus_ncons_full, phase, left_params_reduced );
    convert_full_to_reduced( left_plus_ncons_full, phase, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, phase, v_ncons_res_left, left_flux_reduced,cont_red );

    // ������ ����������� ������ ����� ������ �����
    // ���������� ����� ��������� � ���������� �����������
    debug_info->neighbour_cell = debug_info->current_cell + 1;
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->current_cell_vncons[i_comp] = right_minus_ncons_full[i_comp];
    for ( i_comp = 0; i_comp < n; i_comp++ )
        debug_info->neighbour_cell_vncons[i_comp] = right_plus_ncons_full[i_comp];
    // ������������ ����������� ���������� �������� �� ������ ����������
    convert_full_to_reduced( right_minus_ncons_full, phase, left_params_reduced );
    convert_full_to_reduced( right_plus_ncons_full, phase, right_params_reduced );
    godunov_flux_classical( paramsc, debug_info, left_params_reduced, right_params_reduced, phase, v_ncons_res_right, right_flux_reduced,cont_red );
    
    // ���������� ���������� � ������
    convert_full_to_reduced( center_params_full, phase, center_ncons_params_reduced );
    convert_noncons_to_cons_reduction( paramsc, center_ncons_params_reduced, phase, center_cons_params_reduced );
    for ( i_comp = 0; i_comp < M_REDUCTION; i_comp++ ) {
        solution_cons_reduced[i_comp] = center_cons_params_reduced[i_comp] - dt * ( right_flux_reduced[i_comp] - left_flux_reduced[i_comp] ) / h;
    }
    convert_cons_to_noncons_reduction( paramsc, solution_cons_reduced, phase, solution_ncons_reduced );

    // ������ ������� � ������ ������
    convert_reduced_to_full( solution_ncons_reduced, phase, solution_full );
    
}

// ������ ������������� ����������� ������ �.�. ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// left_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ����� �� ������� (in)
// right_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ������ �� ������� (in)
// phase - ������������� ����, ��� ������� �������������� ����� - ������� ��� ���������� (in)
// v_ncons_res[M_REDUCTION] - ������������ ������ ���������������� ���������� (in)
// flux[M_REDUCTION] - ������������ ������ ������ (out)
void godunov_flux_classical( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                             double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons_res[M_REDUCTION], double flux[M_REDUCTION], double cont_ncons_red[M_REDUCTION] ) {
    
    double cl, cr;                      // �������� ����� ����� � ������ �� �������
    double p_cont, v_cont;              // �������� � �������� �� ���������� �������
   
    cl = calc_sound_velocity_reduced( paramsc, left_ncons_params, phase );
    cr = calc_sound_velocity_reduced( paramsc, right_ncons_params, phase );

    // ������������ ��������� ������� �������� � �������� ���� �� ���������� �������
    calc_contact_pressure_velocity( paramsc, debug_info, left_ncons_params, right_ncons_params, M_REDUCTION,
        cl, cr, phase, &p_cont, &v_cont );

    // ����� �������
    sample_reduced( paramsc, left_ncons_params, right_ncons_params, M_REDUCTION, cl, cr, phase,
        p_cont, v_cont, 0.0, v_ncons_res, cont_ncons_red );
    cont_ncons_red[P] = p_cont;
    cont_ncons_red[V] = v_cont;
    // ������ ������ �������� �� ������� ���������������� ����������
    diff_flux_ncons_reduced( paramsc, v_ncons_res, phase, flux );
        
}

// ���������� ������� ������������ ���������� ������ ������
// params - ��������� � ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// left_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ����� �� ������� (in)
// right_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ������ �� ������� (in)
// phase - ������������� ����, ��� ������� �������������� ����� - ������� ��� ���������� (in)
// v_ncons[M_REDUCTION] - ������������ ������-������� (out)
void get_classical_Riemann_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                                     double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] ) {

    double cl = calc_sound_velocity_reduced( paramsc, left_ncons_params, phase ); // �������� ����� ����� �� �������
    double cr = calc_sound_velocity_reduced( paramsc, right_ncons_params, phase ); // �������� ����� ������ �� �������
    double cont_red[M_REDUCTION];
    double p_cont, v_cont; // �������� � �������� �� ���������� �������
    // ������������ ��������� ������� �������� � �������� ���� �� ���������� �������
    calc_contact_pressure_velocity( paramsc, debug_info, left_ncons_params, right_ncons_params, M_REDUCTION,
        cl, cr, phase, &p_cont, &v_cont );

    // ����� �������
    sample_reduced( paramsc, left_ncons_params, right_ncons_params, M_REDUCTION, cl, cr, phase,
        p_cont, v_cont, 0.0, v_ncons,cont_red );

}

/* ������� ������ ������� � ����� �� ��� �� ������������ ����������� �������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
   + �������� �� ���������� ��������� ���������: Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann
   problem and a high-resolution Godunov method for a model of compressible two-phase flow // Journal of Computational Physics.
   - 2006. - V. 212. - P. 490 - 526. - ������� (9) - (11).

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons_l - ������ ���������������� ���������� ����� �� ������� (in)
   v_ncons_r - ������ ���������������� ���������� ������ �� ������� (in)
   vector_size - ������ ������� ���������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   phase - ����, ��� ������� ���������� ������� - ������� ��� ���������� (in)
   p_cont - �������� �� ���������� ������� (in)
   v_cont - �������� �� ���������� ������� (in)
   s - �������� x/t, ��� �������� ���������� ������� (in)

   v_ncons_res - ������ ���������������� ���������� � ����������� ������������ (out) */
void sample_reduced( struct ParametersCommon *paramsc, double *v_ncons_l, double *v_ncons_r, int vector_size, double cl, double cr,
                     Phase phase, double p_cont, double v_cont, double s, double *v_ncons_res, double cont_ncons_red[M_REDUCTION]  ) {

    double rl, vl, pl;                  /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;                  /* ����������� ���������� ������ �� ������� */
    double g;                           /* ���������� �������� */
    double p01;                          /* �������� � ��������� ��������� */
    double g1, g2, g3, g4, g5, g6, g7;  /* ��������������� ����������, ����������� �� ���������� ��������,
                                           � ������������ � Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics.
                                           - 2nd Edition. - Springer, 1999. - P. 153. */

    /* �������� ����� ���� */
    double shl, stl;    /* �������� "������" � "������" ����� ����� ���������� */
    double sl;          /* �������� ����� ������� ����� */

    /* �������� ������ ���� */
    double shr, str;    /* �������� "������" � "������" ������ ����� ���������� */
    double sr;          /* �������� ������ ������� ����� */

    double cml, cmr;    /* �������� ����� ����� � ������ �� ����������� ������� */
    double c;           /* ��������� �������� ����� ������ ����� ���������� */
    double p_ratio;
    double r, v, p;     /* ���������� �������� �������� ����, ���������, �������� � �������� */

    /* ��������������� ���������� */
    /* ��������� ����� �� ������� */
    rl = v_ncons_l[R];
    vl = v_ncons_l[V];
    pl = v_ncons_l[P];
    /* ��������� ������ �� ������� */
    rr = v_ncons_r[R];
    vr = v_ncons_r[V];
    pr = v_ncons_r[P];
    switch ( phase ) {
        case GAS_PHASE:
            g = paramsc->g2;
            p01 = 0.0;
            break;
        case DISPERSED_PHASE:
            g = paramsc->g1;
            p01 = paramsc->p01;
            break;
        default:
            printf( "\nsample_reduced -> wrong phase identifier.\n\n" );
            exit( EXIT_FAILURE );
    }

    /* ����������� �� ���������� �������� */
    g1 = 0.5 * ( g - 1.0 ) / g;
    g2 = 0.5 * ( g + 1.0 ) / g;
    g3 = 2.0 * g / ( g - 1.0 );
    g4 = 2.0 / ( g - 1.0 );
    g5 = 2.0 / ( g + 1.0 );
    g6 = ( g - 1.0 ) / ( g + 1.0 );
    g7 = 0.5 * ( g - 1.0 );

    if ( s <= v_cont ) {
        /* ��������������� ����� - ����� �� ����������� ������� */
        if ( p_cont <= pl ) {
            /* ����� ����� ���������� */
            cont_ncons_red[R] = rl * pow(( p_cont + p01) / (pl + p01), 1.0 / g );
            shl = vl - cl;
            if ( s <= shl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( ( p_cont + p01 ) / ( pl + p01 ), g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* ��������� ����� �� ����������� ������� */
                    r = rl * pow( ( p_cont + p01 ) / ( pl + p01 ), 1.0 / g );
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* ��������� ������ ����� ����� ���������� */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = ( pl + p01 ) * pow( c / cl, g3 ) - p01;
                }
            }
        }
        else {
            /* ����� ������� ����� */
            p_ratio = ( p_cont + p01 ) / ( pl + p01 );
            cont_ncons_red[R] = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* ��������� �� ����� ������� ������ */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        /* ��������������� ����� - ������ �� ����������� ������� */
        if ( p_cont > pr ) {
            /* ������ ������� ����� */
            p_ratio = ( p_cont + p01 ) / ( pr + p01 );
            cont_ncons_red[R] = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* ��������� �� ������ ������� ������ */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* ������ ����� ���������� */
            cont_ncons_red[R] = rr * pow(( p_cont + p01) / (pr + p01), 1.0 / g );
            shr = vr + cr;
            if ( s >= shr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( ( p_cont + p01 ) / ( pr + p01 ), g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* ��������� ������ �� ����������� ������� */
                   r = rr * pow( ( p_cont + p01 ) / ( pr + p01 ), 1.0 / g );
                   v = v_cont;
                   p = p_cont;
               }
               else {
                    /* ��������� ������ ������ ����� ���������� */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = ( pr + p01 ) * pow( c / cr, g3 ) - p01;
               }
            }
        }
    }
    
    /* ������������ ��������� ������� � ����������� */
    v_ncons_res[R] = r;
    v_ncons_res[V] = v;
    v_ncons_res[P] = p;
    
}

/* ������� ��� ���������� ������� �����, ��������������� ��������� ���������, ���� ����� ��������� ���������� �����
   ����� ��� ������ �� �������, ������ �� ��� ������� ����� ��� ����� ����������

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons - ������ ���������������� ���������� ����� ��� ������ �� ������� (in)
   vector_size - ������ ������� ���������� (in)
   c - �������� ����� ����� ��� ������ �� ������� (in)
   phase - ����, ��� ������� �������� �������� - ������� ��� ���������� (in)
   dir - �����������, ��� �������� �������� �������� - ����� ��� ������ �� ������� (in) */
void draw_adiabatic_curve( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *v_ncons, int vector_size, double c, Phase phase, Direction dir ) {

    int i_pt;                       /* ������� ���������� ����� */
    /* ��� ������ �� ��� �������� (�������) */
    double h = ( ADAIABATIC_CURVE_P_T - ADAIABATIC_CURVE_P_B ) / ( ADIABATIC_CURVE_PNUM - 1 );
    double curr_p;                  /* ������� �������� */
    double curr_v;                  /* ������� �������� */
    int v_index;                    /* ������ ����������, ��������������� �������� ��������������� ���� */
    double f, fd;                   /* �������� ������� � ����������� */
    char fname[MAX_STRING_SIZE];    /* ��� ����� ��� ������ ������ ��� ���������� ������� */
    FILE *adiab_curve;              /* ���������� ����� */    

    /* ������������ ����� ����� � �������� ����� ��� ������ ������ ��� ���������� ������� */
    strcpy_s( fname, paramsc->output_file_directory );
    if ( dir == LEFT ) {
        strcat_s( fname, "\\adiabatic_curve_left.dat" );
    }
    else {
        strcat_s( fname, "\\adiabatic_curve_right.dat" );
    }
    if ( ( fopen_s( &adiab_curve, fname, "wt" ) ) != 0 ) {
        printf( "\ndraw_adiabatic_curve -> can't open file %s for writing\n\n", fname );
        exit( EXIT_FAILURE );
    }
    
    /* ����������� ������� ����������, ��������������� �������� */
    if ( vector_size == M_REDUCTION ) {
        v_index = V;
    }
    else {
        if ( phase == GAS_PHASE )
            v_index = V_GAS;
        else
            v_index = V_DISP;
    }

    for ( i_pt = 0; i_pt < ADIABATIC_CURVE_PNUM; i_pt++ ) {
        curr_p = ADAIABATIC_CURVE_P_B + i_pt * h;   /* ������� �������� �������� */
        calc_F_and_DF( paramsc, curr_p, v_ncons, vector_size, c, phase, &f, &fd );
        
        if ( dir == LEFT ) curr_v = v_ncons[v_index] - f;   
        else curr_v = v_ncons[v_index] + f;
        
        fprintf_s( adiab_curve, "%e %e\n", curr_v, curr_p );
    }

    fclose( adiab_curve );

}


// ���������� �������� �������� ���� ������� ����, ������������ � full_decouple_case_flux (������ NO_GRAD) ��� ������� ������
// case_beta - ����� ���������������� �������� ( 0 - ������� �������� � ������ ������, 1 - ������� �������� ����� 
//             ��� ������ �� �����, ����� ������� ��������� �����, � ����������� �� ����� �������� ���������� ����) (in)
// left_edge_beta - �������� �������� ���� ���������� ���� ����� �� �����, ����� ������� ��������� ����� (in)
// center_beta - �������� �������� ���� ���������� ���� � ������ �������������� ������ (in)
// right_edge_beta - �������� �������� ���� ���������� ���� ������ �� �����, ����� ������� ��������� ����� (in)
// solid_velocity - �������� �������� ���������� ����, ���������� � ���������� ������� ������ ������ � ������� ������� (in)
double full_decouple_case_flux_volume_fraction ( int case_beta, double left_edge_beta, double center_beta, double right_edge_beta, double solid_velocity ){

    switch ( case_beta ){
    
    case 0:
        return center_beta;
        break;
    case 1:
        if ( solid_velocity > 0.0 ){
            return left_edge_beta;
            break;
        }
        else{
            return right_edge_beta;
            break;
        }
    default:
        printf( "\nfull_decouple_case_flux_volume_fraction -> wrong case_beta.\n\n" );
        system ( "Pause" );
        exit( EXIT_FAILURE );
    }

}