/*
 * godunov.h
 *
 * ����� �������� ������� ������� ��� ���������� ������� ��������� �����-��������.
 *
 * Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of
 * compressible two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526.
 * 
 * (c) ����� �����, 2013
 *
 * ������: 4 ���� 2013 �.
 *
 */

#ifndef __GODUNOV_H_
#define __GODUNOV_H_

#include <string.h>
#include "utils.h"
#include "relaxation.h"


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
              double dt, double h, double solution[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, double curr_time, double *configuration_pressure );

// ������ �������� �������� �������� ���� ���������� ���� ����� � ������ �� �������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// left_beta - �������� �������� ���� ���������� ���� ����� �� ������� (in)
// right_beta - �������� �������� ���� ���������� ���� ������ �� ������� (in)
// ���������� ���� �� �������� ������������ Disp_phase_cases
Disp_phase_cases what_case( struct ParametersCommon *paramsc, double left_beta, double right_beta );


// ���� � ������ ���������� ����������� ���������� ����, �� ���������� � ��� ������� ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// solution[M] - ������ ����������� ���������� � �������������� ������ (in/out)
void set_background_state( struct ParametersCommon *paramsc, double solution[M] );

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
                          Disp_phase_cases solver_part, Direction edge, double curr_cell_beta, double godunov_flux[M], int n ) ;


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
void full_decouple_case_sol( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double center_params[M], double left_minus_ncons[M], double left_plus_ncons[M],
                             double right_minus_ncons[M], double right_plus_ncons[M], double dt, double h, double solution[M], int n ) ;


/* �������� ��������� ������� ��������� � ���������� ����������� ������� �������, �� �����
   ������������ ������ "�����", ��������� � �������� �������� ���� �������� � ������� ������.

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   debug_info - ��������� � ���������� ����������� (in)
   left_params_full[M] - ������ ������ ����������� ���������� ����� �� ������� (in)
   right_params_full[M] - ������ ������ ����������� ���������� ������ �� ������� (in)
   beta - �������� �������� ���� � �������������� ������ (in)
   
   flux_full[M] - ������ "�����" (out) */
void full_decouple_case_flux( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_params_full[M], double right_params_full[M],
                              double beta, double flux_full[M] );


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
                               double *p_solid_cont_l, double *p_solid_cont_r, int n, double cont_ncons[M] );

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
                                    double *p_cont_gas, double *v_cont_gas, double *p_cont_disp, double *v_cont_disp );

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
                      double solid_discontinuity_pressures[K_GENERAL_CASE], double *v1, double *v2, double *v21, double *v22 );

/* ������� ������� ���������� ���������������� ������� �� ���������� ����� ��� ����� � ������ �� �������

   c2l - �������� ����� � ������� ���� ����� �� ������� (in)
   c2r - �������� ����� � ������� ���� ������ �� ������� (in)
   c1l - �������� ����� � ���������� ���� ����� �� ������� (in)
   c1r - �������� ����� � ���������� ���� ������ �� ������� (in)

   sound_velocities[K_GENERAL_CASE] - ������ ��������� ����� (out) */
void fill_sound_velocities( double c2l, double c2r, double c1l, double c1r, double sound_velocities[K_GENERAL_CASE] );

/* ������� ������� ������������� ������� �������� ����� � ������ �� ������� � ���������� ����

   p_cont_gas - �������� ���� �� ���������� ������� ��� ����� ������� �������� ���� ���������� ���� (in)
   p_cont_solid - �������� ���������� ���� �� ���������� ������� ��� ����� ������� �������� ���� ���������� ���� (in)

   solid_discontinuity_pressures[K_GENERAL_CASE] - ������ ��������� �������� ����� � ������ �� ������� �������� ����
   ���������� ���� (out) */
void init_solid_discontinuity_pressures( double p_cont_gas, double p_cont_solid,
                                         double solid_discontinuity_pressures[K_GENERAL_CASE] );

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
                                    int vector_size, double cl, double cr, Phase phase, double *p_cont, double *v_cont );

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
                               double cl, double cr, int phase );

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
void calc_F_and_DF( struct ParametersCommon *paramsc, double curr_press, double *v_ncons, int vector_size, double c, int phase, double *F, double *DF );

/* ������ ������� G, ������������ ��������� ���� �� ������ ������� �� ����������� ������� � ����� ��� ������� �������� ���� ���������� ����,
   � �� ����������� �� �������� DG

   Schwendeman D.W., Wahle C.W., Kapila A.K. The Riemann problem and a high-resolution Godunov method for a model of compressible
   two-phase flow // Journal of Computational Physics. - 2006. - V. 212. - P. 490 - 526. - ������� G ������������ ��������� (7).

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   curr_press - �������� � ���������� �������� (in)
   v_ncons[M] - ������ ����������� ���������� (in)
   
   G - �������� ������� (out)
   DG - �������� ����������� (out) */
void calc_G_and_DG( struct ParametersCommon *paramsc, double curr_press, double v_ncons[M], double *G, double *DG );

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
                                        double *v21, double *v22 );

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
                                 double *v_gas_left, double *v_gas_right );

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
                            double p1, double p2, double v_cont, double s, double v_ncons_res[M], double cont_ncons[M] );

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
                          double p2, double v_cont_solid, double v_cont_gas, double v1, double v2, double s, double v_ncons_res[M], double cont_ncons[M] );

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
                      double v_ncons_term[M] );

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
                                    double *v_gas_left, double *v_gas_right );

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
                                     double *v_gas_left, double *v_gas_right );

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
                                  double v_ncons_res_right[M], double solution_full[M], int n );

// ������ ������������� ����������� ������ �.�. ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// left_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ����� �� ������� (in)
// right_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ������ �� ������� (in)
// phase - ������������� ����, ��� ������� �������������� ����� - ������� ��� ���������� (in)
// v_ncons_res[M_REDUCTION] - ������������ ������ ���������������� ���������� (in)
// flux[M_REDUCTION] - ������������ ������ ������ (out)
void godunov_flux_classical( struct ParametersCommon *paramsc, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                             double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons_res[M_REDUCTION], double flux[M_REDUCTION], double cont_ncons_red[M_REDUCTION] );

// ���������� ������� ������������ ���������� ������ ������
// params - ��������� � ����������� ��������������� ������������ (in)
// debug_info - ��������� � ���������� ����������� (in)
// left_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ����� �� ������� (in)
// right_ncons_params[M_REDUCTION] - ����������  ������ ����������� ���������� ������ �� ������� (in)
// phase - ������������� ����, ��� ������� �������������� ����� - ������� ��� ���������� (in)
// v_ncons[M_REDUCTION] - ������������ ������-������� (out)
void get_classical_Riemann_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons_params[M_REDUCTION],
                                     double right_ncons_params[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] );

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
                     Phase phase, double p_cont, double v_cont, double s, double *v_ncons_res, double cont_ncons_red[M_REDUCTION]  );

/* ������� ��� ���������� ������� �����, ��������������� ��������� ���������, ���� ����� ��������� ���������� �����
   ����� ��� ������ �� �������, ������ �� ��� ������� ����� ��� ����� ����������

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons - ������ ���������������� ���������� ����� ��� ������ �� ������� (in)
   vector_size - ������ ������� ���������� (in)
   c - �������� ����� ����� ��� ������ �� ������� (in)
   phase - ����, ��� ������� �������� �������� - ������� ��� ���������� (in)
   dir - �����������, ��� �������� �������� �������� - ����� ��� ������ �� ������� (in) */
void draw_adiabatic_curve( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *v_ncons, int vector_size, double c, Phase phase, Direction dir );

// ���������� �������� �������� ���� ������� ����, ������������ � full_decouple_case_flux (������ NO_GRAD) ��� ������� ������
// case_beta - ����� ���������������� �������� ( 0 - ������� �������� � ������ ������, 1 - ������� �������� ����� 
//             ��� ������ �� �����, ����� ������� ��������� �����, � ����������� �� ����� �������� ���������� ����) (in)
// left_edge_beta - �������� �������� ���� ���������� ���� ����� �� �����, ����� ������� ��������� ����� (in)
// center_beta - �������� �������� ���� ���������� ���� � ������ �������������� ������ (in)
// right_edge_beta - �������� �������� ���� ���������� ���� ������ �� �����, ����� ������� ��������� ����� (in)
// solid_velocity - �������� �������� ���������� ����, ���������� � ���������� ������� ������ ������ � ������� ������� (in)
double full_decouple_case_flux_volume_fraction ( int case_beta, double left_edge_beta, double center_beta, double right_edge_beta, double solid_velocity );

#endif /* __GODUNOV_H_ */