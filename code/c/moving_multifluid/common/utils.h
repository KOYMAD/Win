// utils.h
// ������� �������, ����������� ��� ���������������� ������ ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// ��� ����������� �� ����������� ������
// (c) ����� �����, 2018
// ������: 26 ������� 2018 �.

#ifndef __UTILS_H_
#define __UTILS_H_

#include "struct.h"
#include "utils_bn.h"
#include "io_1d.h"

// ��������� ���������� ����� ��� ������� ��������� ��
#include "cir_1_bn.h"
#include "cir_2_bn.h"
#include "cir_3_bn.h"
#include "cir_4_bn.h"
#include "godunov_bn.h"

// ������ ������� ��������� ��
#include "hll.h"
#include "rusanov_sa_1d2phc.h"
#include "hllc_sa.h"

void initiate_status( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double left, double right, int *status, int *left_ghost, int *right_ghost );

void calc_status (  struct ParametersCommon *paramsc, struct Parameters1d *params1d, double coordinate_of_left_boundary_of_body, double coordinate_of_right_boundary_of_body,
    int *status, int *index_of_left_ghost_cell, int *index_of_right_ghost_cell, double **v_ncons, double v_left[M], double v_right[M]);

void convert_2D_to_1D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector);

void convert_1D_to_2D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector);

// ������������� �������-�������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// files_directory - ����������, � ������� ��������� ��� ������� �����, ��������� ��� ������� (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out)
void init_solution_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, char *files_directory, struct TimeMoment *time_mom, double **initial_solution );

// ��������� ���������� �������
// params1d - ��������� � ����������� ���������� ������ (in)
// v_ncons[M] - ������ ����������� ���������� (in)
// boun_type - ��� ���������� ������� (in)
// boun_v[M] - ������ ����������� ���������� � ��������� ������ (out)
// curr_t - ������� ������ ������� (in)
// curr_inflow_parameters[M] - ������� �������� ���������� ��������(in)
void boundary( struct Parameters1d *params1d, struct Parameters2d *params2d, double v_ncons[M], int boun_type, double boun_v[M], double curr_t, double curr_inflow_parameters[M] ) ;

void current_inflow_parameters( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, int boun_direction, double *curr_inflow_params);

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
                        double dt, double h, double u_next[M], int step_number, int n, bool is_pressure_relaxation_now, int number_of_scalars, double curr_time, double *configuration_pressure, double body_velocity, int i, int *status, double cont_left[M], double cont_right[M] ) ;

// ������������� �������-�������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// files_directory - ����������, � ������� ��������� ��� ������� �����, ��������� ��� ������� (in)
// time_mom - ���������, ������������ ������� ������ ������� (out)
// **initial_solution - ������ �������� � ����������� ���������� - ��������� ������� � ������ ������ (out)
void init_solution_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, char *files_directory, struct TimeMoment *time_mom, double **initial_solution );

// ������������� ��������� ��� �������� � �������� ���������� ����������
// output_file_directory - ����������, � ������� ������ ���� �������� �����
// debug_info - ��������� � ���������� �����������
void init_debug_info( const char *output_file_directory, struct DebugInfo *debug_info );

// ������ ������������ ���������� ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ���������������� ���������� � ��������������� ������
// ���������� ������� ��������
double calc_p_i( const struct ParametersCommon* paramsc, const double v_ncons[M] );

// ������ ����������� ��������� ���������, ����� ����� ������ ��� ������ Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ���������������� ���������� � ��������������� ������
// ���������� ������� ���������
double calc_r_i( const struct ParametersCommon* paramsc, const double v_ncons[M] );

// ������ ��������� ����� ��� ����� ��� �� ������� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// c1 - �������� ����� � ���������� ����
// c2 - �������� ����� � ������� ����
void calc_sound_velocity( const struct ParametersCommon* paramsc, const double v_ncons[M], double* c1, double* c2 );

// ������ ��������� ����� ��� �������� ���� �� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// phase - ������������� ����, ��� ������� �������������� �������� ����� - ������� ��� ����������
// ���������� �������� �����
double calc_sound_velocity_one_phase( const struct ParametersCommon* paramsc, const double v_ncons[M], const Phase phase );

// ������ �������� ����� ��� ����� ���� �� ����������� ������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// phase - ������������� ����, ��� ������� �������������� �������� ����� - ������� ��� ����������
// ���������� �������� �����
double calc_sound_velocity_reduced( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase );

// ������ ����������� ��������� ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ���������������� ���������� � ��������������� ������
// ���������� ������� ��������
double calc_u_i( const struct ParametersCommon* paramsc, const double v_ncons[M] );

// �������������� ������� "��������������" ���������� � ������ �����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// v_ncons - ������ ����������� ����������
// number_of_scalars - ���������� ����������� ��������
// ����������: SUCCEESS            �������� � ����� ����� ������������
//             NEGATIVE_PRESSURE   �������� ���� �� � ����� �� ��� ������������
ReturnCodes convert_cons_to_noncons( const struct ParametersCommon* paramsc, const double v_cons[M], double v_ncons[M], int number_of_scalars  );

/* �������������� ��������������� ����������� ������� "��������������" ���������� � ������ �����������

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_cons[M_REDUCTION] - ������ "��������������" ���������� (in)
   phase - ������������� ����, ��� ������� �������������� �������������� - ������� ��� ���������� (in)

   v_ncons[M_REDUCTION] - ������ ����������� ���������� (out) */
void convert_cons_to_noncons_reduction( struct ParametersCommon *paramsc, double v_cons[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] );

/* ������������ ����������� ������� �� ������� ����������� �������
 
   full_vector[M] - ������ ���������� ������ (in)
   phase - ������������� ����, ��� ������� �� ������� ������� �������� �������������� - ������� ��� ���������� (in)
 
   reduced_vector[M_REDUCTION] - �������������� ���������� ������ (out) */
void convert_full_to_reduced( double full_vector[M], Phase phase, double reduced_vector[M_REDUCTION] );

// �������������� ������� ����������� ���������� � ������ "��������������"
// params - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// v_cons - ������ "��������������" ����������
// number_of_scalars - ���������� ����������� ��������
void convert_noncons_to_cons( const struct ParametersCommon *paramsc, const double v_ncons[M], double v_cons[M], int number_of_scalars );

/* �������������� ����������� ������� ����������� ���������� � ������ "��������������"

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons[M_REDUCTION] - ������ ����������� ���������� (in)
   phase - ������������� ����, ��� ������� �������������� �������������� - ������� ��� ���������� (in)

   v_cons[M_REDUCTION] - ������ "��������������" ���������� (out) */
void convert_noncons_to_cons_reduction( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase, double v_cons[M_REDUCTION] );

// ������������ ����� ������� ����������� ������� �� ������������ �����������
// reduced_vector[M_REDUCTION] - �������������� ���������� ������ (in)
// phase - ������������� ����, ��� ������� �� ��������������� ������� �������� ����� ������� - ������� ��� ���������� (in)
// full_vector[M] - ������ ���������� ������ (out)
void convert_reduced_to_full( double reduced_vector[M], Phase phase, double full_vector[M_REDUCTION] );

// ������ ������� ����������������� "������" �� ������� "��������������" ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// flux - �������������� ������ ����������������� "������"
void diff_flux_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* flux, int number_of_scalars ) ;

// ������ ������� ����������������� "������" �� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// flux - �������������� ������ ����������������� "������"
void diff_flux_ncons( const struct ParametersCommon* paramsc, const double v_ncons[M], array1D* flux, int number_of_scalars ) ;

/* ������ ����������� ������� ����������������� "������" �� ������������ ����������� ������� ����������� ����������

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons - ����������� ���������� ������ ����������� ���������� (in)
   phase - ���� - ������� ��� ����������, ��� ������� �������������� "�����" (in)

   flux - �������������� ����������� ���������� ������ ����������������� "������" (out) */
void diff_flux_ncons_reduced( struct ParametersCommon *paramsc, double *v_ncons, Phase phase, double *flux );

// ������ ����������������� ������� ������ ������ �� ������� "��������������" ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// rhst - ���������������� ������ ������ ������
void rhst_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* rhst, int number_of_scalars );

// ������ ����������������� ������� ������ ������ �� ������� ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_cons - ������ "��������������" ����������
// rhst - ���������������� ������ ������ ������
void rhst_ncons( const struct ParametersCommon* params�, const double v_ncons[M], array1D* rhst, int number_of_scalars );

// ���������� ������ ��������� ������� � ������
// debug_info - ��������� � ���������� �����������
void vector_debug_print( const struct DebugInfo *debug_info );

// ������ ������ ���������� � ���������� ����� � ������ ��������� ���������
// debug_info - ��������� � ���������� ����������� (in)
// index - ������ ����������� �������� � ������� �������� (in)
void debug_print( struct DebugInfo *debug_info, int index );

// ���������� ���������� ������ � ������������� ����
// paramsc - ��������� � ��������� ����������� ��������������� ������������, ����� ����� �������� ��������� ���������� � ������������� ���� (in/out)
// params1d - ��������� � ����������� ���������� ������, ����� ����� �������� ��������� ���������� � ������������� ���� (in/out)
void dimensionalization( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d ) ;

void current_block_number(struct ParametersCommon *params, struct Parameters1d *params1d, struct Parameters2d *params2d, int step_number_x, int step_number_y, int *number_of_block);

#endif // __UTILS_H_