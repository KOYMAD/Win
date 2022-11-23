// hllc_sa.h
// ����� HLLC ���������� �������������� ��������� Saurel-Abgrall
// ����������� ��: Li Q. et al. Difference scheme for two-phase flow // Applied Mathematics and Mechanics. - 2004. - V. 25, No. 5. - P. 536 - 545
// �
// Liang S. et al. Solving seven-equation model for compressible two-phase flow using multiple GPUs // Computers & Fluids. - 2014. - V. 99. - P. 156 - 171
// (c) ����� �����, 2017
// ������: 29 ����� 2017 �.

#ifndef __HLLC_H_
#define __HLLC_H_

#include "relaxation.h"
#include "hll.h"

// ����� HLLC ���������� �������������� ��������� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� �����������
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// slopes_left - ������ �������� � ������ ����� �� ��������������
// slopes_center - ������ �������� � �������������� ������
// slopes_right - ������ �������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ �� ��������� ����
// n - �������� ������ ��������
void hllc_1d(struct ParametersCommon* paramsc, struct Parameters1d* params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
          double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure,  double body_velocity, int i, int *status,double cont_left[M], double cont_right[M] ) ;

// ��������������� �������� ������ HLLC
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// slopes_left - ������ �������� � ������ ����� �� ��������������
// slopes_center - ������ �������� � �������������� ������
// slopes_right - ������ �������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
void Lh_HLLC( struct DebugInfo *debug_info, struct Parameters1d* params1d,  struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars,  int i, int *status, double body_velocity,double cont_left[M], double cont_right[M]  );

// ����� Harten - Lax - van Leer - Contact (HLLC) ������� ������� � ���������� �����
// Liang S. et al. Solving seven-equation model for compressible two-phase flow using multiple GPUs //
// Computers & Fluids. - 2014. - V. 99. - P. 156 - 171.
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� ����� �� �������
// right_ncons - ������ ����������� ���������� ������ �� �������
// flux - ������������ ������ ������
// phi - �������� ���� ���������� ���� �� �����
// s_cont - ������ ��� �������� ����������� �������
void hllc_flux( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M], array1D* flux, double* phi, double* s_cont );

// ������ ������� - ������ ��������� ����� ��� ������ �� ����������� ������� � ������ HLLC
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// v_ncons - ������ ����������� ���������� ����� ��� ������ �� �������
// s - ������ �������� ����� ��� ������ �����
// s_star - ������ ��� s �� ����������
// q_star - ������ ��������� ����� ��� ������ �� ����������� �������
void calc_contact_vector( const struct ParametersCommon *paramsc, const double v_ncons[M], const double s, const double s_star, array1D* q_star );

#endif // __HLLC_H_