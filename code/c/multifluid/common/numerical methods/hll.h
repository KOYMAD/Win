// hll.h
// ����� HLL �� Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) ����� �����, 2015
// ������: 8 ���� 2015 �.

#ifndef __HLL_H_
#define __HLL_H_

#include "struct.h"
#include "utils.h"
#include "eos.h"
#include "relaxation.h"

// ����� HLL ���������� �������������� ��������� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� ����������� (in)
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// slopes_left - ������ �������� � ������ ����� �� ��������������
// slopes_center - ������ �������� � �������������� ������
// slopes_right - ������ �������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ �� ��������� ����
// step_number - ����� �������� ���� �� ������������
// n - �������� ������ ��������
// is_pressure_relaxation_after_this_step - true, ���� ����� ������� ���������� ������ �� ������� ����������� ����� ��������� ���������� ��������; false - �����
//                              ���� � ������ ������ ������ �� ����� ��������� ���������� ��������, �� true �� �������� �� ���������, ��� ��� ������ ����� �������������� � parameters.dat
// number_of_scalars - ����� ����������� ��������
// curr_time - ������� ������ �������
// configuration_presure - ���������������� ��������
void hll( struct ParametersCommon* paramsc, struct Parameters1d* params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
          double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure, double S_center, double S_left, double S_right, int l );

// ��������������� �������� ������ HLL
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
// n - �������� ������ �������� ��� ����� ����������� ��������
// number_of_scalars - ���������� �������������� ���������, ��� �� ���������� ����������� ��������
void Lh_HLL( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars, double S_center,  double S_left, double S_right   );

// ����� Harten - Lax - van Leer (HLL) ������� ������� � ���������� �����
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� ����� �� �������
// right_ncons - ������ ����������� ���������� ������ �� �������
// flux - ������������ ������ ������
// splus - ������ ��� �������� ����� S+
// sminus - ������ ��� �������� ����� S-
// n - �������� ������ ��������
void hll_flux( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M],
               array1D* flux, double* splus, double* sminus, int n, int number_of_scalars  ) ;

// ������ ������ ��� ��������� ���� S+ � S-
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons[M] - ������ ����������� ���������� ����� �� �������
// right_ncons[M] - ������ ����������� ���������� ������ �� �������
// phase - ����, ��� ������� ����������� ��������
// splus - ������ ��� �������� ����� S+
// sminus - ������ ��� �������� ����� S-
void calc_splus_sminus( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M],
                        const Phase phase, double* splus, double* sminus );

// ������������� ������������� dB2 � ������ �����
// left_minus - ����� �����, ������ �������������� ���������� ����� �� �������
// left_plus - ����� �����, ������ �������������� ���������� ������ �� �������
// right_minus - ������ �����, ������ �������������� ���������� ����� �� �������
// right_plus - ������ �����, ������ �������������� ���������� ������ �� �������
// splus_l - ������ ��� �������� ����� S+ �� ����� �����
// sminus_l - ������ ��� �������� ����� S- �� ����� �����
// splus_r - ������ ��� �������� ����� S+ �� ������ �����
// sminus_l - ������ ��� �������� ����� S- �� ������ �����
// dB2dx - ������� �������������
void calc_dB2( const double left_minus[M], const double left_plus[M], const double right_minus[M], const double right_plus[M],
               double splus_l, double sminus_l, double splus_r, double sminus_r, double* dB2dx, double S_center, double S_left, double S_right );

/*
// �������������� �������� ������ HLL
// params - ��������� � ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ��������������� ���������
void Lr( const struct Parameters1d* params, const double left_ncons[M], const double center_ncons[M],
         const double right_ncons[M], const double dt, const double h, double solution_ncons[M] );

// �������� ���������� �������� ������ HLL
// params - ��������� � ����������� ��������������� ������������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ����������� ��������������� ��������� ��������
void Lrv( const struct Parameters1d* params, const double center_ncons[M], const double dt, const double h, double solution_ncons[M] );

// �������� ���������� �������� ������ HLL
// ---
// ���������� �� ������:
// ������ �.�. ��������� ������������� ����������� ������� � ������� ����������� ���������� ���� // ������� ���. � 2009. � �. 16, � 2. � �. 62 � 70.
// �������� �������� ������ � \science\utkin\docs\���������� ��������.docx
// ---
// params - ��������� � ����������� ��������������� ������������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ��������� � ��������� ���������� ��������
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� � ������� ��������������� ����������
void Lrp_Ivanov( const struct Parameters1d* params, const double center_ncons[M], double solution_ncons[M] );
*/
#endif // __HLL_H