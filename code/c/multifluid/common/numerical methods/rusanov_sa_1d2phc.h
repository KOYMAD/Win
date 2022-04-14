// rusanov_sa_1d2phc.h
// ����� �������� ���������� �������������� ��������� Saurel-Abgrall, 1D ������
// ����������� ��: Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) ����� �����, 2018
// ������: 4 ����� 2018 �.

#ifndef __RUSANOV_SA_1D2PHC_H_
#define __RUSANOV_SA_1D2PHC_H_

#include "relaxation.h"

// ����� �������� ���������� �������������� ��������� Saurel-Abgrall, 1D ������
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
// configuration_pressure - ���������������� ��������
void rusanov_1d( struct ParametersCommon* paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
                 double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                 double dt, double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure ) ;

// ��������������� �������� ������ ��������, 1D ������
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
void Lh_rusanov_1d( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
                    const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                    const double dt, const double h, double solution_ncons[M] );

// ����� �������� ������� ������� � ���������� �����, 1D ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� ����� �� �������
// right_ncons - ������ ����������� ���������� ������ �� �������
// flux - ������������ ������ ������
// smax - ������ ��� �������� S+
void rusanov_flux_1d( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M], array1D* flux, double *splus );

// ������ ������ ��� �������� ����� S+ � ������ ����� ���
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons[M] - ������ ����������� ���������� ����� �� �������
// right_ncons[M] - ������ ����������� ���������� ������ �� �������
// splus - ������ ��� �������� ����� S+
void calc_splus_rusanov( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M], double* splus );

#endif // __RUSANOV_SA_1D2PHC_H_