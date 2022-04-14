// relaxation.h
// ���������� ��������� � �������� ��� �� ��������� �������

// ���������� ��������� ����������� �� ������:
// Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467.

// ���������� �������� ����������� �� ������:
// ������ �.�. ��������� ������������� ����������� ������� � ������� ����������� ���������� ���� // ������� ���. � 2009. � �. 16, � 2. � �. 62 � 70.
// �������� �������� ������ � \science\utkin\docs\���������� ��������.docx

// (c) ����� �����, 2017
// ������: 29 ����� 2017 �.

#ifndef __RELAXATION_H_
#define __RELAXATION_H_

#include "struct.h"
#include "utils.h"
#include "eos.h"
#include "io.h"
#include "source_terms.h"

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
         const double right_ncons[M], const double dt, const double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure);

// �������� ���������� �������� ��� ������� ������� ��������� ���� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ��������������� ��������� ��������
void Lrv( const struct ParametersCommon* paramsc, const double center_ncons[M], const double dt, const double h, double solution_ncons[M], int n );

// �������� ���������� �������� ��� ������� ������� ��������� ���� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� �����������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ��������� � ��������� ���������� ��������
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� � ������� ��������������� ����������
void Lrp_Ivanov( const struct ParametersCommon* paramsc, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int n );

// �������� ���������� �������� ��� ������� ������� ��������� ���� Baer-Nunziato � ������ ��������������� Schwendeman
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� c ����������� ��������������� ������������, ��������� 1d ������
// debug_info - ��������� � ���������� �����������
// center_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ��������� � ��������� ���������� ��������
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� � ������� ��������������� ����������
// step_number - ����� ������� ������
void Lrp_compaction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, struct DebugInfo *debug_info, const double center_ncons[M], double solution_ncons[M], int step_number, int n, double *configuration_pressure);

#endif // __RELAXATION_H_