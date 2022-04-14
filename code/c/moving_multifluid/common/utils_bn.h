// utils_bn.h
// ������� �������, ����������� ��� ���������������� ������ ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// (c) ����� �����, 2013 - 2015
// ������: 16 ������� 2013 �.

#ifndef __UTILS_BN_H_
#define __UTILS_BN_H_

#include "struct.h"
#include "math_utils.h"

// ������ ������� ���������� ������� ��������� �����-�������� � "��������������" ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// v_cons[M] - ������ "��������������" ���������� (in)
// A_cons[M][M] - ������� ���������� ������� ��������� �����-�������� (out)
void calc_A_cons( struct ParametersCommon *paramsc, double v_cons[M], double A[M][M] );

// ������ ������� ���������� ������� ��������� �����-�������� � ����������� ����������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// v_ncons[M] - ������ ����������� ���������� (in)
// A_ncons[M][M] - ������� ���������� ������� ��������� �����-�������� � ����������� ���������� (out)
void calc_A_ncons( struct ParametersCommon *paramsc, double v_ncons[M], double A_ncons[M][M] );

void calc_omega_cons( const struct ParametersCommon* paramsc, const double v_cons[M], double omega[M][M] );

// ������ ������� �� ����������� �������� ���������� ������� ��������� �����-�������� � ����������� ����������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// v_ncons[M] - ������ ����������� ���������� (in)
// omega[M][M] - ������� �� ����������� �������� ���������� ������� ��������� �����-�������� � ����������� ���������� (out)
void calc_omega_ncons( struct ParametersCommon *paramsc, double v_ncons[M], double omega[M][M] );

// ������ �������, �������� ������� �� ����������� �������� ���������� ������� ��������� �����-��������
// � ����������� ����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// v_ncons[M] - ������ ����������� ���������� (in)
// omega_inverse[M][M] - �������, �������� ������� �� ����������� �������� ���������� ������� ��������� �����-��������
// � ����������� ���������� (out)
void calc_omega_ncons_inverse( struct ParametersCommon *paramsc, double v_ncons[M], double omega_inverse[M][M] );

// ������ ������������ ������� �� ������� ����������� ����� ������� ��������� �����-�������� �� ���������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// v_ncons[M] - ������ ����������� ���������� (in)
// sign_lambda[M][M] - ������������ ������� �� ������� ����������� ����� ������� ��������� �����-�������� �� ��������� (out)
void calc_sign_lambda( struct ParametersCommon *paramsc, double v_ncons[M], double sign_lambda[M][M] );

void calc_abs_lambda( const struct ParametersCommon* paramsc, const double v_cons[M], double lambda[M][M] );

/* ������ ������������ ������� � ������������ ������� ������� ��������� �����-�������� �� ���������

   paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
   v_ncons[M] - ������ ����������� ���������� (in)

   lambda[M][M] - ������������ ������� � ������������ ������� ������� ��������� �����-�������� �� ��������� (out) */
void calc_lambda( struct ParametersCommon *paramsc, double v_ncons[M], double lambda[M][M] );

// �������� ������������ �������� ���������������� ����������
// params - ��������� � ��������� ����������� ��������������� ������������ (in) 
// v_ncons[M] - ������ ���������������� ���������� (in)
void check_params_correctness( struct ParametersCommon *paramsc, double v_ncons[M] );

// ��������, ����� �� ������� ������� ������������ ������ �����, ������ � ����� � ����� ������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// v_ncons[M] - ������ ���������������� ���������� (in)
// ���������� true, ���� ����������� ����������� �����; false - �����
bool check_matrix_decomposition( struct ParametersCommon *paramsc, double v_ncons[M] );

#endif // __UTILS_BN_H_