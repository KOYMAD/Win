/*
 * cir_1.h
 *
 * ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������. ������ �����������.
 *
 * ����������� �.�., ��������� �.�., ������� �.�. �������������� ������� ���������� ������� ��������������� ������
 * ���������. - �.: ���������, 2001. - �. 62. - ������� (2.3.18), (2.3.19).
 *
 * (c) ����� �����, 2013
 *
 * ������: 24 ��� 2013 �.
 *
 */

#ifndef __CIR1_H_
#define __CIR1_H_

#include "struct.h"

// ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������, ������ �����������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// left_params[M] - ������ ����������� ���������� � ������ ����� �� �������������� (in)
// center_params[M] - ������ ����������� ���������� � �������������� ������ (in)
// right_params[M] - ������ ����������� ���������� � ������ ������ �� �������������� (in)
// dt - ��������� ��� (in)
// h - ���������������� ��� (in)
// solution[M] - ������ ����������� ���������� � �������������� ������ �� ��������� ���� (out)
void cir_1( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M],
            double dt, double h, double solution[M] );

// ������ �������, �������� � ������������� �������-������� �� �����
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// ncons_params[M] - ������ ����������� ���������� (in)
// S[M][M] - ������� ������� (out)
void cir_util_cir1( struct ParametersCommon *paramsc, double ncons_params[M], double S[M][M] );

/*  ������ ������������� �������-������� �� ����� ������
    left_params[M] - ������ ����������� ���������� � ������ ����� �� ����� (in)
    right_params[M] - ������ ����������� ���������� � ������ ������ �� ����� (in)
    S[M][M] - �������, ������������ �� ���������� � �������������� ������ (in)
    edge_sol[M] - ������������� �������-������� � ����������� ���������� �� ����� (out) */
void calc_edge_solution_cir1( double left_params[M], double right_params[M], double S[M][M], double edge_sol[M] );

/*  ������ "������" ����� ����� ������
    edge_sol[M] - ������������� �������-������� � ����������� ���������� �� ����� (in)
    A_ncons[M][M] - ������� ������� � ����������� ���������� � �������������� ������ (in)
    flux[M] - "�����" ����� ����� ������ (out) */
void calc_flux_cir1( double edge_sol[M], double A_ncons[M][M], double flux[M] );

#endif /* __CIR1_H_ */