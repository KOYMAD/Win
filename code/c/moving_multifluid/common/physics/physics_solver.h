/*
 * physics_solver.h
 *
 * ���� ���� "������" � ������, ������� ��������� �������������� � �������������� ����������.
 *
 * (c) ����� �����, 2014
 *
 * ������: 7 ������ 2014 �.
 *
 */

#ifndef __PHYSICS_SOLVER_H
#define __PHYSICS_SOLVER_H

#include "struct.h"

#include "memory.h"

// ������� ������� ������������ ���������������� ���������, ����������� ��������� �������������� � ���� "������"
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ��������������� ����������� 1dphc (in)
// dt - ������� ��� �� ������� (in)
// v_ncons[M] - ������ ����������� ���������� � ������, �� ������ - ���������� � ���������� ����� "������" (in/out)
void physics_solver( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double dt, double *v_ncons, double *v_ncons_next, double *v_ncons_prev, int *number_of_block, double *beta, double *B, double *C, int n, int ignition_flag, int initial_ignition, int cell_number, int left_ghost  ) ;

// ������ ������� ������ ������
// params� - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ��������������� ����������� 1dphc (in)
// v_ncons[M] - ������� ������ ����������� ���������� � ������ (in)
// right_hand_side_terms[M] - ������� ������ ������ ������ (out)
void calc_right_hand_side_terms( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double *v_ncons, double *right_hand_side_terms, int *number_of_block, int n, int ignition_flag, int initial_ignition ) ;

// ������� ���� ������� �����-����� (���������-���������) 2-�� �������
// params� - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ��������������� ����������� 1dphc (in)
// v_cons[M] - ������� ������ �������������� ���������� � ������ (in / out)
// right_hand_side_terms[M] - ������� ������ ������ ������ (in)
// sub_dt - ��� �� ������� (in)
// i_step_current - ����� �������� ���� �� ������������ (in)
void runge_kutta_predictor_2_nd_order (struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, Direction2d dir, double v_cons[M], double right_hand_side_terms[M], double sub_dt, int *number_of_block, int ignition_flag, int initial_ignition );

// ������� ������ ������� ��������� (7) �� ������� ����� ������ ������� ������� - ������ 1d2phc
// params� - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ��������������� ����������� 1dphc (in)
// v_cons[M1D] - ������� ������ �������������� ���������� � ������ (in / out)
// sub_dt - ��� �� ������� (in)
// number_of_block - ����� ����� (in)
// n - �������� ������ �������� (in)
void impicit_euler_full_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double v_cons[M], double sub_dt, int number_of_block, int n );

// ������� ������ ������� ��������� (9) �� ������� ����� ������ ������� ������� - ������ 2d2phc
// params� - ��������� � ��������� ����������� ��������������� ������������ (in)
// params2d - ��������� � ��������������� ����������� 2d2phc (in)
// v_cons[M] - ������� ������ �������������� ���������� � ������ (in / out)
// sub_dt - ��� �� ������� (in)
// number_of_block - ����� ����� (in)
// n - �������� ������ ��������
void impicit_euler_full_2d( struct ParametersCommon *paramsc, struct Parameters2d *params2d, double v_cons[M], double sub_dt, int *number_of_block, int n );

#endif /* __PHYSICS_SOLVER_H */