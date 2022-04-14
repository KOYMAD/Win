/*
 * exact_solution_1d2phc.h
 *
 * ���������� ������� ������� ������ � ������� ������� ��� ������� ��������� ����� - ��������.
 *
 * (c) ����� �����, 2013
 *
 * ������: 24 ������� 2013 �.
 *
 */

#ifndef __EXACT_SOLUTION_1D2PHC_H
#define __EXACT_SOLUTION_1D2PHC_H

#include "godunov_bn.h"
#include "utils.h"
#include "memory.h"
#include "grid.h"

// ������ ������� ���������� �������� �� ������� [-0.5;0.5], ����� ����������, ����� �� ��������� ��������� �������
#define LEFT_BOUN_EX_SOL    -0.5
#define RIGHT_BOUN_EX_SOL   0.5

// ���������� ������� ������� ������ � ������� ������� ��� ������� ��������� �����-��������
// output_directory - ����������, ���� ����� ������� ���� � ������ �������� (in) 
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ (in)
// debug_info - ��������� � ���������� ����������� (in)
void build_exact_sol( char *output_directory, struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, int n );

#endif /* __EXACT_SOLUTION_1D2PHC_H */