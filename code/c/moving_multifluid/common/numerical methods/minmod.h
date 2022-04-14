// minmod.h
// ��������� ������� ������������� ������ �� ���� �������-��������� ����������� �������� �������������� ����������
// � �������������� ������������ minmod.
// (c) ����� �����, 2016
// ������: 22 ������� 2016 �.

#ifndef __MINMOD_H_
#define __MINMOD_H_

#include "utils.h"
#include "math_utils.h"

// ������������� �������������� ���������� � �������������� ������������ minmod
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
// u_prev - ������� ����������� ���������� �� ���� ������� ��������� �������
// xc - ������ ��������� ������� ����� �����
// slopes - ������� �������� �� ���� ���������� ������� ��������� �������
void reconstruction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, double **u_prev, double *xc, double **slopes, int number_of_scalars );

#endif // __MINMOD_H_