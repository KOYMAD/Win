// utils.cc
// ������� �������, ����������� ��� ���������������� ������ ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// (c) ����� �����, 2013 - 2015
// ������: 16 ������� 2013 �.

#ifndef __UTILS_1D2PHC_H_
#define __UTILS_1D2PHC_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "math_utils.h"
#include "io_1d.h"
#include "utils.h"


// ������ ���� �������������� �� �������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
// *x - ���������� ����� �����
// **v_ncons - ������� ����������� ���������� � ������� �����
// ���������� ��� �������������� �� �������
double get_time_step( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, const double* x, double **v_ncons, int print ) ;

// �������� ���������� ����� � �������
// params1d - ��������� � ����������� ���������� ������
// *x - ���������� ����� �����
// **v_ncons - ������� ����������� ���������� � ������� �����
// initial_total_mass - �������� ��������� ����� �������� � �������
// mass_diff - ������������� ����������� ��������� ����� �������� � ������� ������������ ��������
void check_mass( const struct Parameters1d* params1d, const double* x, double** v_ncons, const double initial_total_mass, double* mass_diff );

// ������ ������ ����� �������� � �������
// params1d - ��������� � ����������� ���������� ������
// *x - ���������� ����� �����
// **v_cons - ������� "��������������" ���������� � ������� �����
// ���������� ��������� ����� ���� � ���������� ���� � �������
double get_total_mass( const struct Parameters1d* params1d, const double* x, double** v_cons );

void calculate_DCJ(struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *Dcj);


#endif // __UTILS_1D2PHC_H_