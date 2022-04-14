// eos.h
// ��������� ��������� ���
// (c) ����� �����, 2015
// ������: 1 ������� 2015 �.

#ifndef __EOS_H
#define __EOS_H

#include "struct.h"

// ���������� ���������� ������� ������� ���� ��� ������� �������� � ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// p - ��������
// r - ���������
double e_gas( const struct ParametersCommon* paramsc, const double p, const double r );

// ���������� �������� ������� ���� ��� ������� ���������� ������� � ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// e - ���������� �������
// r - ���������
double p_gas( const struct ParametersCommon* paramsc, const double e, const double r );

// ���������� ���������� ������� ���������� ���� ��� ������� �������� � ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// p - ��������
// r - ���������
double e_disp( const struct ParametersCommon* paramsc, const double p, const double r );

// ���������� �������� ���������� ���� ��� ������� ���������� ������� � ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// e - ���������� �������
// r - ���������
double p_disp( const struct ParametersCommon* paramsc, const double e, const double r );

// ���������� ����������� ���������� ���� ��� ������� �������� � ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// p - ��������
// r - ���������
double T_disp( const struct ParametersCommon* paramsc, const double p, const double r );

// ���������� ����������� ������� ���� ��� ������� �������� � ���������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// p - ��������
// r - ���������
double T_gas( const struct ParametersCommon* paramsc, const double p, const double r );

#endif // __EOS_H