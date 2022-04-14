#ifndef __SOURCE_TERMS_H_
#define __SOURCE_TERMS_H_

#include "eos.h"

// ���������� � ���� ������� ��������� ��������� � ���������� �����
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ��������������� ����������� ��� 1d2phc ������ (in)
// params2d - ��������� � ��������������� ����������� ��� 2d2phc ������ (in)
// number_of_block[2] - ������ � ������� �����. � 1d ������ ����� �������� ������ number_of_block[0] (in)
// v_ncons[M] - ������ ��������� ����������� ���������� � ����� (out)
void calc_initial_parameters_values(const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, int *number_of_block, double v_ncons[M]);

// ������ ���� ���������� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// v_ncons[M] - ������� ������ ����������� ���������� � ������ (in)
// ���������� ���� ������
double calc_friction_force( struct ParametersCommon *paramsc, double v_ncons[M] );

// ������ ������������ �������������
// Houim R.W., Oran E.S. A multiphase model for compressible granular-gaseous flows: formulation and initial tests // J. Fluid Mech. - 2016. - V. 789. - P. 166 - 220.
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// v_ncons[M] - ������� ������ ����������� ���������� � ������ (in)
// ���������� ����������� �������������
double calc_Cd( struct ParametersCommon *paramsc, double v_ncons[M] );

// ������ �������� ���������������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ��������������� ����������� ��� 1d2phc ������ (in)
// v_ncons[M] - ������� ������ ����������� ���������� � ������ (in)
// ���������� �������� ���������������
double calc_compaction_rate( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block );

// ������ configuration pressure Beta
// params1d - ��������� � ��������������� ����������� ��� 1d2phc ������ (in)
// v_ncons[M] - ������� ������ ����������� ���������� � ������ (in)
// ���������� configuration pressure Beta
double calc_configuration_pressure ( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block );

double compaction_energy( const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, const struct Parameters2d *params2d, Direction2d dir, double v_ncons[M], int *number_of_block );

// ������ ���������� �������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// v_ncons[M] - ������ ����������� ���������� � ������ ������ (in)
double calc_chemical_reaction(const struct ParametersCommon *paramsc, double v_ncons[M], int ignition_flag);

// ������ �������������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// v_ncons[M] - ������ ����������� ���������� � ������ ������ (in)
double calc_heat_transfer(const struct ParametersCommon *paramsc, double v_ncons[M]);

#endif // __SOURCE_TERMS_H_