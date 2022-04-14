// minmod.cc
// ��������� ������� ������������� ������ �� ���� �������-��������� ����������� �������� �������������� ����������
// � �������������� ������������ minmod.
// (c) ����� �����, 2016
// ������: 22 ������� 2016 �.

#include "minmod.h"

// ������������� �������������� ���������� � �������������� ������������ minmod
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
// u_prev - ������� ����������� ���������� �� ���� ������� ��������� �������
// xc - ������ ��������� ������� ����� �����
// slopes - ������� �������� �� ���� ���������� ������� ��������� �������
void reconstruction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, double **u_prev, double *xc, double **slopes, int number_of_scalars ) {

    for ( int i = 1; i < params1d->cells_number - 1; i++ ) { // ���� �� ���� ���������� ������� ��������� �������
        // �������������� � �������������� �����
        double u_left[M];
        double u_center[M];
        double u_right[M];
        convert_noncons_to_cons( paramsc, u_prev[i-1], u_left, number_of_scalars );
        convert_noncons_to_cons( paramsc, u_prev[i], u_center, number_of_scalars );
        convert_noncons_to_cons( paramsc, u_prev[i+1], u_right, number_of_scalars );
        for ( int j = 0; j < M; j++ ) { // ���� �� ����������� ������� �������������� ����������
            double slope_cand_1 = ( u_center[j] - u_left[j] ) / ( xc[i] - xc[i-1] ); // ����� �������� �� �������� �����������
            double slope_cand_2 = ( u_right[j] - u_center[j] ) / ( xc[i+1] - xc[i] ); // ������ �������� �� �������� �����������
            slopes[i][j] = 0.5 * ( sign( slope_cand_1, paramsc->eps_general ) + sign( slope_cand_2, paramsc->eps_general ) ) * 
                min( fabs( slope_cand_1 ), fabs( slope_cand_2 ) );
        }
    }
    // ��� �������� � ��������� ������� ������ ����������� ������� ���� - ������ �������
    for ( int j = 0; j < M; j++ ) { // ���� �� ����������� ������� �������������� ����������
        slopes[0][j] = 0.0;
        slopes[params1d->cells_number-1][j] = 0.0;
    }

}