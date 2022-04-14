// cir_4.cc
// ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������. ��������� �����������.
// ����������� �.�., ��������� �.�., ������� �.�. �������������� ������� ���������� ������� ��������������� ������
// ���������. - �.: ���������, 2001. - �. 67. - ������� (2.3.44), (2.3.45).
// (c) ����� �����, 2014
// ������: 08 ������� 2014 �.

#include "cir_4_bn.h"

#include "utils.h"
#include "math_utils.h"

// ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������, ��������� �����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_params - ������ ����������� ���������� � ������ ����� �� ��������������
// center_params - ������ ����������� ���������� � �������������� ������
// right_params - ������ ����������� ���������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution - ������ ����������� ���������� � �������������� ������ �� ��������� ����
void cir_4( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M],
            double dt, double h, double solution[M] ) {

    // "������" ����� ����� ������
    array1D left_flux( M );
    array1D right_flux( M );

    // ������ ���������������� ������������ �� ����� � ������ ������ ������
    array1D left_rhst( M );
    array1D right_rhst( M );

    // ��������� ��� � �������������� ����������
    double left_cons_params[M];
    convert_noncons_to_cons( paramsc, left_params, left_cons_params, 0 );
    double center_cons_params[M];
    convert_noncons_to_cons( paramsc, center_params, center_cons_params, 0 );
    double right_cons_params[M];
    convert_noncons_to_cons( paramsc, right_params, right_cons_params, 0 );

    // ������ "�������" ����� ����� � ������ ����� ������
    calc_flux_cir4( paramsc, left_cons_params, center_cons_params, &left_flux );
    calc_flux_cir4( paramsc, center_cons_params, right_cons_params, &right_flux );

    // ������ ������ ���������������� ������������ �� ����� � ������ ������
    calc_rhst_term_cir4( paramsc, left_cons_params, center_cons_params, &left_rhst );
    calc_rhst_term_cir4( paramsc, center_cons_params, right_cons_params, &right_rhst );

    double solution_cons[M];
    // ���������� ������� "��������������" ����������
    for ( int iComp = 0; iComp < M; iComp++ )
        solution_cons[iComp] = center_cons_params[iComp] - dt * ( right_flux[iComp] - left_flux[iComp] ) / h +
            0.5 * ( left_rhst[iComp] + right_rhst[iComp] ) / h;

    convert_cons_to_noncons( paramsc, solution_cons, solution, 0 );
}

// ������ "������" ����� ����� ������
// params - ��������� � ��������� ����������� ��������������� ������������
// left_params - ������ "��������������" ���������� ����� �� �������
// right_params - ������ "��������������" ���������� ������ �� �������
// flux - "�����" ����� ����� ������
void calc_flux_cir4( const struct ParametersCommon* paramsc, const double left_params[M], const double right_params[M], array1D* flux ) {

    vector<double> aver_flux( M ); // ����������� ���������������� "�����"
    vector<double> left_diff_flux( M ); // ������ ����������������� "������" �� ���������� ����� �� �������
    vector<double> right_diff_flux( M ); // ������ ����������������� "������" �� ���������� ������ �� �������
    vector<double> charact_part( M ); // ������������������ ������������ "������"
    double aver_params[M]; // ����������� ��������� �� �����
    double aver_ncons_params[M]; // ����������� ������ ����������� ���������� �� �����

    // ������ ����������� ���������� �� �����
    for ( int iComp = 0; iComp < M; iComp++ )
        aver_params[iComp] = 0.5 * ( left_params[iComp] + right_params[iComp] );
    
    // ������ ������, ���������������� ������������ "������"
    diff_flux_cons( paramsc, left_params, &left_diff_flux, 0 );
    diff_flux_cons( paramsc, right_params, &right_diff_flux, 0 );
    for ( int iComp = 0; iComp < M; iComp++ )   
        aver_flux[iComp] = 0.5 * ( left_diff_flux[iComp] + right_diff_flux[iComp] );

    // ������ ������, ������������������ ������������
    // ������ ������������� �� ������ ������������ ����� ������� �����
    convert_cons_to_noncons( paramsc, aver_params, aver_ncons_params, 0 );
    double c1, c2;
    calc_sound_velocity(paramsc, aver_ncons_params, &c1, &c2 );
    double eigenvalues[M-1] = { fabs( aver_ncons_params[V_DISP] ), fabs( aver_ncons_params[V_GAS] ),
                                fabs( aver_ncons_params[V_DISP] + c1 ), fabs( aver_ncons_params[V_DISP] - c1 ),
                                fabs( aver_ncons_params[V_GAS] + c2 ), fabs( aver_ncons_params[V_GAS] - c2 ) };
    double max_eigenvalue = 0.0;
    for ( int iComp = 0; iComp < M - 1; iComp++ ) {
        if ( eigenvalues[iComp] > max_eigenvalue )
            max_eigenvalue = eigenvalues[iComp];
    }
    for ( int iComp = 0; iComp < M; iComp++ )
        charact_part[iComp] = 0.5 * max_eigenvalue * ( left_params[iComp] - right_params[iComp] );

    for ( int iComp = 0; iComp < M; iComp++ )
        (*flux)[iComp] = aver_flux[iComp] + charact_part[iComp];
}

// ������ ������ ���������������� ������ � ������ ����� ��� ����� �� ������ ��������� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_params - ������ "��������������" ���������� ����� �� �������
// right_params - ������ "��������������" ���������� ������ �� �������
// rhst - ������� ����� ���������������� ������ � ������ �����
void calc_rhst_term_cir4( const struct ParametersCommon* paramsc, const double left_params[M], const double right_params[M], array1D* rhst ) {

    double aver_params[M]; // ����������� ��������� �� �����

    // ������ ����������� ���������� �� �����
    for ( int iComp = 0; iComp < M; iComp++ )
        aver_params[iComp] = 0.5 * ( left_params[iComp] + right_params[iComp] );

    // ������ ����������������� ������� ������ ������
    rhst_cons( paramsc, aver_params, rhst, 0 );

    // ��������� �� �������� �������� ���� ���������� ����
    for ( int iComp = 0; iComp < M; iComp++ )
        (*rhst)[iComp] *= right_params[B_DISP] - left_params[B_DISP];

}