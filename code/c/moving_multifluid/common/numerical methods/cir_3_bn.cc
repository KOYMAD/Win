/*
 * cir_3.cc
 *
 * ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������. ������ �����������.
 *
 * ����������� �.�., ��������� �.�., ������� �.�. �������������� ������� ���������� ������� ��������������� ������
 * ���������. - �.: ���������, 2001. - �. 62. - ������� (2.3.22), (2.3.23) + ������������� ������� S �� ����� ��� S ( 0.5 * ( u_left + u_right ) ).
 *
 * (c) ����� �����, 2013
 *
 * ������: 31 ��� 2013 �.
 *
 */

#include "cir_3_bn.h"

#include "utils.h"
#include "math_utils.h"

// ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������, ������ �����������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// left_params[M] - ������ ����������� ���������� � ������ ����� �� �������������� (in)
// center_params[M] - ������ ����������� ���������� � �������������� ������ (in)
// right_params[M] - ������ ����������� ���������� � ������ ������ �� �������������� (in)
// dt - ��������� ��� (in)
// h - ���������������� ��� (in)
// solution[M] - ������ ����������� ���������� � �������������� ������ �� ��������� ���� (out)
void cir_3( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M], double dt, double h, double solution[M] ) {

    double A[M][M]; // ������� �������

    // ������������� �������-������� �� ����� � ������ ������ ������
    double left_edge_sol[M];
    double right_edge_sol[M];

    // "������" ����� ����� ������
    double left_flux[M];
    double right_flux[M];
    
    if( !check_matrix_decomposition( paramsc, center_params ) ) {
        printf( "\ncir_3 -> matrix decomposition is incorrect\n\n" );
        exit( EXIT_FAILURE );
    }

    // ������ ������������� �������-������� �� ������ ������
    calc_edge_solution_cir3( paramsc, left_params, center_params, left_edge_sol );
    calc_edge_solution_cir3( paramsc, center_params, right_params, right_edge_sol );

    calc_A_ncons( paramsc, center_params, A ); // ������ ������� ������� �� ���������� � �������������� ������

    // ����������� "�������" ����� ����� � ������ ����� ������
    calc_flux_cir3( left_edge_sol, A, left_flux );
    calc_flux_cir3( right_edge_sol, A, right_flux );

    for ( int i = 0; i < M; i++ )
        solution[i] = center_params[i] - dt * ( right_flux[i] - left_flux[i] ) / h;

}

// ������ ������������� �������-������� �� ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// left_params[M] - ������ ����������� ���������� � ������ ����� �� ����� (in)
// right_params[M] - ������ ����������� ���������� � ������ ������ �� ����� (in)
// edge_sol[M] - ������������� �������-������� � ����������� ���������� �� ����� (out)
void calc_edge_solution_cir3( struct ParametersCommon *paramsc, double left_params[M], double right_params[M], double edge_sol[M] ) {

    double S[M][M]; // �������, �������� � ������������� �������-������� �� �����
    double edge_params[M]; // ������������ ����������� ���������� �� �����

    for ( int i = 0; i < M; i++ )
        edge_params[i] = 0.5 * ( left_params[i] + right_params[i] );

    cir_util_cir3( paramsc, edge_params, S ); // ������ �������, �������� � ������������� �������-������� �� �����

    // ������ ������������� �������-������� �� �����
    for ( int i = 0; i < M; i++ ) {
        edge_sol[i] = 0.5 * ( left_params[i] + right_params[i] );
        for ( int j = 0; j < M; j++ )
            edge_sol[i] += 0.5 * S[i][j] * ( left_params[j] - right_params[j] );
    }

}

// ������ �������, �������� � ������������� �������-������� �� �����
// params - ��������� � ����������� ��������������� ������������ (in)
// ncons_params[M] - ������ ����������� ���������� (in)
// m[M][M] - ������� ������� (out)
void cir_util_cir3( struct ParametersCommon *paramsc, double ncons_params[M], double m[M][M] ) {

    double omega_ncons[M][M]; // ������� �� ����������� �������� ���������� ������� ��������� ��
    double omega_ncons_inverse[M][M]; // �������� � ������� �� ����������� �������� ���������� ������� ��������� ��
    double sign_lambda[M][M]; // ������������ ������� �� ������� ����������� ����� ���������� ������� ��������� ��
    double m_tmp[M][M]; // ��������� ������������ ������ omega_ncons � sign_lambda

    calc_omega_ncons( paramsc, ncons_params, omega_ncons );
    calc_sign_lambda( paramsc, ncons_params, sign_lambda );
    calc_omega_ncons_inverse( paramsc, ncons_params, omega_ncons_inverse );
    if ( !check_inverse_matrix( omega_ncons, omega_ncons_inverse, paramsc->eps_general ) ) {
        printf( "Error: the inverse matrix is incorrect.\n" );
        exit( EXIT_FAILURE );
    }
    mult_matrixes( omega_ncons, sign_lambda, m_tmp , M);
    mult_matrixes( m_tmp, omega_ncons_inverse, m, M );

}

/*  ������ "������" ����� ����� ������
    edge_sol[M] - ������������� �������-������� � ����������� ���������� �� ����� (in)
    A_ncons[M][M] - ������� ������� � ����������� ���������� � �������������� ������ (in)
    flux[M] - "�����" ����� ����� ������ (out) */
void calc_flux_cir3( double edge_sol[M], double A[M][M], double flux[M] ) {

    int i, j;

    for ( i = 0; i < M; i++ ) {
        flux[i] = 0.0;
        for ( j = 0; j < M; j++ )
            flux[i] += A[i][j] * edge_sol[j];
    }

}

	