// cir_1.cc
// ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������. ������ �����������.
// ����������� �.�., ��������� �.�., ������� �.�. �������������� ������� ���������� ������� ��������������� ������
// ���������. - �.: ���������, 2001. - �. 62. - ������� (2.3.18), (2.3.19).
// (c) ����� �����, 2013
// ������: 24 ��� 2013 �.

#include "cir_1_bn.h"

#include "utils.h"
#include "math_utils.h"

// ����� �������-��������-��� �������������� ���������� ������� ��������� ����� - ��������, ������ �����������
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// left_params[M] - ������ ����������� ���������� � ������ ����� �� �������������� (in)
// center_params[M] - ������ ����������� ���������� � �������������� ������ (in)
// right_params[M] - ������ ����������� ���������� � ������ ������ �� �������������� (in)
// dt - ��������� ��� (in)
// h - ���������������� ��� (in)
// solution[M] - ������ ����������� ���������� � �������������� ������ �� ��������� ���� (out)
void cir_1( struct ParametersCommon *paramsc, double left_params[M], double center_params[M], double right_params[M],
            double dt, double h, double solution[M] ) {

    double A[M][M]; // ������� �������
    double S[M][M]; // ������� ��� ������������� �������-������� �� ����� ������

    // ������������� �������-������� �� ����� � ������ ������ ������
    double left_edge_sol[M];
    double right_edge_sol[M];

    // "������" ����� ����� ������
    double left_flux[M];
    double right_flux[M];
    
    if( !check_matrix_decomposition( paramsc, center_params ) ) {
        printf( "Error: the matrix decomposition is incorrect.\n" );
        exit( EXIT_FAILURE );
    }

    // ������ ������� ��� ������������� �������-������� �� ����� ������ �� ���������� � �������������� ������
    cir_util_cir1( paramsc, center_params, S );

    // ������ ������������� �������-������� �� ������ ������
    calc_edge_solution_cir1( left_params, center_params, S, left_edge_sol );
    calc_edge_solution_cir1( center_params, right_params, S, right_edge_sol );

    // ������ ������� ������� �� ���������� � �������������� ������
    calc_A_ncons( paramsc, center_params, A );

    // ����������� "�������" ����� ����� � ������ ����� ������
    calc_flux_cir1( left_edge_sol, A, left_flux );
    calc_flux_cir1( right_edge_sol, A, right_flux );

    for ( int i = 0; i < M; i++ )
        solution[i] = center_params[i] - dt * ( right_flux[i] - left_flux[i] ) / h;

}

// ������ �������, �������� � ������������� �������-������� �� �����
// paramsc - ��������� � ������ ����������� ��������������� ������������ (in)
// ncons_params[M] - ������ ����������� ���������� (in)
// S[M][M] - ������� ������� (out)
void cir_util_cir1( struct ParametersCommon *paramsc, double ncons_params[M], double S[M][M] ) {

    /* ������� �� ����������� �������� ���������� ������� ��������� �����-�������� */
    double omega_ncons[M][M];

    /* �������� � ������� �� ����������� �������� ���������� ������� ��������� �����-�������� */
    double omega_ncons_inverse[M][M];
	
    /* ������������ ������� �� ������� ����������� ����� ���������� ������� ��������� �����-�������� */
    double sign_lambda[M][M];

    /* ��������� ������������ ������ omega_ncons � sign_lambda */
    double m_tmp[M][M];

    calc_omega_ncons( paramsc, ncons_params, omega_ncons );
    calc_sign_lambda( paramsc, ncons_params, sign_lambda );
    calc_omega_ncons_inverse( paramsc, ncons_params, omega_ncons_inverse );
    if ( !check_inverse_matrix( omega_ncons, omega_ncons_inverse, paramsc->eps_general ) ) {
        printf( "Error: the inverse matrix is incorrect.\n" );
        exit( EXIT_FAILURE );
    }
    mult_matrixes( omega_ncons, sign_lambda, m_tmp, M );
    mult_matrixes( m_tmp, omega_ncons_inverse, S, M );

}

/*  ������ ������������� �������-������� �� ����� ������
    left_params[M] - ������ ����������� ���������� � ������ ����� �� ����� (in)
    right_params[M] - ������ ����������� ���������� � ������ ������ �� ����� (in)
    S[M][M] - �������, ������������ �� ���������� � �������������� ������ (in)
    edge_sol[M] - ������������� �������-������� � ����������� ���������� �� ����� (out) */
void calc_edge_solution_cir1( double left_params[M], double right_params[M], double S[M][M], double edge_sol[M] ) {

    int i, j;

    /* ������ ������������� �������-������� �� ����� */
    for ( i = 0; i < M; i++ ) {
        edge_sol[i] = 0.5 * ( left_params[i] + right_params[i] );
        for ( j = 0; j < M; j++ )
            edge_sol[i] += 0.5 * ( S[i][j] * ( left_params[j] - right_params[j] ) );
    }

}

/*  ������ "������" ����� ����� ������
    edge_sol[M] - ������������� �������-������� � ����������� ���������� �� ����� (in)
    A_ncons[M][M] - ������� ������� � ����������� ���������� � �������������� ������ (in)
    flux[M] - "�����" ����� ����� ������ (out) */
void calc_flux_cir1( double edge_sol[M], double A[M][M], double flux[M] ) {

    int i, j;

    for ( i = 0; i < M; i++ ) {
        flux[i] = 0.0;
        for ( j = 0; j < M; j++ )
            flux[i] += A[i][j] * edge_sol[j];
    }

}

	