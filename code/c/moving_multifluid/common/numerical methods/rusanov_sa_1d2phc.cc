// rusanov_sa_1d2phc.cc
// ����� �������� ���������� �������������� ��������� Saurel-Abgrall, 1D ������
// ����������� ��: Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) ����� �����, 2018
// ������: 4 ����� 2018 �.

#include "rusanov_sa_1d2phc.h"

// ����� �������� ���������� �������������� ��������� Saurel-Abgrall, 1D ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� �����������
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// slopes_left - ������ �������� � ������ ����� �� ��������������
// slopes_center - ������ �������� � �������������� ������
// slopes_right - ������ �������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ �� ��������� ����
// configuration_pressure - ���������������� ��������
void rusanov_1d( struct ParametersCommon* paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
                 double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                 double dt, double h, double solution_ncons[M], int step_number, int n, double curr_time, double *configuration_pressure ) {

    double solution_ncons_Lh[M]; // ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������

    Lh_rusanov_1d( paramsc, left_ncons, center_ncons, right_ncons, slopes_left, slopes_center, slopes_right,
        dt, h, solution_ncons_Lh ); // �������� ���������������� ���������
    if ( paramsc->pressure_relaxation == true )
        Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_pressure  ); // �������� ��������������� ���������
    else {
        for ( int i = 0; i < n; i++ )
            solution_ncons[i] = solution_ncons_Lh[i];
    }

}

// ��������������� �������� ������ ��������, 1D ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// slopes_left - ������ �������� � ������ ����� �� ��������������
// slopes_center - ������ �������� � �������������� ������
// slopes_right - ������ �������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
void Lh_rusanov_1d( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
                    const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                    const double dt, const double h, double solution_ncons[M] ) {

    // ������������� �������� �������
    double left_minus_ncons[M];
    double left_minus_cons[M];
    double left_plus_ncons[M];
    double left_plus_cons[M];
    double right_minus_ncons[M];
    double right_minus_cons[M];
    double right_plus_ncons[M];
    double right_plus_cons[M];
    for ( int j = 0; j < M; j++ ) {
        // ������������� �������� �� ����� ������� �����
        convert_noncons_to_cons( paramsc, left_ncons, left_minus_cons, 0 );
        left_minus_cons[j] += 0.5 * h * slopes_left[j];
        convert_cons_to_noncons( paramsc, left_minus_cons, left_minus_ncons, 0 );
        // ������������� �������� �� ����� ������� ������
        convert_noncons_to_cons( paramsc, center_ncons, left_plus_cons, 0 );
        left_plus_cons[j] -= 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, left_plus_cons, left_plus_ncons, 0 );
        // ������������� �������� �� ������ ������� �����
        convert_noncons_to_cons( paramsc, center_ncons, right_minus_cons, 0 );
        right_minus_cons[j] += 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, right_minus_cons, right_minus_ncons, 0 );
        // ������������� �������� �� ������ ������� ������
        convert_noncons_to_cons( paramsc, right_ncons, right_plus_cons, 0 );
        right_plus_cons[j] -= 0.5 * h * slopes_right[j];
        convert_cons_to_noncons( paramsc, right_plus_cons, right_plus_ncons, 0 );
    }

    // ������ ������ ����� ����� ������� ������
    array1D left_flux( M );
    double splus_left;
    rusanov_flux_1d( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_left );

    // ������ ������ ����� ������ ������� ������
    array1D right_flux( M );
    double splus_right;
    rusanov_flux_1d( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_right );

    // ������������� ���������������� ������������ � ������ ������
    array1D rhst( M ); // ���������������� ������ ������ ������
    rhst_ncons( paramsc, center_ncons, &rhst, 0 );
    double dB1dx = 0.5 * ( right_plus_ncons[B_DISP] - left_minus_ncons[B_DISP] ) / h; // ������������� ��������� �������� ���� ���������� ����
 
    // ������������ ���������� �������-������� �� ��������� ���� �� ������� ��� ������ ���������� (�������� ���� ���������� ����)
    double center_cons[M]; // ������ �������������� ���������� � ������� �������������� ������
    convert_noncons_to_cons( paramsc, center_ncons, center_cons, 0 );
    double solution_cons[M]; // ������-������� �������������� ���������� �� ��������� ���� �� �������
    for ( int i = 1; i < M; i++ )
        solution_cons[i] = center_cons[i] - dt * ( right_flux[i] - left_flux[i] ) / h + dt * dB1dx * rhst[i];

    // ������ �������� �������� ���� ���������� ����
    double u_i = calc_u_i( paramsc, center_ncons );
    double t1 = - 0.5 * u_i * ( right_plus_ncons[B_DISP] - left_minus_ncons[B_DISP] );
    double t2 = 0.5 * ( splus_right * ( right_plus_ncons[B_DISP] - right_minus_ncons[B_DISP] ) - splus_left * ( left_plus_ncons[B_DISP] - left_minus_ncons[B_DISP] ) );
    solution_cons[B_DISP] = center_ncons[B_DISP] + dt * ( t1 + t2 ) / h;
    
    convert_cons_to_noncons( paramsc, solution_cons, solution_ncons, 0 );

}

// ����� �������� ������� ������� � ���������� �����, 1D ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� ����� �� �������
// right_ncons - ������ ����������� ���������� ������ �� �������
// flux - ������������ ������ ������
// smax - ������ ��� �������� S+
void rusanov_flux_1d( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M], array1D* flux, double *splus ) {

    // ������ �������� ����������������� ������ �� ���������� ����� � ������ �� �������
    array1D left_diff_flux( M );
    diff_flux_ncons( paramsc, left_ncons, &left_diff_flux, 0 );
    array1D right_diff_flux( M );
    diff_flux_ncons( paramsc, right_ncons, &right_diff_flux, 0 );

    // ������ ������ ��� �������� ����� S+ � ������ ����� ���
    // ������������ ������������� ���� (������� 2009 �., ���. 329)
    calc_splus_rusanov( paramsc, left_ncons, right_ncons, splus );

    // ������ �������� �������������� ���������� ����� � ������ �� �������
    double left_cons[M];
    convert_noncons_to_cons( paramsc, left_ncons, left_cons, 0 );
    double right_cons[M];
    convert_noncons_to_cons( paramsc, right_ncons, right_cons, 0 );
	
    for ( int i = 1; i < M; i++ )
        (*flux)[i] = 0.5 * ( left_diff_flux[i] + right_diff_flux[i] - (*splus) * ( right_cons[i] - left_cons[i] ) );

}

// ������ ������ ��� �������� ����� S+ � ������ ����� ���
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons[M] - ������ ����������� ���������� ����� �� �������
// right_ncons[M] - ������ ����������� ���������� ������ �� �������
// splus - ������ ��� �������� ����� S+
void calc_splus_rusanov( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M], double* splus ) {

    double c_gas_left, c_disp_left; 
    calc_sound_velocity( paramsc, left_ncons, &c_gas_left, &c_disp_left ); // �������� ����� ����� �� ������� ��� ����� ���
    double c_gas_right, c_disp_right; 
    calc_sound_velocity( paramsc, right_ncons, &c_gas_right, &c_disp_right ); // �������� ����� ������ �� ������� ��� ����� ���

    double t1 = max( max( max( fabs( left_ncons[V_GAS] + c_gas_left ), fabs( right_ncons[V_GAS] + c_gas_right ) ), fabs( left_ncons[V_DISP] + c_disp_left ) ),
        fabs( right_ncons[V_DISP] + c_disp_right ) );
    double t2 = max( max( max( fabs( left_ncons[V_GAS] - c_gas_left ), fabs( right_ncons[V_GAS] - c_gas_right ) ), fabs( left_ncons[V_DISP] - c_disp_left ) ),
        fabs( right_ncons[V_DISP] - c_disp_right ) );

    *splus = max( t1, t2 );

}