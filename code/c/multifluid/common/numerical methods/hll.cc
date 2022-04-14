// hll.cc
// ����� HLL �� Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
// and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467
// (c) ����� �����, 2015
// ������: 8 ���� 2015 �.

#include "hll.h"

// ����� HLL ���������� �������������� ��������� Saurel-Abgrall
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// debug_info - ��������� � ���������� ����������� (in)
// left_ncons - ������ ����������� ���������� � ������ ����� �� ��������������
// center_ncons - ������ ����������� ���������� � �������������� ������
// right_ncons - ������ ����������� ���������� � ������ ������ �� ��������������
// slopes_left - ������ �������� � ������ ����� �� ��������������
// slopes_center - ������ �������� � �������������� ������
// slopes_right - ������ �������� � ������ ������ �� ��������������
// dt - ��������� ���
// h - ���������������� ���
// solution_ncons - ������ ����������� ���������� � �������������� ������ �� ��������� ����
// step_number - ����� �������� ���� �� ������������
// n - �������� ������ ��������
// is_pressure_relaxation_after_this_step - true, ���� ����� ������� ���������� ������ �� ������� ����������� ����� ��������� ���������� ��������; false - �����
//                              ���� � ������ ������ ������ �� ����� ��������� ���������� ��������, �� true �� �������� �� ���������, ��� ��� ������ ����� �������������� � parameters.dat
// number_of_scalars - ����� ����������� ��������
// curr_time - ������� ������ �������
// configuration_presure - ���������������� ��������
void hll( struct ParametersCommon* paramsc, struct Parameters1d* params1d, struct DebugInfo *debug_info, double left_ncons[M], double center_ncons[M],
          double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
          double dt, double h, double solution_ncons[M], int step_number, int n, bool is_pressure_relaxation_after_this_step, int number_of_scalars, double curr_time, double *configuration_presure, double S_center, double S_left, double S_right, int l ) {
    double solution_ncons_Lh[M]; // ������ ����������� ���������� � �������������� ������ ����� �������� ���������������� ���������
if (l == 1049)
    printf("\n %lf", center_ncons[P_GAS]);
    Lh_HLL( paramsc, left_ncons, center_ncons, right_ncons, slopes_left, slopes_center, slopes_right,
        dt, h, solution_ncons_Lh, n, number_of_scalars, S_center, S_left, S_right );
    if (l == 1049)
    printf("\n %lf", solution_ncons_Lh[P_GAS]);// �������� ���������������� ���������
    if ( paramsc->pressure_relaxation == true && is_pressure_relaxation_after_this_step == true ){
        Lr( paramsc, params1d, debug_info, left_ncons, solution_ncons_Lh, right_ncons, dt, h, solution_ncons, step_number, n, curr_time, configuration_presure); // �������� ��������������� ���������
        if (l == 1049)
    printf("\n %lf", solution_ncons[P_GAS]);
    }
    else {
        for ( int i = 0; i < n; i++ )
            solution_ncons[i] = solution_ncons_Lh[i];
    }

}

// ��������������� �������� ������ HLL
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
// n - �������� ������ �������� ��� ����� ����������� ��������
// number_of_scalars - ���������� �������������� ���������, ��� �� ���������� ����������� ��������
void Lh_HLL( const struct ParametersCommon* paramsc, const double left_ncons[M], const double center_ncons[M],
             const double right_ncons[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
             const double dt, const double h, double solution_ncons[M], int n, int number_of_scalars, double S_center, double S_left, double S_right  ) {

    // ������������� �������� �������
    double left_minus_ncons[M];
    double left_minus_cons[M];
    double left_plus_ncons[M];
    double left_plus_cons[M];
    double right_minus_ncons[M];
    double right_minus_cons[M];
    double right_plus_ncons[M];
    double right_plus_cons[M];
    for ( int j = 0; j < n; j++ ) {

        // ������������� �������� �� ����� ������� �����
        convert_noncons_to_cons( paramsc, left_ncons, left_minus_cons, number_of_scalars );
        left_minus_cons[j] += 0.5 * h * slopes_left[j];
        convert_cons_to_noncons( paramsc, left_minus_cons, left_minus_ncons, number_of_scalars );
        // ������������� �������� �� ����� ������� ������
        convert_noncons_to_cons( paramsc, center_ncons, left_plus_cons, number_of_scalars );
        left_plus_cons[j] -= 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, left_plus_cons, left_plus_ncons, number_of_scalars );
        // ������������� �������� �� ������ ������� �����
        convert_noncons_to_cons( paramsc, center_ncons, right_minus_cons, number_of_scalars );
        right_minus_cons[j] += 0.5 * h * slopes_center[j];
        convert_cons_to_noncons( paramsc, right_minus_cons, right_minus_ncons, number_of_scalars );
        // ������������� �������� �� ������ ������� ������
        convert_noncons_to_cons( paramsc, right_ncons, right_plus_cons, number_of_scalars );
        right_plus_cons[j] -= 0.5 * h * slopes_right[j];
        convert_cons_to_noncons( paramsc, right_plus_cons, right_plus_ncons, number_of_scalars );
    }

    // ������ ������ ����� ����� ������� ������

    array1D left_flux( M );
    double splus_l, sminus_l; // ������ ������ ��� ��������� ���� S+ � S- ��� ������ ������� �� ����� �����
    hll_flux( paramsc, left_minus_ncons, left_plus_ncons, &left_flux, &splus_l, &sminus_l, n, number_of_scalars);

    // ������ ������ ����� ������ ������� ������
    array1D right_flux( M );
    double splus_r, sminus_r; // ������ ������ ��� ��������� ���� S+ � S- ��� ������ ������� �� ������ �����
    hll_flux( paramsc, right_minus_ncons, right_plus_ncons, &right_flux, &splus_r, &sminus_r, n, number_of_scalars );
    //printf( "\n %lf %lf", left_flux[P_GAS],right_flux[P_GAS]);
    // ������������� ���������������� ������������ � ������ ������
    array1D rhst( M ); // ���������������� ������ ������ ������
    array1D additional( M ); // ������� � ������ ����� ��� ����������� �������
    rhst_ncons(h, splus_l, sminus_l, splus_r, sminus_r, paramsc, center_ncons, &rhst, &additional, number_of_scalars, S_center, S_left, S_right );
    double dB2dx; // ������������� ��������� �������� ���� 
    calc_dB2( left_minus_cons, left_plus_cons, right_minus_cons, right_plus_cons, splus_l, sminus_l, splus_r, sminus_r, &dB2dx, S_center, S_left, S_right);
    // ������������ ���������� �������-������� �� ��������� ���� �� ������� ��� ������ ���������� (�������� ���� ���������� ����)
    double center_cons[M]; // ������ �������������� ���������� � ������� �������������� ������
    convert_noncons_to_cons( paramsc, center_ncons, center_cons, number_of_scalars );
    double solution_cons[M]; // ������-������� �������������� ���������� �� ��������� ���� �� �������
    for ( int i = 1; i < n; i++ ){
        solution_cons[i] = center_cons[i] - dt * ( right_flux[i] - left_flux[i] ) / h + dt *( dB2dx / h * rhst[i] + additional[i] );
    }

 
    // ������ �������� �������� ���� ���������� ����
    double u_i = calc_u_i( paramsc, center_ncons );
    double t1 = ( u_i * ( splus_r * right_minus_cons[B_DISP] - sminus_r * right_plus_cons[B_DISP] ) + splus_r * sminus_r * ( right_plus_cons[B_DISP] - right_minus_cons[B_DISP] ) )
        / ( splus_r - sminus_r );
    double t2 = ( u_i * ( splus_l * left_minus_cons[B_DISP] - sminus_l * left_plus_cons[B_DISP] ) + splus_l * sminus_l * ( left_plus_cons[B_DISP] - left_minus_cons[B_DISP] ) )
        / ( splus_l - sminus_l );        
        solution_cons[B_DISP] = center_ncons[B_DISP] - dt *( t1 - t2 ) / h ;
        if (solution_cons[B_DISP] > 0.94)
            solution_cons[B_DISP] = 0.94;
  
    convert_cons_to_noncons( paramsc, solution_cons, solution_ncons, number_of_scalars);

   
}

// ����� Harten - Lax - van Leer (HLL) ������� ������� � ���������� �����
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons - ������ ����������� ���������� ����� �� �������
// right_ncons - ������ ����������� ���������� ������ �� �������
// flux - ������������ ������ ������
// splus - ������ ��� �������� ����� S+
// sminus - ������ ��� �������� ����� S-
// n - �������� ������ ��������
void hll_flux( const struct ParametersCommon *paramsc, const double left_ncons[M], const double right_ncons[M],
               array1D* flux, double* splus, double* sminus, int n, int number_of_scalars ) {

    // ������ �������� ����������������� ������ �� ���������� ����� � ������ �� �������
    array1D left_diff_flux( M );
    diff_flux_ncons( paramsc, left_ncons, &left_diff_flux, number_of_scalars);
    array1D right_diff_flux( M );
    diff_flux_ncons( paramsc, right_ncons, &right_diff_flux, number_of_scalars );
    //printf( "\n %lf %lf", left_diff_flux[P_GAS],right_diff_flux[P_GAS]);
    // ������ ������ ��� ��������� ���� S+ � S- ��� ������� ����
    double splus_g, sminus_g;
    calc_splus_sminus( paramsc, left_ncons, right_ncons, GAS_PHASE, &splus_g, &sminus_g );

    // ������ ������ ��� ��������� ���� S+ � S- ��� ���������� ����
    double splus_s, sminus_s;
    calc_splus_sminus( paramsc, left_ncons, right_ncons, DISPERSED_PHASE, &splus_s, &sminus_s );

    // ������ ������ ��� ��������� ���� S+ � S- ��� ������ �������
    *splus = max( splus_g, splus_s );
    *sminus = min( sminus_g, sminus_s );

    // ������ �������� �������������� ���������� ����� � ������ �� �������
    double left_cons[M];
    convert_noncons_to_cons( paramsc, left_ncons, left_cons, number_of_scalars);
    double right_cons[M];
    convert_noncons_to_cons( paramsc, right_ncons, right_cons, number_of_scalars );
	
    for ( int i = 1; i < n; i++ )
        (*flux)[i] = ( (*splus) * left_diff_flux[i] - (*sminus) * right_diff_flux[i] +
            (*splus) * (*sminus) * ( right_cons[i] - left_cons[i] ) ) / ( (*splus) - (*sminus) );

}

// ������ ������ ��� ��������� ���� S+ � S-
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// left_ncons[M] - ������ ����������� ���������� ����� �� �������
// right_ncons[M] - ������ ����������� ���������� ������ �� �������
// phase - ����, ��� ������� ����������� ��������
// splus - ������ ��� �������� ����� S+
// sminus - ������ ��� �������� ����� S-
void calc_splus_sminus( const struct ParametersCommon* paramsc, const double left_ncons[M], const double right_ncons[M],
                        const Phase phase, double* splus, double* sminus ) {

    double c_l = calc_sound_velocity_one_phase( paramsc, left_ncons, phase ); // �������� ����� ����� �� �������
    double c_r = calc_sound_velocity_one_phase( paramsc, right_ncons, phase ); // �������� ����� ������ �� �������

    double v_l, v_r;
    if ( phase == GAS_PHASE ) {
        v_l = left_ncons[V_GAS];
        v_r = right_ncons[V_GAS];
    }
    else {
        v_l = left_ncons[V_DISP];
        v_r = right_ncons[V_DISP];
    }

    // � ������:
    // Saurel R., Abgrall R. A multiphase Godunov method for compressible multifluid
    // and multiphase flows // Journal of Computational Physics. - 1999. - V. 150. - P. 425 - 467,
    // �������� ���:
    // *splus = max( 0.0, max( v_l + c_l, v_r + c_r ) );
    // *sminus = min( 0.0, min( v_l - c_l, v_r - c_r ) );
    // �� � ������������ ������, �� ������� ��� ������ ���������:
    // Davis S.F. Simplified second-order Godunov-type methods // SIAM J. Sci. Stat. Comput. - 1988.
    // - V. 9, No. 3. - P. 445 - 473,
    // ��� � � ����� ���� (������� 2009 �., ���. 328), �������� ���:
    *splus = max( v_l + c_l, v_r + c_r );
    *sminus = min( v_l - c_l, v_r - c_r );
    // � ���������� ������ ������� ���� �� ������, ����� - ��������

}

// ������������� ������������� dB2 � ������ �����
// left_minus - ����� �����, ������ �������������� ���������� ����� �� �������
// left_plus - ����� �����, ������ �������������� ���������� ������ �� �������
// right_minus - ������ �����, ������ �������������� ���������� ����� �� �������
// right_plus - ������ �����, ������ �������������� ���������� ������ �� �������
// splus_l - ������ ��� �������� ����� S+ �� ����� �����
// sminus_l - ������ ��� �������� ����� S- �� ����� �����
// splus_r - ������ ��� �������� ����� S+ �� ������ �����
// sminus_l - ������ ��� �������� ����� S- �� ������ �����
// dB2dx - ������� �������������

// ������ ����������� �� S -- done
void calc_dB2( const double left_minus[M], const double left_plus[M], const double right_minus[M], const double right_plus[M],
               double splus_l, double sminus_l, double splus_r, double sminus_r, double* dB2dx, double S_center, double S_left, double S_right) {

    *dB2dx = ( ( splus_r * right_minus[B_DISP] - sminus_r * right_plus[B_DISP] ) / ( splus_r - sminus_r )
        - ( splus_l * left_minus[B_DISP] - sminus_l * left_plus[B_DISP] ) / ( splus_l - sminus_l ) );



}

