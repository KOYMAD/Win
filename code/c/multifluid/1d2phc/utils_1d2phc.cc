// utils.cc
// ������� �������, ����������� ��� ���������������� ������ ������ ��������� ���� Baer-Nunziato � Saurel-Abgrall
// (c) ����� �����, 2013 - 2015
// ������: 16 ������� 2013 �.

#include "utils_1d2phc.h"



// ������ ���� �������������� �� �������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
// *x - ���������� ����� �����
// **v_ncons - ������� ����������� ���������� � ������� �����
// ���������� ��� �������������� �� �������
double get_time_step( struct DebugInfo *debug_info, const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, const double* x, double **v_ncons, int print ) {

    double new_step = INFINITY1; // �������������� ��� �� �������
   
    if ( paramsc->constant_time_step ) {
        // ������ � ���������� ����� �� �������
        return paramsc->dt;
    }
    else {
        // ������ � ������������ ������� ���� �� �������
        for ( int i = 0; i < params1d->cells_number; i++ ) {
            
            double h = x[i+1] - x[i]; // ������� ������ ������
       //     printf("i = %d\n", i);
            double c1, c2; // �������� ����� ������� � ���������� ����
            calc_sound_velocity(debug_info,  paramsc, v_ncons[i], &c1, &c2 );
            
            double curr_gas_step, curr_solid_step; // ���� �� ������� �� �������� ������������ ��� ������� � ���������� ���
            // ������������� ���� �� ������� �� ������ ���������� ������� ����
            curr_gas_step = h / ( fabs( v_ncons[i][V_GAS] ) + c2 );
            if ( curr_gas_step < new_step )
                new_step = curr_gas_step;
            // ������������� ���� �� ������� �� ������ ���������� ���������� ����
            curr_solid_step = h / ( fabs( v_ncons[i][V_DISP] ) + c1 );
            if ( curr_solid_step < new_step )
                new_step = curr_solid_step;


        }
        return paramsc->cfl * new_step;
    }

}

// �������� ���������� ����� � �������
// params1d - ��������� � ����������� ���������� ������
// *x - ���������� ����� �����
// **v_ncons - ������� ����������� ���������� � ������� �����
// initial_total_mass - �������� ��������� ����� �������� � �������
// mass_diff - ������������� ����������� ��������� ����� �������� � ������� ������������ ��������
void check_mass( const struct Parameters1d* params1d, const double* x, double** v_ncons, const double initial_total_mass, double* mass_diff ) {

    double current_total_mass = get_total_mass( params1d, x, v_ncons ); // ������� ��������� ����� �������� � �������
        
    // ����������� ��������� ����� � �������
    if ( initial_total_mass > 0.0 )
        *mass_diff = fabs( initial_total_mass - current_total_mass ) / initial_total_mass;
    else {
        printf( "check_mass -> initial substance mass is zero\n" );
        exit( EXIT_FAILURE );
    }
    
}

// ������ ������ ����� �������� � �������
// params1d - ��������� � ����������� ���������� ������
// *x - ���������� ����� �����
// **v_cons - ������� "��������������" ���������� � ������� �����
// ���������� ��������� ����� ���� � ���������� ���� � �������
double get_total_mass( const struct Parameters1d* params1d, const double* x, double** v_cons ) {

    double total_mass = 0.0;

    for ( int i = 0; i < params1d->cells_number; i++ ) {
        total_mass += v_cons[i][R_DISP] * ( x[i+1] - x[i] );
        total_mass += v_cons[i][R_GAS] * ( x[i+1] - x[i] );
    }

    return total_mass;

}


void calculate_DCJ(struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *Dcj){
	
	double density_mixture_background_block2 = params1d->block_values[1][R_DISP] * params1d->block_values[1][B_DISP] + params1d->block_values[1][R_GAS] * (1.0 - params1d->block_values[1][B_DISP]);
	double pressure_mixture_background_block2 = params1d->block_values[1][P_DISP] * params1d->block_values[1][B_DISP] + params1d->block_values[1][P_GAS] * (1.0 - params1d->block_values[1][B_DISP]);
	double internal_energy_gas = e_gas(paramsc, params1d->block_values[1][P_GAS], params1d->block_values[1][R_GAS]);
	double internal_energy_disp = e_disp(paramsc, params1d->block_values[1][P_DISP], params1d->block_values[1][R_DISP]);
	double internal_energy_mixture_background_block2 = 1.0 / density_mixture_background_block2 * ( params1d->block_values[1][R_DISP] * params1d->block_values[1][B_DISP] * internal_energy_disp + params1d->block_values[1][R_GAS] * (1.0 - params1d->block_values[1][B_DISP]) * internal_energy_gas);

	double vm0 = 1.0 / density_mixture_background_block2;

	double mu = (paramsc->g2 - 1.0) / (1.0 + paramsc->g2);

	double a4 = (1.0 + mu) * pressure_mixture_background_block2;
	double a3 = -4.0 * mu * ( internal_energy_mixture_background_block2 + pressure_mixture_background_block2 * vm0);
	double a2 = mu * (1.0 + mu) * vm0 * (2.0 * internal_energy_mixture_background_block2 + pressure_mixture_background_block2 * vm0);

	double Discriminant = a3 * a3 - 4.0 * a4 * a2;

	double vcj_1 = (-a3 + sqrt(Discriminant)) / 2.0 / a4;
	double vcj_2 = (-a3 - sqrt(Discriminant)) / 2.0 / a4;

	*Dcj = sqrt( mu * vcj_2 * vm0 * vm0 * ( (1 - mu) * pressure_mixture_background_block2 * vm0 + 2.0 * internal_energy_mixture_background_block2 ) / (vcj_2 - mu * vm0) );
	
}



