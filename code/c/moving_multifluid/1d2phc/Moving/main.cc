#include "structures.h"
#include "constants.h"
#include "read_info.h"
#include "memory.h"
#include "initial_and_boundary_data.h"
#include "time_integration.h"
#include "body_dynamic.h"
#include "write_info.h"

int main()
{
    struct Conservative_vector *conservative_variables; // �������������� ���������� � ��������� �������
    int *status; // ������ �������� �����
    struct Parameters parameters; // ��������� �������
    read_parameters ( &parameters ); // ���������� ���������� �������
    dimless_parameters ( &parameters ); // ���������������� ���������� �������
    get_memory_for_1D_conservative_vector_array ( parameters.number_of_cells, &conservative_variables );
    get_memory_for_1D_int_array ( parameters.number_of_cells, &status );

    // ���������� ����� � ������ ������ ����
    double coordinate_of_left_boundary_of_body = parameters.initial_coordinate_of_left_boundary_of_body;
    double coordinate_of_right_boundary_of_body = parameters.initial_coordinate_of_right_boundary_of_body;
    initiate_status ( parameters, coordinate_of_left_boundary_of_body,
        coordinate_of_right_boundary_of_body, status ); // ����������� �������� �����
    struct TimeMoment time_mom = { 0, 0.0 }; // ���������, ������������ ������� ������������ ������ �������
    int zones_count = 0; // ������� ����� ������ � ������������� ���������� ����
    // ������� ��������� ������ � ��������� �������
    initiate_data( &parameters, conservative_variables, status, &time_mom, &zones_count ); // �������������
    printf( "\n> The solution is initiated successfully.\n" );
    bool condition_to_enter_the_cycle;
    is_time_evolution( &parameters, time_mom, &condition_to_enter_the_cycle ); // ������� ����� � ���� �� �������
    double body_velocity = parameters.body_velosity; // �������� ����
    // �������� ����� �� ������� �� ������� � ���������� ����� ��������������
    double time_step_on_current_step, time_step_on_previous_step = 0.0;
    while ( condition_to_enter_the_cycle ) // �������� ���� �� �������
    {
        // ������ ���������� ������
        write_piston_trajectory ( &parameters, coordinate_of_left_boundary_of_body, time_mom );
        time_step_on_current_step = parameters.big;
        calc_body_velosity ( &parameters, &body_velocity, time_step_on_previous_step,
            status, conservative_variables, time_mom ); // ����������� �������� ����
         // ������ �������� �������� ���� �� �������
        calc_time_step ( &parameters, status, conservative_variables, body_velocity, &time_step_on_current_step );
        // ������ �������������� ����������, �������� ���� � �������� �������� �� ����� ���� �� �������
	time_integration ( &parameters, conservative_variables, time_step_on_current_step, status, body_velocity, 
            &coordinate_of_left_boundary_of_body, &coordinate_of_right_boundary_of_body );
        time_mom.curr_t += time_step_on_current_step; // ���������� ������� �������
        time_mom.steps_num++; // ���������� ������ ����
        time_step_on_previous_step = time_step_on_current_step; // ������������� �������� ����������� ���� �� ������� 
        is_time_evolution( &parameters, time_mom, &condition_to_enter_the_cycle ); // ��������� �� ���������� ������
        printf( "%d, %e\n", time_mom.steps_num, time_mom.curr_t * parameters.time_diml );
        if ( ( parameters.exit_time_cycle == FINAL_TIME && time_mom.curr_t >= parameters.t_fin *
            zones_count / parameters.number_of_output_files ) || ( parameters.exit_time_cycle ==
            ITERATIONS && (time_mom.steps_num%parameters.step_of_output_files) == 0 ) )
        {
            write_results( &parameters, status, conservative_variables, time_mom );
            zones_count++;
        }
    }
    printf( "\n> Calculation is succesfully over.\n" );
    return 0;
}