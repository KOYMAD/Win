#ifndef __STRUCTURES_H_
#define __STRUCTURES_H_

#include "constants.h"

struct Primitive_vector // ������ ����������� ����������
{
    double density; // ���������
    double velosity; // ��������
    double pressure; // ��������
};

struct Conservative_vector // ������ �������������� ����������
{
    double specific_mass; // �������� �����
    double specific_momentum; // �������� �������
    double specific_energy; // �������� �������
};

struct Flux_vector // ������ �������
{
    double mass_flux; // ����� �����
    double momentum_flux; // ����� ��������
    double energy_flux; // ����� �������
};

// ���� ��� ������� ��������� ������
struct Block
{
    int number_of_left_cell; // ����� ������ �����
    int number_of_right_cell; // ������ ������ �����
    struct Primitive_vector initial_primitive_vector; // ��������� ������ � �����
};

struct Parameters // ��������� �������
{
    // �����
    int number_of_cells; // ����� �����
    double coordinate_of_left_boundary; // ���������� ����� �������
    double coordinate_of_right_boundary; // ���������� ������ �������

    // ���������� ���������� ����
    double initial_coordinate_of_left_boundary_of_body; // ���������� ����� ������� ����
    double initial_coordinate_of_right_boundary_of_body; // ���������� ������ ������� ����

    // �������� ����
    double body_velosity; // �������� ����
    // ������� ����������� �������, �������� � ����� ����
    double body_cross_section_devided_to_mass; 

    // ������� ������ �� ����� �� �������: 
    int exit_time_cycle; // 1 - �� ���������� time_cycle_iterations ��������, 
                         // 2 - �� ���������� ��������� �������� ������� ������� �������� t_fin
    int time_cycle_iterations; // ����� ������� � ����
    double t_fin; // �������� ������ �������

    // ����� ������� ���� �� �������:
    int time_step_method; // 1 - ������������, �� ������� ������������ �������-���������-����
                          // 2 - �����������, dt
    double cfl; // ����� �������-���������-����
    double dt; // ��� �� �������
    
    // ��������� ����������������
    double length_diml;
    double time_diml;
    double mass_diml;

    // ������ �������� ������
    int number_of_output_files; // ����� ������������ �� ������ ������ ������������� ���������� ���� ��� exit_time_cycle == 2
    int step_of_output_files; // ��� ������ � ���� ������������� ���������� ���� ��� exit_time_cycle == 1

    // ��������� �����
    int flux_calculation_method;
    int max_iter_num;
    double p_max_ratio;
    int approximation_order; // ������� ������������� ���������� ����� �� ������������
    
    // ��������� �������
    struct Block blocks[2]; // ������ ������

    // ����� ���������
    double eps_general;
    double eps_spatial;
    double big;

    // ��������� ��������� ����
    double g; // ���������� ��������

    // ���� ��������� ������� 1 - ������, 2 - ��������� �� �� � ������ ���� M
    // 3 - ������� ������������� �������� �������
    int type_of_left_border;
    int type_of_right_border;
    double Mach_number; // ����� ���� ��

    // ���� ����������� ������
    bool is_debug_regime; // true - ������� ����� ���������� ������ ������
    bool is_�ontinue; // true - �������� ����������� ������ ������� �� � ���������� ������� �������
    bool is_calculate_drag_force; // true - ������� ����� ������� � ������ ���, ����������� �� �������
    bool is_initiate_parameters_in_outer_cells; // true - ������� ������ ����������������
    // ����������� � ������������ � ������� ��������� ������, ����� ������ -big
};

// ���������, ������������ ������� ������ ������������� ������� � ����� ����
struct TimeMoment
{
    int steps_num;  // ������� ���������� �����
    double curr_t;  // ������� �����
};

#endif /* __STRUCTURES_H_ */