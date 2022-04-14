// io_1d.h
// ������� �����/������, ����������� ��� ���������� ������
// (c) ����� �����, 2013 - 2018
// ������: 17 ��� 2012 �.

#ifndef __IO_1D_H_
#define __IO_1D_H_

#include <stdio.h>
#include <string.h>

#include "struct.h"

// ��������� ��� ���������� ����� ����� ����������� ������
#define TEN_TO_FIFTH    100000
#define TEN_TO_FOURTH   10000
#define TEN_TO_THIRD    1000
#define TEN_TO_SECOND   100
#define TEN             10

// ���������� ��������� � ����������� ���������� ������
// params - ���������� ����������������� ����� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
void fill_parameters_1d( FILE *params, struct ParametersCommon* paramsc, struct Parameters1d* params1d );

// ���������� ��������� � ����������� ���������� ������ � ���������� ���������� �������� �����, � ����� ��������� ���������� ��������� ������.
// �������� ������������ ������� ����������.
// paramsc - ��������� � ��������� ����������� ��������������� ������������
// params1d - ��������� � ����������� ���������� ������
// input_file_directory - ����������, � ������� ��������� ���� � ����������� ������
// output_file_directory - ����������, � ������� ������ ���� �������� ����� � ������������
// file_num - ������� ����� ����� Parameters1d.dat ��� ����������
void fill_parameters_struct_1d( struct ParametersCommon* paramsc, struct Parameters1d* params1d, const char* input_file_directory, const char* output_file_directory,
                                const int file_num );

// ������� ���������� ������ �����, � ������� ����������� �������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ ��������������� ������������ (in)
// nodes_coord - ������ ��������� ����� (in)
// sensors_cells - ������ ������� �����, � ������� ����������� ������� (out)
void get_sensors_cells( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *nodes_coord, int *sensors_cells );

// ������������ ������������ ����� ����������� �������� � �������� �����,
// � ����� �������� ������ ��� �������� �� ������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ (in)
// nodes_coord - ������ ��������� ����� (in)
// sensors_cells - ������ ������� �����, � ������� ����������� ������� (out)
// sensors_files[MAX_SENSORS_NUM] - ������ �������� ������������ ��� ������ ������ � �������� (out)
void prepare_sensors( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *nodes_coord, int *sensors_cells, FILE *sensors_files[MAX_SENSORS_NUM] );

/* ������ ����� � ����������� ��� ����������� ����������� �������

   files_directory - ����������, � ������� ���������� ������ (in)
   time_moment - ������ � ����������� � ������� ������� ��� ����� �����, ��� ������� � ��������� ��� ��� ������� ���� � �������������� ������������ (in)
   filename - ��� ���������� ����������� ����� (in) */
void write_restart_info( char *files_directory, char *time_moment, char *filename );


// ������ ������� � ����
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// params1d - ��������� � ����������� ���������� ������ (in)
// output_directory - ����������, ���� ����� ������� ���� � ������ �������� (in)
// cells_number - ����� ����� (in) 
// *xc - ������ ��������� ������� ����� (in)
// **v_ncons - ��������� ������ ����������� ���������� � ������� ���� ����� ����� (in)
// file_number - ����� ����� � �������������� ������������ (in)
// output_filename - ��� �������� ����� � �������������� ������������ (out)
// n - �������� ������ �������� (in)
void write_solution_1d2phc( struct ParametersCommon *paramsc, struct Parameters1d* params1d, char *output_directory, int cells_number, double *xc, double **v_ncons,
                     int file_number, char *output_filename, double *beta, double *B, double *C, int n, int *ignition_flag, int *initial_ignition, int *status );

// ������ � ���� � ������������ ������������ ���������� ��� Tecplot
// params - ��������� � ��������� ����������� ��������������� ������������ (in)
// cells_num - ����� ����� (in)
// file_to_write - ���������� �����, ��������� ��� ������ (in)
// file_number - ����� ����� � �������������� ������������ (in)
void write_solution_tecplot_header( struct ParametersCommon *paramsc, struct Parameters1d *params1d, int cells_num, FILE *file_to_write, int file_number ) ;

// ������ ������ �� �������������� ����������� � ������ ��������� ���������� �������� � ��������� ������ ��������
// paramcs - ��������� � ��������� ����������� ��������������� ������������
// paramc1d - ��������� � ����������� ���������� ������
// time_mom - ���������� � ������� ������� �������
// debug_info - ���������� ����������
// u_prev - ������ ����������� ���������� �� ������� ����
void write_statistics( const struct ParametersCommon *paramsc, struct Parameters1d *params1d, const struct TimeMoment *time_mom,
                       const struct DebugInfo *debug_info, double **u_prev );

// ������ � ���� �������� ������� ���������� � ���������� ������ � ������ ������ �������
// paramcs - ��������� � ��������� ����������� ��������������� ������������
// paramc1d - ��������� � ����������� ���������� ������
// curr_t - ������� ������ �������
// file_to_write - ���������� ����� ��� ������
// file_name - ��� ������������ �����
// v_ncons - ������ ����������� ���������� �� ������� ����
// number_of_variables - ���������� ���������� � ������� v_ncons
void write_time_dependent_information_in_a_cell (const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, double curr_t, char *output_directory, char *file_name, double *v_ncons, int number_of_variables) ;

// ������ ��������� ��������
// paramsc - ��������� � ��������� ����������� ��������������� ������������ (in)
// v_ncons - ������� ���������������� ���������� �� ���� ������� ��������� ������� (in)
// time_mom - ���������, ������������ ������� ������ ������� (in)
// sensors_cells - ������ ������� �����, � ������� ����������� ������� (in)
void write_sensors( struct ParametersCommon *paramsc, double **v_ncons, struct TimeMoment *time_mom, FILE *sensors_files[MAX_SENSORS_NUM],
                    int sensors_cells[MAX_SENSORS_NUM], int n, double body_velocity );

#endif /* __IO_1D_H_ */