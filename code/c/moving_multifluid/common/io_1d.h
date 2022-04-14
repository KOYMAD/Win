// io_1d.h
// Функции ввода/вывода, специфичные для одномерной задачи
// (c) Уткин Павел, 2013 - 2018
// Создан: 17 мая 2012 г.

#ifndef __IO_1D_H_
#define __IO_1D_H_

#include <stdio.h>
#include <string.h>

#include "struct.h"

// константы для заполнения имени файла лидирующими нулями
#define TEN_TO_FIFTH    100000
#define TEN_TO_FOURTH   10000
#define TEN_TO_THIRD    1000
#define TEN_TO_SECOND   100
#define TEN             10

// Заполнение структуры с параметрами одномерной задачи
// params - дескриптор конфигурационного файла задачи
// paramsc - структура с основными параметрами вычислительного эксперимента
// params1d - структура с параметрами одномерной задачи
void fill_parameters_1d( FILE *params, struct ParametersCommon* paramsc, struct Parameters1d* params1d );

// Заполнение структуры с параметрами одномерной задачи в результате считывания входного файла, а также обработки аргументов командной строки.
// Проверка корректности задания параметров.
// paramsc - структура с основными параметрами вычислительного эксперимента
// params1d - структура с параметрами одномерной задачи
// input_file_directory - директория, в которой находится файл с параметрами задачи
// output_file_directory - директория, в которую должны быть записаны файлы с результатами
// file_num - текущий номер файла Parameters1d.dat для считывания
void fill_parameters_struct_1d( struct ParametersCommon* paramsc, struct Parameters1d* params1d, const char* input_file_directory, const char* output_file_directory,
                                const int file_num );

// Функция определяет номера ячеек, в которых расположены датчики
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с параметрами одномерной задачи вычислительного эксперимента (in)
// nodes_coord - массив координат узлов (in)
// sensors_cells - массив номеров ячеек, в которых расположены датчики (out)
void get_sensors_cells( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *nodes_coord, int *sensors_cells );

// Установление соответствия между положениями датчиков и номерами ячеек,
// а также открытие файлов для датчиков на запись
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с параметрами одномерной задачи (in)
// nodes_coord - массив координат узлов (in)
// sensors_cells - массив номеров ячеек, в которых расположены датчики (out)
// sensors_files[MAX_SENSORS_NUM] - массив файловых дескрипторов для записи данных с датчиков (out)
void prepare_sensors( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double *nodes_coord, int *sensors_cells, FILE *sensors_files[MAX_SENSORS_NUM] );

/* Запись файла с информацией для возможности продолжения расчета

   files_directory - директория, в которой происходит расчет (in)
   time_moment - строка с информацией о моменте времени или числе шагов, для которых в последний раз был записан файл с промежуточными результатами (in)
   filename - имя последнего записанного файла (in) */
void write_restart_info( char *files_directory, char *time_moment, char *filename );


// Запись решения в файл
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// params1d - структура с параметрами одномерной задачи (in)
// output_directory - директория, куда будет записан файл с точным решением (in)
// cells_number - число ячеек (in) 
// *xc - массив координат центров ячеек (in)
// **v_ncons - двумерный массив примитивных переменных в центрах всех ячеек сетки (in)
// file_number - номер файла с промежуточными результатами (in)
// output_filename - имя текущего файла с промежуточными результатами (out)
// n - реальный размер векторов (in)
void write_solution_1d2phc( struct ParametersCommon *paramsc, struct Parameters1d* params1d, char *output_directory, int cells_number, double *xc, double **v_ncons,
                     int file_number, char *output_filename, double *beta, double *B, double *C, int n, int *ignition_flag, int *initial_ignition, int *status );

// Запись в файл с результатами заголовочной информации для Tecplot
// params - структура с основными параметрами вычислительного эксперимента (in)
// cells_num - число ячеек (in)
// file_to_write - дескриптор файла, открытого для записи (in)
// file_number - номер файла с промежуточными результатами (in)
void write_solution_tecplot_header( struct ParametersCommon *paramsc, struct Parameters1d *params1d, int cells_num, FILE *file_to_write, int file_number ) ;

// Запись файлов со статистической информацией о работе алгоритма релаксации давлений и алгоритма метода Годунова
// paramcs - структура с основными параметрами вычислительного эксперимента
// paramc1d - структура с параметрами одномерной задачи
// time_mom - информация о текущем моменте времени
// debug_info - отладочная информация
// u_prev - вектор примитивных переменных на прошлом шаге
void write_statistics( const struct ParametersCommon *paramsc, struct Parameters1d *params1d, const struct TimeMoment *time_mom,
                       const struct DebugInfo *debug_info, double **u_prev );

// Запись в файл значений искомых параметров в конкретной ячейке в данный момент времени
// paramcs - структура с основными параметрами вычислительного эксперимента
// paramc1d - структура с параметрами одномерной задачи
// curr_t - текущий момент времени
// file_to_write - дескриптор файла для записи
// file_name - имя открываемого файла
// v_ncons - вектор примитивных переменных на текущем шаге
// number_of_variables - количество переменных в векторе v_ncons
void write_time_dependent_information_in_a_cell (const struct ParametersCommon *paramsc, const struct Parameters1d *params1d, double curr_t, char *output_directory, char *file_name, double *v_ncons, int number_of_variables) ;

// Запись показаний датчиков
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// v_ncons - вектора неконсервативных переменных во всех ячейках расчетной области (in)
// time_mom - структура, определяющая текущий момент времени (in)
// sensors_cells - массив номеров ячеек, в которых расположены датчики (in)
void write_sensors( struct ParametersCommon *paramsc, double **v_ncons, struct TimeMoment *time_mom, FILE *sensors_files[MAX_SENSORS_NUM],
                    int sensors_cells[MAX_SENSORS_NUM], int n, double body_velocity );

#endif /* __IO_1D_H_ */