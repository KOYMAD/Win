// io.h
// Общие функции ввода/вывода, не зависящие от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 26 августа 2018 г.

#ifndef __IO_H_
#define __IO_H_

#include "struct.h"

// Заполнение структуры с основными параметрами вычислительного эксперимента
// params - дескриптор конфигурационного файла задачи
// paramsc - структура с основными параметрами вычислительного эксперимента
void fill_parameters_common( FILE *params, struct ParametersCommon* paramsc );

// Запись файла с информацией для возможности продолжения расчета
// files_directory - директория, в которой происходит расчет (in)
// time_moment - строка с информацией о моменте времени или числе шагов, для которых
// в последний раз был записан файл с промежуточными результатами (in)
// filename - имя последнего записанного файла (in)
void write_restart_info( char *files_directory, char *time_moment, char *filename );



/* Считывание файла с решением для реализации возможности продолжения прерванного расчета 

   params - структура с параметрами вычислительного эксперимента (in)
   restart_info_file - дескриптор файла с информацией о последнем доступном файле с промежуточными результатами (in)

   **initial_solution - массив векторов в примитивных переменных - начальных условий в каждой ячейке (out) */
void read_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, FILE *restart_info_file, double **initial_solution );

// Попытка считать файл с информацией о рестарте и, в случае успеха, получение дескриптора
// файла с последними записанными распределениями
// paramsc - структура с основными параметрами вычислительного эксперимента (in) 
// files_directory - директория, в которой происходит расчет (in)
// time_mom - структура, определяющая текущий момент времени (out)
// file_to_read - дескриптор файла с распределениями для рестарта (out)
void try_to_restart( struct ParametersCommon *paramsc, char *files_directory, struct TimeMoment *time_mom, FILE *file_to_read );

#endif /* __IO_COMMON_H_ */