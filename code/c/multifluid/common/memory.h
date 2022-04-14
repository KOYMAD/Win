// memory.h
// Функции для работы с памятью
// (c) Уткин Павел, 2013
// Создан: 12 октября 2013 г.

#ifndef __MEMORY_H_
#define __MEMORY_H_

#include <stdio.h>
#include <stdlib.h>

// Выделение памяти под одномерный массив элементов типа double
// size - размер массива (in) 
// **a - указатель на массив (out)
void get_memory_for_1D_double_array( int size, double **a );

void get_memory_for_1D_int_array( int size, int **a ) ;

// Выделение памяти под двумерный массив элементов типа double
// size1 - первая размерность массива (in)
// size2 - вторая размерность массива (in)
// ***a - указатель на массив (out)
void get_memory_for_2D_double_array( int size1, int size2, double ***a );

// Выделение памяти под двумерный массив элементов типа int
// size1 - первая размерность массива (in)
// size2 - вторая размерность массива (in)
// ***a - указатель на массив (out)
void get_memory_for_2D_int_array( int size1, int size2, int ***a );

// Выделение памяти под трехмерный массив элементов типа double
// size1 - первая размерность массива (in)
// size2 - вторая размерность массива (in)
// size3 - третья размерность массива (in)
// ****a - указатель на массив (out)
void get_memory_for_3D_double_array( int size1, int size2, int size3, double ****a );

#endif /* __MEMORY_H_ */