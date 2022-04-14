// memory.cc
// Функции для работы с памятью
// (c) Уткин Павел, 2013
// Создан: 12 октября 2013 г.

#include "memory.h"

// Выделение памяти под одномерный массив элементов типа double
// size - размер массива (in) 
// **a - указатель на массив (out)
void get_memory_for_1D_double_array( int size, double **a ) {

    *a = ( double * )malloc( size * sizeof( double ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_1D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }

}

// Выделение памяти под одномерный массив элементов типа int
// size - размер массива (in) 
// **a - указатель на массив (out)
void get_memory_for_1D_int_array( int size, int **a ) {

    *a = ( int * )malloc( size * sizeof( int ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_1D_int_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }

}

// Выделение памяти под двумерный массив элементов типа double
// size1 - первая размерность массива (in)
// size2 - вторая размерность массива (in)
// ***a - указатель на массив (out)
void get_memory_for_2D_double_array( int size1, int size2, double ***a ) {
    
    *a = ( double ** )malloc( size1 * sizeof( double * ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_2D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( int i = 0; i < size1; i++ )
        get_memory_for_1D_double_array( size2, &((*a)[i]) );

}

// Выделение памяти под двумерный массив элементов типа int
// size1 - первая размерность массива (in)
// size2 - вторая размерность массива (in)
// ***a - указатель на массив (out)
void get_memory_for_2D_int_array( int size1, int size2, int ***a ) {
    
    *a = ( int ** )malloc( size1 * sizeof( int * ) );
    if ( NULL == *a ) 
    {
        printf( "get_memory_for_2D_int_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( int i = 0; i < size1; i++ )
        get_memory_for_1D_int_array( size2, &((*a)[i]) );

}

// Выделение памяти под трехмерный массив элементов типа double
// size1 - первая размерность массива (in)
// size2 - вторая размерность массива (in)
// size3 - третья размерность массива (in)
// ****a - указатель на массив (out)
void get_memory_for_3D_double_array( int size1, int size2, int size3, double ****a ) {
    
    *a = ( double *** )malloc( size1 * sizeof( double ** ) );
    if ( NULL == *a ) {
        printf( "get_memory_for_3D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( int j = 0; j < size1; j++ )
        get_memory_for_2D_double_array( size2, size3, &((*a)[j]) );
}
