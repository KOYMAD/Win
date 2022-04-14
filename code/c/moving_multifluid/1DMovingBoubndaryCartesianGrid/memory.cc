#include "memory.h"

// ��������� ������ ��� ���������� ������ ��������� ���� Conservative_vector
// size - ������ ������� (in) 
// **a - ��������� �� ������ (out)
void get_memory_for_1D_conservative_vector_array( int size, struct Conservative_vector **a )
{
    *a = ( Conservative_vector * )malloc( size * sizeof( Conservative_vector ) );
    if ( NULL == *a )
    {
        printf( "get_memory_for_1D_Conservative_vector_array -> Can't allocate memory.\n" );
        system("Pause");
    }
}

// ��������� ������ ��� ���������� ������ ��������� ���� int
// size - ������ ������� (in) 
// **a - ��������� �� ������ (out)
void get_memory_for_1D_int_array( int size, int **a )
{
    *a = ( int * )malloc( size * sizeof( int ) );
    if ( NULL == *a )
    {
        printf( "get_memory_for_1D_int_array -> Can't allocate memory.\n" );
        system("Pause");
    }
}

// ��������� ������ ��� ���������� ������ ��������� ���� Conservative_vector
// size - ������ ������� (in) 
// **a - ��������� �� ������ (out)
void get_memory_for_1D_flux_vector_array( int size, struct Flux_vector **a )
{
    *a = ( Flux_vector * )malloc( size * sizeof( Flux_vector ) );
    if ( NULL == *a )
    {
        printf( "get_memory_for_1D_Flux_vector_array -> Can't allocate memory.\n" );
        system("Pause");
    }
}