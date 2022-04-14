#ifndef __MEMORY_H
#define __MEMORY_H

#include "structures.h"
#include "stdio.h"
#include "stdlib.h"

void get_memory_for_1D_conservative_vector_array( int size, struct Conservative_vector **a );

void get_memory_for_1D_int_array( int size, int **a );

void get_memory_for_1D_flux_vector_array( int size, struct Flux_vector **a );

#endif /* __MEMORY_H */