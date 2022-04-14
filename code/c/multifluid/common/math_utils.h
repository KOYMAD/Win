// math_utils.h
// Mathematical util functions.
// (c) Pavel Utkin, 2013
// Created 24.05.2013

#ifndef __MATH_UTILS_H
#define __MATH_UTILS_H

#include "struct.h"
#include "memory.h"

void solve_linear_system( const double A[M][M], const double b[M], const int n, const double eps, double x[M] ) ;

bool ludcmp( double A[M][M], const int n, const double eps, int permut[M] );

void lubksb( double A[M][M], int permut[M], const int n, double b[M] );

double norm(  double *a, const int n );

void calc_inverse_matrix( double A[M][M], const int n, const double eps, double B[M][M] ) ;

// Two square matrixes multiplication, C = A * B
// A - the first matrix
// B - the second matrix
// C - the result
void mult_matrixes( double A[M][M], double B[M][M], double C[M][M], const int n ) ;

// signum function
// x - arguement
// eps - accuracy of the comparison with zero
// Returns: sign( x )
int sign( double x, double eps );

/* ѕроверка, €вл€етс€ ли матрица B обратной к матрице A

   A[M][M] - исходна€ матрица (in)
   B[M][M] - обратна€ матрица (in)
   eps - точность сравнени€ с нулем и единицей (in)
   
   ¬озвращает true, если €вл€етс€; false - иначе */
bool check_inverse_matrix( double A[M][M], double B[M][M], double eps );

void mult_matrix_vector( const int n, double A[M][M], const double x[M], double b[M] ) ;

/* ѕокомпонентное сравнение двух векторов

   n - размер вектора (in) 
   a[K_GENERAL_CASE] - первый вектор (in)
   b[K_GENERAL_CASE] - второй вектор (in)
   eps - точность сравнени€ двух чисел (in)

   ¬озвращает true, если векторы совпадают; false - в противном случае */
bool compare_vectors( int n, double a[K_GENERAL_CASE], double b[K_GENERAL_CASE], double eps );

/* ‘ункци€ max, возвращающа€ наибольшее из двух чисел

   a - первое число (in)
   b - второе число (in)
   
   ¬озвращает наибольшее из двух чисел */
double max( double a, double b );

/* ‘ункци€ min, возвращающа€ наименьшее из двух чисел

   a - первое число (in)
   b - второе число (in)
   
   ¬озвращает наименьшее из двух чисел */
double min( double a, double b );

void copy_1D_int_vector( int a[M], int copy_a[M], int n );

void copy_1D_double_vector( const double a[M], double copy_a[M], int n );

// ‘ункци€ ’евисайда H(a)
// ¬озвращает 1.0, если a > 0, и 0.0 в остальных случа€х
double Heviside_function (double a);

#endif // __MATH_UTILS_H