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

/* ��������, �������� �� ������� B �������� � ������� A

   A[M][M] - �������� ������� (in)
   B[M][M] - �������� ������� (in)
   eps - �������� ��������� � ����� � �������� (in)
   
   ���������� true, ���� ��������; false - ����� */
bool check_inverse_matrix( double A[M][M], double B[M][M], double eps );

void mult_matrix_vector( const int n, double A[M][M], const double x[M], double b[M] ) ;

/* �������������� ��������� ���� ��������

   n - ������ ������� (in) 
   a[K_GENERAL_CASE] - ������ ������ (in)
   b[K_GENERAL_CASE] - ������ ������ (in)
   eps - �������� ��������� ���� ����� (in)

   ���������� true, ���� ������� ���������; false - � ��������� ������ */
bool compare_vectors( int n, double a[K_GENERAL_CASE], double b[K_GENERAL_CASE], double eps );

/* ������� max, ������������ ���������� �� ���� �����

   a - ������ ����� (in)
   b - ������ ����� (in)
   
   ���������� ���������� �� ���� ����� */
double max( double a, double b );

/* ������� min, ������������ ���������� �� ���� �����

   a - ������ ����� (in)
   b - ������ ����� (in)
   
   ���������� ���������� �� ���� ����� */
double min( double a, double b );

void copy_1D_int_vector( int a[M], int copy_a[M], int n );

void copy_1D_double_vector( const double a[M], double copy_a[M], int n );

// ������� ��������� H(a)
// ���������� 1.0, ���� a > 0, � 0.0 � ��������� �������
double Heviside_function (double a);

#endif // __MATH_UTILS_H