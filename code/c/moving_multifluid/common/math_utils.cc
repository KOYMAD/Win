/*
 * math_utils.cc
 *
 * Математические рабочие функции.
 *
 * (c) Уткин Павел, 2013
 *
 * Создан: 17 мая 2013 г.
 *
 */

#include "math_utils.h"


// Solution of linear equations system A * x = b by means of LU-decomposition
// Press W.H. et al. Numerical Recipes in C. - Cambridge University Press, 1997. - P. 46 - 48.
// A - initial square matrix
// b - right-hand side vector
// n - actual matrix size
// eps - small value instead of zero (is needed by the original algoritm)
// x - solution
void solve_linear_system( const double A[M][M], const double b[M], const int n, const double eps, double x[M] ) {

    int permut[M]; // array with the information about permutations performed during LU-decomposition process
    double A_c[M][M]; // initial matrix copy to make transformations with

    // creation of the matrix copy for future manipulations with
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            A_c[i][j] = A[i][j];

    // LU-decomposition of A_c matrix
    if ( !ludcmp( A_c, n, eps, permut ) ) {
        printf( "\nsolve_linear_system -> matrix in singular\n\n" );
        exit( EXIT_FAILURE );
    }

    // initialization of solution vector with the right-hand side vector
    for ( int i = 0; i < n; i++ )
        x[i] = b[i];

    // solution of the system using the performed LU-decomposition
    lubksb( A_c, permut, n, x );

}

// LU-decomposition of matrix A using Crout algorithm
// Press W.H. et al. Numerical Recipes in C. - Cambridge University Press, 1997. - P. 46 - 48. 
// A - initial matrix, is transformed at the end
// n - actual matrix size
// eps - small value instead of zero (is needed by the original algoritm)
// permut - array with the information about permutations performed during LU-decomposition process
// Returns: - false - if the initial matrix is singular;
//          - true - otherwise
bool ludcmp( double A[M][M], const int n, const double eps, int permut[M] ) {

    double scale_factor[M];

    double big; // for the scaling and row permutation
    double tmp, sum;
    int imax; // for rows permutation

    // scaling multipliers calcualtion for each matrix row
    for ( int i = 0; i < n; i++ ) { // rows cycle
        big = 0.0;
        for ( int j = 0; j < n; j++ ) { // columns cycle
            if ( ( tmp = fabs( A[i][j] ) ) > big )
                big = tmp;
        }
        if ( big == 0.0 ) return false;
        scale_factor[i] = 1.0 / big;
    }

    // Crout algorithm
    for ( int j = 0; j < n; j++ ) { // columns cycle
        // calculation of the elements above the main diagonal
        for ( int i = 0; i < j; i++ ) {
            sum = A[i][j];
            for ( int k = 0; k < i; k++ )
                sum -= A[i][k] * A[k][j];
            A[i][j] = sum;
        }
        big = 0.0;
        // calculation of the rest elements
        for ( int i = j; i < n; i++ ) {
            sum = A[i][j];
            for ( int k = 0; k < j; k++ )
                sum -= A[i][k] * A[k][j];
            A[i][j] = sum;
            if ( ( tmp = scale_factor[i] * fabs( sum ) ) >= big ) {
                big = tmp;
                imax = i;
            }
        }
        // check the necessity of rows permutation
        if ( j != imax ) {
            for ( int k = 0; k < n; k++ ) {
                tmp = A[imax][k];
                A[imax][k] = A[j][k];
                A[j][k] = tmp;
            }
            scale_factor[imax] = scale_factor[j];
        }
        permut[j] = imax;
        if ( A[j][j] == 0.0 ) A[j][j] = eps;
        if ( j != n - 1 ) {
            tmp = 1.0 / A[j][j];
            for ( int i = j + 1; i < n; i++ )
                A[i][j] *= tmp;
        }
    }

    return true;

}

// Solution of linear system A * x = b where the preliminary LU-decomposition of matrix A is performed
// Press W.H. et al. Numerical Recipes in C. - Cambridge University Press, 1997. - P. 46 - 48.
// A - LU-decomposed matrix of the system
// permut - array with the information about permutations performed during LU-decomposition process
// n - actual matrix size
// b - input - right-hand side vector, output - vector-solution
void lubksb( double A[M][M], int permut[M], const int n, double b[M]) {

    int ip;
    int ii = 0;
    double sum;

    for ( int i = 0; i < n; i++ ) {
        ip = permut[i];
        sum = b[ip];
        b[ip] = b[i];
        if ( ii != 0 ) {
            for ( int j = ii - 1; j < i; j++ )
                sum -= A[i][j] * b[j];
        }
        else if ( sum != 0.0 ) {
            ii = i + 1;
        }
        b[i] = sum;
    }
    
    for ( int i = n - 1; i >= 0; i-- ) {
        sum = b[i];
        for ( int j = i + 1; j < n; j++ )
            sum -= A[i][j] * b[j];
        b[i] = sum / A[i][i];
    }

}

// Euclidean vector norm calculation
// a - vector
// n - actual vector size
// Returns: the Euclidean norm of the input vector
double norm( double *a, const int n ) {

    int i;
    double norm = 0.0;

    for ( i = 0; i < n; i++ )
        norm += pow( a[i], 2.0 );

    return sqrt( norm );

}

// Calculation the inverse matrix
// Press W.H. et al. Numerical Recipes in C. - Cambridge University Press, 1997. - P. 48.
// A - initial matrix
// n - actual matrix size
// eps - small value instead of zero (is needed by the original algoritm)
// B - inverse matrix
void calc_inverse_matrix( double A[M][M], const int n, const double eps, double B[M][M] ) {

    int permut[M]; // array with the information about permutations performed during LU-decomposition process
    double A_c[M][M]; // initial matrix copy to make transformations with

    // creation of the matrix copy for future manipulations with
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            A_c[i][j] = A[i][j];

    // LU-decomposition of A_c matrix
    if ( !ludcmp( A_c, n, eps, permut ) ) {
        printf( "\ncalc_inverse_matrix -> matrix in singular\n\n" );
        exit( EXIT_FAILURE );
    }

    // columns of inverse matrix calculation
    for ( int iCol = 0; iCol < n; iCol++ ) {
        double *col;
        get_memory_for_1D_double_array(n, &col);
        col[iCol] = 1.0;
        lubksb( A_c, permut, n, col );
        for ( int iRow = 0; iRow < n; iRow++ ) {
            B[iRow][iCol] = col[iRow];
        }
        free(col);
    }

}

/* Перемножение двух квадратных матриц, C = A * B

   A[M][M] - первая матрица (in)
   B[M][M] - вторая матрица (in)

   C[M][M] - результат (out) */
void mult_matrixes( double A[M][M], double B[M][M], double C[M][M], const int n ) {

    int i, j, k;
	
    for ( i = 0; i < n; i++ ) {
        for ( j = 0; j < n; j++ ) {
            C[i][j] = 0.0;
                for ( k = 0; k < n; k++ )
                    C[i][j] += A[i][k] * B[k][j];
	}
    }

}

/* Функция signum

   x - аргумент (in)
   eps - точность сравнения с нулем (in)
   
   Возвращает sign( x ) */
int sign( double x, double eps ) {

    if ( fabs( x ) < eps )
        return 0;
    else if ( x > 0.0 )
        return 1;
    else
        return -1;

}

/* Проверка, является ли матрица B обратной к матрице A

   A[M][M] - исходная матрица (in)
   B[M][M] - обратная матрица (in)
   eps - точность сравнения с нулем и единицей (in)
   
   Возвращает true, если является; false - иначе */
bool check_inverse_matrix( double A[M][M], double B[M][M], double eps ) {

    int i, j;

    double C[M][M]; /* результат умножения матриц A и B */

    mult_matrixes( A, B, C, M );

    for ( i = 0; i < M; i++ ) {
        for ( j = 0; j < M; j++ ) {
            if ( i != j ) {
                /* вне диагональные элементы */
                if ( fabs( C[i][j] ) > eps )
                    return false;
            }
            else {
                /* диагональные элементы */
                if ( fabs( C[i][j] - 1.0 ) > eps )
                    return false;
            }
        }
    }

    return true;

}

// Multiplication of square matrix on a vector, A * x = b
// n - actual matrix size
// A - initial matrix
// x - vector
// b - result
void mult_matrix_vector( const int n, double A[M][M], const double x[M], double b[M] ) {

    for ( int i = 0; i < n; i++ ) {
        b[i] = 0.0;
        for ( int j = 0; j < n; j++ ) {
            b[i] += A[i][j] * x[j];
	}
    }

}

/* Покомпонентное сравнение двух векторов

   n - размер вектора (in) 
   a[K_GENERAL_CASE] - первый вектор (in)
   b[K_GENERAL_CASE] - второй вектор (in)
   eps - точность сравнения двух чисел (in)

   Возвращает true, если векторы совпадают; false - в противном случае */
bool compare_vectors( int n, double a[K_GENERAL_CASE], double b[K_GENERAL_CASE], double eps ) {

    int i;

    for ( i = 0; i < n; i++ ) {
        if ( fabs( a[i] - b[i] ) > eps )
            return false;
    }

    return true;

}

/* Функция max, возвращающая наибольшее из двух чисел

   a - первое число (in)
   b - второе число (in)
   
   Возвращает наибольшее из двух чисел */
double max( double a, double b ) {

    return ( a < b ) ? b : a;

}

/* Функция min, возвращающая наименьшее из двух чисел

   a - первое число (in)
   b - второе число (in)
   
   Возвращает наименьшее из двух чисел */
double min( double a, double b ) {

    return ( a < b ) ? a : b;

}

void copy_1D_int_vector( int a[M], int copy_a[M], int n ){
    
    for (int i = 0; i < n; i++)
        copy_a[i] = a[i];

} 

void copy_1D_double_vector( const double a[M], double copy_a[M], int n ){
    
    for (int i = 0; i < n; i++)
        copy_a[i] = a[i];

} 

// Функция Хевисайда H(a)
// Возвращает 1.0, если a > 0, и 0.0 в остальных случаях
double Heviside_function (double a){
    
    if (a > 0.0)
        return 1.0;
    else
        return 0.0;

}