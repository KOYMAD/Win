// minmod.cc
// Повышение порядка аппроксимации метода за счет кусочно-линейного восполнения сеточных консервативных переменных
// с использованием ограничителя minmod.
// (c) Уткин Павел, 2016
// Создан: 22 августа 2016 г.

#include "minmod.h"

// Реконструкция консервативных переменных с использованием ограничителя minmod
// paramsc - структура с основными параметрами вычислительного эксперимента
// params1d - структура с параметрами одномерной задачи
// u_prev - вектора примитивных переменных во всех ячейках расчетной области
// xc - массив координат центров ячеек сетки
// slopes - вектора наклонов во всех внутренних ячейках расчетной области
void reconstruction( const struct ParametersCommon* paramsc, const struct Parameters1d* params1d, double **u_prev, double *xc, double **slopes, int number_of_scalars ) {

    for ( int i = 1; i < params1d->cells_number - 1; i++ ) { // цикл по всем внутренним ячейкам расчетной области
        // преобразование в консервативную форму
        double u_left[M];
        double u_center[M];
        double u_right[M];
        convert_noncons_to_cons( paramsc, u_prev[i-1], u_left, number_of_scalars );
        convert_noncons_to_cons( paramsc, u_prev[i], u_center, number_of_scalars );
        convert_noncons_to_cons( paramsc, u_prev[i+1], u_right, number_of_scalars );
        for ( int j = 0; j < M; j++ ) { // цикл по компонентам вектора консервативных переменных
            double slope_cand_1 = ( u_center[j] - u_left[j] ) / ( xc[i] - xc[i-1] ); // левый кандидат на значение производной
            double slope_cand_2 = ( u_right[j] - u_center[j] ) / ( xc[i+1] - xc[i] ); // правый кандидат на значение производной
            slopes[i][j] = 0.5 * ( sign( slope_cand_1, paramsc->eps_general ) + sign( slope_cand_2, paramsc->eps_general ) ) * 
                min( fabs( slope_cand_1 ), fabs( slope_cand_2 ) );
        }
    }
    // для простоты в граничных ячейках кладем производные равными нулю - первый порядок
    for ( int j = 0; j < M; j++ ) { // цикл по компонентам вектора консервативных переменных
        slopes[0][j] = 0.0;
        slopes[params1d->cells_number-1][j] = 0.0;
    }

}