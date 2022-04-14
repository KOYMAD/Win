/*
 * grid.cc
 *
 * Функции для работы с расчетной сеткой.
 *
 * (c) Уткин Павел, 2013
 *
 * Создан: 17 мая 2013 г.
 *
 */

#include "grid.h"

/* Определение координат узлов и центров ячеек сетки

   left_boundary_x - координата левой границы расчетной области (in)
   right_boundary_x - координата правой границы расчетной области (in)
   cells_num - количество ячеек сетки (in)

   *xc - массив координат центров ячеек (out)
   *x - массив координат узлов сетки (out) */
void build_grid( double left_boundary_x, double right_boundary_x, int cells_num, double *xc, double *x ) {

    int i;
    double h;   /* шаг сетки */

   /* пока реализовано только равномерное распределение узлов сетки */
   h = ( right_boundary_x - left_boundary_x ) / cells_num;
   
   /* координаты узлов */
   for ( i = 0; i < cells_num + 1; i++ ) {
       x[i] = left_boundary_x + i * h;
   }

   /* координаты центров ячеек */
   for ( i = 0; i < cells_num; i++ ) {
       xc[i] = 0.5 * ( x[i] + x[i+1] );
   }

}
