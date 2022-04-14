/*
 * grid.h
 *
 * ‘ункции дл€ работы с расчетной сеткой.
 *
 * (c) ”ткин ѕавел, 2013
 *
 * —оздан: 17 ма€ 2013 г.
 *
 */

#ifndef __GRID_H_
#define __GRID_H_


/* ќпределение координат узлов и центров €чеек сетки

   left_boundary_x - координата левой границы расчетной области (in)
   right_boundary_x - координата правой границы расчетной области (in)
   cells_num - количество €чеек сетки (in)

   *xc - массив координат центров €чеек (out)
   *x - массив координат узлов сетки (out) */
void build_grid( double left_boundary_x, double right_boundary_x, int cells_num, double *xc, double *x );

#endif /* __GRID_H_ */