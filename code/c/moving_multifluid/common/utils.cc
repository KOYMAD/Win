// utils.cc
// Рабочие функции, специфичные для рассматриваемого класса систем уравнений типа Baer-Nunziato и Saurel-Abgrall
// вне зависимости от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 26 августа 2018 г.

#include "utils.h"
void initiate_status( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double left, double right, int *status, int *left_ghost, int *right_ghost )
{
    double grid_step = ( params1d->right_boundary_x - params1d->left_boundary_x ) /
        params1d->cells_number; // шаг сетки
    for ( int i = 0 ; i < params1d->cells_number ; i++ )
    {
        if ( left >= i * grid_step ) // если левая граница тела правее левой границы ячейки
        {
            if ( left >= ( i + 1 ) * grid_step ) // если левая граница тела правее правой границы ячейки
            {
                if ( right >= i * grid_step ) // если правая граница тела правее левой границы ячейки
                {
                    if ( right >= ( i + 1 ) * grid_step ) // если правая граница тела правее правой границы ячейки
                        status[i] = INNER;
                    else // если правая граница тела левее правой границы ячейки
                    {
                        printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                        system ( "Pause" );
                    }
                }
                else // если правая граница тела левее левой границы ячейки
                {
                    printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                    system ( "Pause" );
                }
            }
            else // если левая граница тела левее правой границы ячейки
            {
                if ( right >= i * grid_step ) // если правая граница тела правее левой границы ячейки
                {
                    if ( right >= ( i + 1 ) * grid_step ){ // если правая граница тела правее правой границы ячейки
                        status[i] = GHOST;
                        *left_ghost = i;
                    }
                    else // если правая граница тела левее правой границы ячейки
                    {
                        printf( "\ninitiate_status -> all body is inside current cell\n" );
                        system ( "Pause" );
                    }
                }
                else // если правая граница тела левее левой границы ячейки
                {
                    printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                    system ( "Pause" );
                }
            }
        }
        else // если левая граница тела левее левой границы ячейки
        {
            if ( left >= ( i + 1 ) * grid_step ) // если левая граница тела правее правой границы ячейки
            {
                printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                system ( "Pause" );
            }
            else // если левая граница тела левее правой границы ячейки
            {
                if ( right >= i * grid_step ) // если правая граница тела правее левой границы ячейки
                {
                    if ( right >= ( i + 1 ) * grid_step ) // если правая граница тела правее правой границы ячейки
                        status[i] = OUTER;
                    else{ // если правая граница тела левее правой границы ячейки
                        status[i] = GHOST;
                        *right_ghost = i;
                    }

                }
                else // если правая граница тела левее левой границы ячейки
                {
                    if ( right >= ( i + 1 ) * grid_step ) // если правая граница тела правее правой границы ячейки
                    {
                        printf( "\ninitiate_status -> impossible coordinates of body boubdaries\n" );
                        system ( "Pause" );
                    }
                    else // если правая граница тела левее правой границы ячейки
                        status[i] = INNER;
                }
            }
        }
    }
    for ( int i = 0 ; i < params1d->cells_number ; i++ )
    {
        if ( GHOST == status[i] )
        {
            if ( i != 0 && INNER == status[i - 1] && OUTER == status[i + 1] )
                status[i - 1] = BOUNDARY;
            else if ( i != params1d->cells_number - 1 && INNER == status[i + 1] 
            && OUTER == status[i - 1] )
                status[i + 1] = BOUNDARY;
            else
            {
                printf( "\ndefine_status -> wrong status definition\n" );
                system ( "Pause" );
            }
        }
    }
}

void calc_status ( struct ParametersCommon *paramsc, struct Parameters1d *params1d, double coordinate_of_left_boundary_of_body, double coordinate_of_right_boundary_of_body,
    int *status, int *index_of_left_ghost_cell, int *index_of_right_ghost_cell, double **v_ncons,double v_left[M], double v_right[M] )
{
    double grid_step = ( params1d->right_boundary_x - params1d->left_boundary_x )
        / params1d->cells_number;

    // левая граница тела сдвинулась направо в ячейку i + 1
    if ( coordinate_of_left_boundary_of_body >= ( *index_of_left_ghost_cell + 1 ) * grid_step )
    {
        printf("\n right \n %d %lf %lf %lf", *index_of_left_ghost_cell, grid_step, coordinate_of_left_boundary_of_body, v_left[P_GAS]);

        status[*index_of_left_ghost_cell + 1] = GHOST;
        status[*index_of_left_ghost_cell] = BOUNDARY;
        status[*index_of_left_ghost_cell - 1] = INNER;
        v_ncons[ *index_of_left_ghost_cell ][P_GAS] = v_left[P_GAS];
        v_ncons[ *index_of_left_ghost_cell ][V_GAS] = v_left[V_GAS];
        v_ncons[ *index_of_left_ghost_cell ][R_GAS] = v_left[R_GAS];
        v_ncons[ *index_of_left_ghost_cell ][P_DISP] = v_left[P_DISP];
        v_ncons[ *index_of_left_ghost_cell ][V_DISP] = v_left[V_DISP];
        v_ncons[ *index_of_left_ghost_cell ][R_DISP] = v_left[R_DISP];
        v_ncons[ *index_of_left_ghost_cell ][B_DISP] =  0.8 * v_left[B_DISP];
        *index_of_left_ghost_cell += 1;
    }
    else if ( coordinate_of_left_boundary_of_body >= *index_of_left_ghost_cell * grid_step ) // осталась в той же ячейке
    {}
    else // сдвинулась налево в ячейку i - 1
    {
        status[*index_of_left_ghost_cell - 2] = BOUNDARY;
        status[*index_of_left_ghost_cell - 1] = GHOST;
        status[*index_of_left_ghost_cell] = OUTER;
        v_ncons[ *index_of_left_ghost_cell ][P_GAS] = 0;
        v_ncons[ *index_of_left_ghost_cell ][V_GAS] = 0;
        v_ncons[ *index_of_left_ghost_cell ][R_GAS] = 0;
        *index_of_left_ghost_cell -= 1;
    }

    // правая граница тела сдвинулась направо в ячейку i + 1
    if ( coordinate_of_right_boundary_of_body >= ( *index_of_right_ghost_cell + 1 ) * grid_step ) // сдвинулась направо в ячейку i + 1
    {
        status[*index_of_right_ghost_cell + 2] = BOUNDARY;
        status[*index_of_right_ghost_cell + 1] = GHOST;
        status[*index_of_right_ghost_cell] = OUTER;
        v_ncons[ *index_of_right_ghost_cell ][P_GAS] = 1e5;
        v_ncons[ *index_of_right_ghost_cell ][V_GAS] = 0;
        v_ncons[ *index_of_right_ghost_cell ][R_GAS] = 1.2;
        v_ncons[ *index_of_right_ghost_cell ][P_DISP] = 1e5;
        v_ncons[ *index_of_right_ghost_cell ][V_DISP] = 0;
        v_ncons[ *index_of_right_ghost_cell ][R_DISP] = 1578;
        *index_of_right_ghost_cell += 1;
    }
    else if ( coordinate_of_right_boundary_of_body >= *index_of_right_ghost_cell * grid_step ) // осталась в той же ячейке
    {}
    else // сдвинулась налево в ячейку i - 1
    {
        status[*index_of_right_ghost_cell - 1] = GHOST;
        status[*index_of_right_ghost_cell] = BOUNDARY;
        status[*index_of_right_ghost_cell + 1] = INNER;
        v_ncons[ *index_of_right_ghost_cell ][P_GAS] = v_right[P_GAS];
        v_ncons[ *index_of_right_ghost_cell ][V_GAS] = v_right[V_GAS];
        v_ncons[ *index_of_right_ghost_cell ][R_GAS] = v_right[R_GAS];
        v_ncons[ *index_of_right_ghost_cell ][P_DISP] = v_right[P_DISP];
        v_ncons[ *index_of_right_ghost_cell ][V_DISP] = v_right[V_DISP];
        v_ncons[ *index_of_right_ghost_cell ][R_DISP] = v_right[R_DISP];
        v_ncons[ *index_of_right_ghost_cell ][B_DISP] = v_right[B_DISP];
        *index_of_right_ghost_cell -= 1;
    }
}
void convert_2D_to_1D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector){
    u_1D_vector[B_DISP] = u_2D_vector[B_DISP_2D];
    u_1D_vector[R_DISP] = u_2D_vector[R_DISP_2D];
    u_1D_vector[P_DISP] = u_2D_vector[P_DISP_2D];
    u_1D_vector[R_GAS] = u_2D_vector[R_GAS_2D];
    u_1D_vector[P_GAS] = u_2D_vector[P_GAS_2D];
    switch (dir){
        case X_DIRECTION:
            u_1D_vector[V_DISP] = u_2D_vector[V_DISP_2D];
            u_1D_vector[V_GAS] = u_2D_vector[V_GAS_2D];
            break;
        case Y_DIRECTION:
            u_1D_vector[V_DISP] = u_2D_vector[U_DISP_2D];
            u_1D_vector[V_GAS] = u_2D_vector[U_GAS_2D];
            break;
        default:
            printf( "\nconvert_2D_to_1D -> wrong direction - only x- or y- are possible\n\n" );
            exit( EXIT_FAILURE );
    }
}

void convert_1D_to_2D_vector(double *u_2D_vector, Direction2d dir, double *u_1D_vector){
    u_2D_vector[B_DISP_2D] = u_1D_vector[B_DISP];
    u_2D_vector[R_DISP_2D] = u_1D_vector[R_DISP];
    u_2D_vector[P_DISP_2D] = u_1D_vector[P_DISP];
    u_2D_vector[R_GAS_2D] = u_1D_vector[R_GAS];
    u_2D_vector[P_GAS_2D] = u_1D_vector[P_GAS];
    switch (dir){
        case X_DIRECTION:
            u_2D_vector[V_DISP_2D] = u_1D_vector[V_DISP];
            u_2D_vector[V_GAS_2D] = u_1D_vector[V_GAS];
            break;
        case Y_DIRECTION:
            u_2D_vector[U_DISP_2D] = u_1D_vector[V_DISP];
            u_2D_vector[U_GAS_2D] = u_1D_vector[V_GAS];
            break;
        default:
            printf( "\nconvert_2D_to_1D -> wrong direction - only x- or y- are possible\n\n" );
            exit( EXIT_FAILURE );
    }
}


// Инициализация вектора-решения для одномерной задачи
// paramsc - структура с основными параметрами вычислительного эксперимента (in)
// files_directory - директория, в которой находятся все внешние файлы, требуемые для расчета (in)
// time_mom - структура, определяющая текущий момент времени (out)
// **initial_solution - массив векторов в примитивных переменных - начальных условий в каждой ячейке (out)
void init_solution_1d( struct ParametersCommon *paramsc, struct Parameters1d *params1d, char *files_directory, struct TimeMoment *time_mom, double **initial_solution ) {

    FILE *restart_info_file = NULL; // дескриптор файла с информацией о последнем доступном файле с промежуточными результатами
    int i_cell_global = 0; // глобальный номер ячейки
    
    /* попытка считать файл с информацией о рестарте, определяется значение params->isContinue */
/*    try_to_restart( paramsc, files_directory, time_mom, restart_info_file );

    if ( paramsc->isContinue ) {
        // запускаемый расчет является продолжением предыдущего
        read_solution( paramsc, params1d, restart_info_file, initial_solution );
    }
    else { */
        /* инициализация полей в расчетной области */
        for ( int i_block = 0; i_block < params1d->ic_blocks_number; i_block++ ) { // цикл по блокам
            for ( int i_cell = params1d->cell_begin[i_block]; i_cell <= params1d->cell_end[i_block]; i_cell++ ) { // цикл по ячейкам внутри блока
                initial_solution[i_cell_global][B_DISP] = params1d->block_values[i_block][B_DISP];
                initial_solution[i_cell_global][R_DISP] = params1d->block_values[i_block][R_DISP];
                initial_solution[i_cell_global][V_DISP] = params1d->block_values[i_block][V_DISP];
                initial_solution[i_cell_global][P_DISP] = params1d->block_values[i_block][P_DISP];
                initial_solution[i_cell_global][R_GAS] = params1d->block_values[i_block][R_GAS];
                initial_solution[i_cell_global][V_GAS] = params1d->block_values[i_block][V_GAS];
                initial_solution[i_cell_global][P_GAS] = params1d->block_values[i_block][P_GAS];
                for (int i_scalar = 0; i_scalar < params1d->number_of_scalars; i_scalar++)
                    initial_solution[i_cell_global][Z0 + i_scalar] = params1d->scalar_values[i_block][i_scalar];
                i_cell_global++;
            }
        }
 //   }
    
}


// Обработка граничного условия
// params1d - структура с параметрами одномерной задачи (in)
// v_ncons[M] - вектор примитивных переменных (in)
// boun_type - тип граничного условия (in)
// boun_v[M] - вектор примитивных переменных в фиктивной ячейке (out)
// curr_t - текущий момент времени (in)
// curr_inflow_parameters[M] - текущие значения параметров втекания(in)
void boundary( struct Parameters1d *params1d, struct Parameters2d *params2d, double v_ncons[M], int boun_type, double boun_v[M], double curr_t, double curr_inflow_parameters[M] ) {

    switch ( boun_type ) {
        case WALL:
            boun_v[B_DISP] = v_ncons[B_DISP]; // объемная доля дисперсной фазы
	    boun_v[R_DISP] = v_ncons[R_DISP]; // плотность дисперсной фазы
	    boun_v[V_DISP] = - v_ncons[V_DISP]; // скорость дисперсной фазы
            boun_v[P_DISP] = v_ncons[P_DISP]; // давление дисперсной фазы
            boun_v[R_GAS] = v_ncons[R_GAS]; // плотность газовой фазы
            boun_v[V_GAS] = - v_ncons[V_GAS]; // скорость газовой фазы
            boun_v[P_GAS] = v_ncons[P_GAS]; // давление газовой фазы
	    break;
        case FREE:
            boun_v[B_DISP] = v_ncons[B_DISP];   /* объемная доля дисперсной фазы */
	    boun_v[R_DISP] = v_ncons[R_DISP];   /* плотность дисперсной фазы */ 
	    boun_v[V_DISP] = v_ncons[V_DISP];   /* скорость дисперсной фазы */
            boun_v[P_DISP] = v_ncons[P_DISP];   /* давление дисперсной фазы */
            boun_v[R_GAS] = v_ncons[R_GAS];   /* плотность газовой фазы */
            boun_v[V_GAS] = v_ncons[V_GAS];   /* скорость газовой фазы */
            boun_v[P_GAS] = v_ncons[P_GAS];   /* давление газовой фазы */
            break;
        case INFLOW:
            boun_v[B_DISP] = curr_inflow_parameters[B_DISP]; /* объемная доля дисперсной фазы */
	    boun_v[R_DISP] = curr_inflow_parameters[R_DISP]; /* плотность дисперсной фазы */ 
	    boun_v[V_DISP] = curr_inflow_parameters[V_DISP]; /* скорость дисперсной фазы */
            boun_v[P_DISP] = curr_inflow_parameters[P_DISP]; /* давление дисперсной фазы */
            boun_v[R_GAS] = curr_inflow_parameters[R_GAS]; /* плотность газовой фазы */
            boun_v[V_GAS] = curr_inflow_parameters[V_GAS]; /* скорость газовой фазы */
            boun_v[P_GAS] = curr_inflow_parameters[P_GAS]; /* давление газовой фазы */
            break;
        case COMPLEX_INFLOW:
            if ( curr_t < params1d->change_inflow_time ) {
                boun_v[B_DISP] = params1d->inflow_values[B_DISP]; /* объемная доля дисперсной фазы */
	        boun_v[R_DISP] = params1d->inflow_values[R_DISP]; /* плотность дисперсной фазы */ 
	        boun_v[V_DISP] = params1d->inflow_values[V_DISP]; /* скорость дисперсной фазы */
                boun_v[P_DISP] = params1d->inflow_values[P_DISP]; /* давление дисперсной фазы */
                boun_v[R_GAS] = params1d->inflow_values[R_GAS]; /* плотность газовой фазы */
                boun_v[V_GAS] = params1d->inflow_values[V_GAS]; /* скорость газовой фазы */
                boun_v[P_GAS] = params1d->inflow_values[P_GAS]; /* давление газовой фазы */
            }
            else {
                boun_v[B_DISP] = params1d->changed_inflow_values[B_DISP]; /* объемная доля дисперсной фазы */
	        boun_v[R_DISP] = params1d->changed_inflow_values[R_DISP]; /* плотность дисперсной фазы */ 
	        boun_v[V_DISP] = params1d->changed_inflow_values[V_DISP]; /* скорость дисперсной фазы */
                boun_v[P_DISP] = params1d->changed_inflow_values[P_DISP]; /* давление дисперсной фазы */
                boun_v[R_GAS] = params1d->changed_inflow_values[R_GAS]; /* плотность газовой фазы */
                boun_v[V_GAS] = params1d->changed_inflow_values[V_GAS]; /* скорость газовой фазы */
                boun_v[P_GAS] = params1d->changed_inflow_values[P_GAS]; /* давление газовой фазы */
            }
            break;
        default:
            printf( "\nboundary -> wrong boundary condition.\n" );
            exit( EXIT_FAILURE );
            break;
    }
    for (int i = 0; i < params1d->number_of_scalars; i++){
        if ( boun_type == INFLOW )
            boun_v[Z0 + i] = curr_inflow_parameters[Z0 + i];
        else
            boun_v[Z0 + i] = v_ncons[Z0 + i];
    }
}

void current_inflow_parameters( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d, int boun_direction, double *curr_inflow_params){

    if ( paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC ){
        
        for (int i = 0; i < M1D + params1d->number_of_scalars; i++)
            curr_inflow_params[i] = params1d->inflow_values[i];

    }
    else if ( paramsc->program_name == TWOD2PHC ){
        if (boun_direction == LEFT_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_left_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_left_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_left_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_left_bc[P_DISP_2D];
            curr_inflow_params[R_GAS] = params2d->inflow_values_left_bc[R_GAS_2D];
            curr_inflow_params[V_GAS] = params2d->inflow_values_left_bc[V_GAS_2D];
            curr_inflow_params[P_GAS] = params2d->inflow_values_left_bc[P_GAS_2D]; 
        }
        else if (boun_direction == RIGHT_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_right_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_right_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_right_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_right_bc[P_DISP_2D];
            curr_inflow_params[R_GAS] = params2d->inflow_values_right_bc[R_GAS_2D];
            curr_inflow_params[V_GAS] = params2d->inflow_values_right_bc[V_GAS_2D];
            curr_inflow_params[P_GAS] = params2d->inflow_values_right_bc[P_GAS_2D];    
        }
        else if (boun_direction == UP_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_up_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_up_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_up_bc[U_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_up_bc[P_DISP_2D];
            curr_inflow_params[R_GAS] = params2d->inflow_values_up_bc[R_GAS_2D];
            curr_inflow_params[V_GAS] = params2d->inflow_values_up_bc[U_GAS_2D];
            curr_inflow_params[P_GAS] = params2d->inflow_values_up_bc[P_GAS_2D];
        }
        else if (boun_direction == DOWN_BOUNDARY){
            curr_inflow_params[B_DISP] = params2d->inflow_values_down_bc[B_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_down_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_down_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_down_bc[P_DISP_2D];
            curr_inflow_params[R_DISP] = params2d->inflow_values_down_bc[R_DISP_2D];
            curr_inflow_params[V_DISP] = params2d->inflow_values_down_bc[V_DISP_2D];
            curr_inflow_params[P_DISP] = params2d->inflow_values_down_bc[P_DISP_2D]; 
        }
        else{
            
            printf_s("Wrong boun_direction number - %d. Only %d, %d and %d, %d are allowed.\n", boun_direction, LEFT_BOUNDARY, RIGHT_BOUNDARY, UP_BOUNDARY, DOWN_BOUNDARY);
            exit(EXIT_FAILURE);

        }
    }
    else{
        
        printf_s("Wrong program name number - %d. Only %d, %d and %d are allowed.\n", paramsc->program_name, ONED2PHC, ONED3PHC, TWOD2PHC);
        exit(EXIT_FAILURE);

    }

}

// Расчет одного шага по времени в одной ячейке
// paramsc - структура с основными параметрами вычислительного эксперимента
// debug_info - структура с отладочной информацией
// u_left - вектор примитивных переменных в ячейке слева от рассчитываемой
// u_center - вектор примитивных переменных в рассчитываемой ячейке
// u_right - вектор примитивных переменных в ячейке справа от рассчитываемой
// slopes_left - вектор наклонов в ячейке слева от рассчитываемой
// slopes_center - вектор наклонов в рассчитываемой ячейке
// slopes_right - вектор наклонов в ячейке справа от рассчитываемой
// dt - шаг интегрирования по времени
// h - размер ячейки
// u_next - вектор примитивных переменных в рассчитываемой ячейке на следующем временном шаге
// n - реальный размер векторов
// is_pressure_relaxation_now - true, если после решения одномерной задачи по данному направлению нужно применять релаксацию давлений; false - иначе
//                              если в данной задачу вообще не нужно применять релаксацию давлений, то true не повлияет на результат, так как данная опция контролируется в parameters.dat
// number_of_scalars - число лагранжевых скаляров в данной задаче
// curr_time - текущий момент времени
// configuration_pressure - конфигурационное давление
void calc_step_in_cell( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct DebugInfo *debug_info, double u_left[M], double u_center[M],
                        double u_right[M], double slopes_left[M], double slopes_center[M], double slopes_right[M],
                        double dt, double h, double u_next[M], int step_number, int n, bool is_pressure_relaxation_now, int number_of_scalars, double curr_time, double *configuration_pressure, double body_velocity, int i, int *status, double cont_left[M], double cont_right[M] ) {
    
    switch ( paramsc->numerical_method ) {
        case CIR1:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR_DISP scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_1( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case CIR2:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR_GAS scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_2( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case CIR3:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR3 scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_3( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case CIR4:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> CIR4 scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            cir_4( paramsc, u_left, u_center, u_right, dt, h, u_next );
            break;
        case GODUNOV:
            if ( paramsc->media_model != BAER_NUNZIATO ) {
                printf( "\ncalc_step_in_cell -> Godunov scheme is realized for BN equations only\n" );
                exit( EXIT_FAILURE );
            }
            godunov( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, is_pressure_relaxation_now, curr_time, configuration_pressure  );
            break;
        case HLL:
            hll( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, is_pressure_relaxation_now, params1d->number_of_scalars, curr_time, configuration_pressure, body_velocity, i, status,  cont_left, cont_right );
            break;
        case RUSANOV:
            rusanov_1d( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, curr_time, configuration_pressure  );
            break;
        case HLLC:
            hllc_1d( paramsc, params1d, debug_info, u_left, u_center, u_right, slopes_left, slopes_center, slopes_right, dt, h, u_next, step_number, n, is_pressure_relaxation_now, params1d->number_of_scalars, curr_time, configuration_pressure, body_velocity, i, status,  cont_left, cont_right );
            break;
        default:
            printf( "\ncalc_step_in_cell -> wrong scheme number\n" );
            exit( EXIT_FAILURE );
        }
    //printf("\n%lf %lf %lf \n", u_next[P_GAS]);
}


// Инициализация структуры для хранения и передачи отладочной информации
// output_file_directory - директория, в которую должны быть записаны файлы
// debug_info - структура с отладочной информацией
void init_debug_info( const char *output_file_directory, struct DebugInfo *debug_info ) {

    (*debug_info).current_cell = -1;
    (*debug_info).current_cell_x = -1.0;
    (*debug_info).neighbour_cell = -1;
    for ( int i_component = 0; i_component < M; i_component++ ) {
        (*debug_info).current_cell_vncons[i_component] = -1.0;
        (*debug_info).neighbour_cell_vncons[i_component] = -1.0;
    }
    
    for ( int i_relax_case = 0; i_relax_case < RELAX_CASES_NUM; i_relax_case++ )
        (*debug_info).relaxation_cases[i_relax_case] = 0;
    for ( int i_godunov_case = 0; i_godunov_case < GODUNOV_CASES_NUM; i_godunov_case++ )
        (*debug_info).godunov_cases[i_godunov_case] = 0;

    char string[MAX_STRING_SIZE];

    strcpy_s( string, output_file_directory );
    strcat_s( string, "\\relaxation.dat" );
    if ( ( fopen_s( &(*debug_info).relaxation_out, string, "wt" ) ) != 0 ) {
        printf( "\nread_Parameters1d -> can't open file %s for format writing\n\n", string );
        exit( EXIT_FAILURE );
    }

    strcpy_s( string, output_file_directory );
    strcat_s( string, "\\godunov.dat" );
    if ( ( fopen_s( &(*debug_info).godunov_out, string, "wt" ) ) != 0 ) {
        printf( "\nread_Parameters1d -> can't open file %s for format writing\n\n", string );
        exit( EXIT_FAILURE );
    }

}

// Расчет усредненного межфазного давления
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор неконсервативных переменных в рассматриваемой ячейке
// Возвращает искомое давление
double calc_p_i( const struct ParametersCommon* paramsc, const double v_ncons[M] ) {

    switch ( paramsc->media_model ) {
        case BAER_NUNZIATO:
            return v_ncons[P_GAS];
        case SAUREL_ABGRALL:
            return ( v_ncons[B_DISP] * v_ncons[P_DISP] + ( 1.0 - v_ncons[B_DISP] ) * v_ncons[P_GAS] );
        default:
            printf( "calc_p_i -> wrong value of params->media_model\n" );
            exit( EXIT_FAILURE );
    }

}

// Расчет усредненной межфазной плотности, имеет смысл только для модели Saurel-Abgrall
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор неконсервативных переменных в рассматриваемой ячейке
// Возвращает искомую плотность
double calc_r_i( const struct ParametersCommon* paramsc, const double v_ncons[M] ) {

    switch ( paramsc->media_model ) {
        case SAUREL_ABGRALL:
            return ( v_ncons[B_DISP] * v_ncons[R_DISP] + ( 1.0 - v_ncons[B_DISP] ) * v_ncons[R_GAS] );
        case BAER_NUNZIATO:
        default:
            printf( "calc_r_i -> wrong value of params->media_model\n" );
            exit( EXIT_FAILURE );
    }

}

// Расчет скоростей звука для обеих фаз по полному вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// c1 - скорость звука в дисперсной фазе
// c2 - скорость звука в газовой фазе
void calc_sound_velocity( const struct ParametersCommon* paramsc, const double v_ncons[M], double* c1, double* c2 ) {

    if (v_ncons[P_DISP] < 0.0 || v_ncons[P_GAS] < 0.0){
        if (v_ncons[P_DISP] < 0.0)
            printf("\ncalc_sound_velocity -> pressure of solid phase is negative\n");
            printf("\ncalc_sound_velocity -> pressure of solid phase is %lf\n", v_ncons[P_DISP]);
            printf("\ncalc_sound_velocity -> volume_fraction of solid phase is %lf\n", v_ncons[B_DISP]);
        if (v_ncons[P_GAS] < 0.0)
            printf("\ncalc_sound_velocity -> pressure of gas phase is negative\n");
        exit(EXIT_FAILURE);
    }

    *c1 = sqrt( paramsc->g1 * ( v_ncons[P_DISP] + paramsc->p01 ) / v_ncons[R_DISP] );
    *c2 = sqrt( paramsc->g2 * ( v_ncons[P_GAS] + paramsc->p02 ) / v_ncons[R_GAS] * 
        (1.0 + paramsc->b_virial * v_ncons[R_GAS] - paramsc->b_virial * v_ncons[R_GAS] * paramsc->b_virial * v_ncons[R_GAS] / paramsc->g2 / (1.0 + paramsc->b_virial * v_ncons[R_GAS]) ) );

}

// Расчет скоростей звука для заданной фазы по вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// phase - идентификатор фазы, для которой рассчитывается скорость звука - газовая или дисперсная
// Возвращает скорость звука
double calc_sound_velocity_one_phase( const struct ParametersCommon* paramsc, const double v_ncons[M], const Phase phase ) {

    if ( phase == GAS_PHASE ){
        if ( v_ncons[P_GAS] < 0.0){
            printf("\ncalc_sound_velocity_one_phase -> pressure of gas phase is negative %lf \n", v_ncons[P_GAS]);
            exit(EXIT_FAILURE);
        }
        return sqrt( paramsc->g2 * ( v_ncons[P_GAS] + paramsc->p02 ) / v_ncons[R_GAS] * 
        (1.0 + paramsc->b_virial * v_ncons[R_GAS] - paramsc->b_virial * v_ncons[R_GAS] * paramsc->b_virial * v_ncons[R_GAS] / paramsc->g2 / (1.0 + paramsc->b_virial * v_ncons[R_GAS]) ) );
    }
    else{
        if ( v_ncons[P_DISP] < 0.0){
            printf("\ncalc_sound_velocity_one_phase -> pressure of solid phase is negative\n");
            exit(EXIT_FAILURE);
        }
        return sqrt( paramsc->g1 * ( v_ncons[P_DISP] + paramsc->p01 ) / v_ncons[R_DISP] );
    }
}

// Расчет скорости звука для одной фазы по однофазному вектору переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// phase - идентификатор фазы, для которой рассчитывается скорость звука - газовая или дисперсная
// Возвращает скорость звука
double calc_sound_velocity_reduced( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase ) {

    if ( phase == GAS_PHASE )
        return sqrt( paramsc->g2 * v_ncons[P] / v_ncons[R] );
    else
        return sqrt( paramsc->g1 * ( v_ncons[P] + paramsc->p01 ) / v_ncons[R] );

}

// Расчет усредненной межфазной скорости
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор неконсервативных переменных в рассматриваемой ячейке
// Возвращает искомую скорость
double calc_u_i( const struct ParametersCommon* paramsc, const double v_ncons[M] ) {

    double sum1, sum2; 
    
    switch ( paramsc->media_model ) {
        case BAER_NUNZIATO:
            return v_ncons[V_DISP];
        case SAUREL_ABGRALL:
            sum1 = v_ncons[B_DISP] * v_ncons[R_DISP] * v_ncons[V_DISP] + ( 1.0 - v_ncons[B_DISP] ) * v_ncons[R_GAS] * v_ncons[V_GAS];
            sum2 = calc_r_i( paramsc, v_ncons );
            return ( sum1 / sum2 );
        default:
            printf( "calc_u_i -> wrong value of params->media_model\n" );
            exit( EXIT_FAILURE );
    }

}

// Преобразование вектора "консервативных" переменных в вектор примитивных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// v_ncons - вектор примитивных переменных
// Возвращает: SUCCEESS            давление в обеих фазах положительно
//             NEGATIVE_PRESSURE   давление хотя бы в одной из фаз отрицательно
ReturnCodes convert_cons_to_noncons( const struct ParametersCommon* paramsc, const double v_cons[M], double v_ncons[M], int number_of_scalars ) {

    double b2 = 1 - v_cons[B_DISP]; // объемная доля газовой фазы
    ReturnCodes return_code = SUCCESS; // код возврата функции

    v_ncons[B_DISP] = v_cons[B_DISP]; // b1, объемная доля дисперсной фазы
    v_ncons[R_DISP] = v_cons[R_DISP] / v_cons[B_DISP]; // r1, плотность дисперсной фазы
    v_ncons[V_DISP] = v_cons[V_DISP] / v_cons[R_DISP]; // v1, скорость дисперсной фазы
    v_ncons[P_DISP] = ( paramsc->g1 - 1.0 ) * ( v_cons[P_DISP] / v_cons[B_DISP] - 0.5 * pow( v_cons[V_DISP], 2.0 ) / v_cons[B_DISP] / v_cons[R_DISP] ) -
        paramsc->g1 * paramsc->p01; // p1, давление дисперсной фазы
    v_ncons[R_GAS] = v_cons[R_GAS] / b2; // r2, плотность газовой фазы
    v_ncons[V_GAS] = v_cons[V_GAS] / v_cons[R_GAS]; // v2, скорость газовой фазы
    v_ncons[P_GAS] = ( paramsc->g2 - 1.0 ) * (1.0 + paramsc->b_virial * v_cons[R_GAS] / b2) * ( v_cons[P_GAS] / b2 - 0.5 * pow( v_cons[V_GAS], 2.0 ) / b2 / v_cons[R_GAS] ) -
        paramsc->g2 * paramsc->p02; // p2, давление газовой фазы 

    for (int i = 0; i < number_of_scalars; i++)
        v_ncons[Z0 + i] = v_cons[Z0 + i] / v_cons[R_DISP];

    return return_code;

}

/* Преобразование редуцированного однофазного вектора "консервативных" переменных в вектор примитивных

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_cons[M_REDUCTION] - вектор "консервативных" переменных (in)
   phase - идентификатор фазы, для которой осуществляется преобразование - газовая или дисперсная (in)

   v_ncons[M_REDUCTION] - вектор примитивных переменных (out) */
void convert_cons_to_noncons_reduction( struct ParametersCommon *paramsc, double v_cons[M_REDUCTION], Phase phase, double v_ncons[M_REDUCTION] ) {

    /* плотность */
    v_ncons[R] = v_cons[R];
    
    /* скорость */
    v_ncons[V] = v_cons[V] / v_cons[R];
    
    /* давление */
    switch ( phase ) {
        case GAS_PHASE:
            v_ncons[P] = ( paramsc->g2 - 1.0 ) * ( v_cons[P] - 0.5 * pow( v_cons[V], 2.0 ) / v_cons[R] );
            break;
        case DISPERSED_PHASE:
            v_ncons[P] = ( paramsc->g1 - 1.0 ) * ( v_cons[P] - 0.5 * pow( v_cons[V], 2.0 ) / v_cons[R] ) -
                paramsc->g1 * paramsc->p01;
            break;
        default:
            printf( "\nconvert_cons_to_noncons_reduction -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

/* Формирование однофазного вектора по полному двухфазному вектору
 
   full_vector[M] - полный двухфазный вектор (in)
   phase - идентификатор фазы, для которой по полному вектору строится редуцированный - газовая или дисперсная (in)
 
   reduced_vector[M_REDUCTION] - редуцированный однофазный вектор (out) */
void convert_full_to_reduced( double full_vector[M], Phase phase, double reduced_vector[M_REDUCTION] ) {

    switch ( phase ) {
        case GAS_PHASE:
            reduced_vector[R] = full_vector[R_GAS];
            reduced_vector[V] = full_vector[V_GAS];
            reduced_vector[P] = full_vector[P_GAS];
            break;
        case DISPERSED_PHASE:
            reduced_vector[R] = full_vector[R_DISP];
            reduced_vector[V] = full_vector[V_DISP];
            reduced_vector[P] = full_vector[P_DISP];
            break;
        default:
            printf( "\nconvert_full_to_reduced -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

// Преобразование вектора примитивных переменных в вектор "консервативных"
// params - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// v_cons - вектор "консервативных" переменных
// number_of_scalars - количество лагранжевых скаляров
void convert_noncons_to_cons( const struct ParametersCommon *paramsc, const double v_ncons[M], double v_cons[M], int number_of_scalars ) {

    double b2 = 1 - v_ncons[B_DISP]; // объемная доля газовой фазы 

    v_cons[B_DISP] = v_ncons[B_DISP]; // b1, объемная доля дисперсной фазы 
    v_cons[R_DISP] = v_ncons[B_DISP] * v_ncons[R_DISP]; // b1 * r1
    v_cons[V_DISP] = v_cons[R_DISP] * v_ncons[V_DISP]; // b1 * r1 * v1
    v_cons[P_DISP] = v_ncons[B_DISP] * v_ncons[R_DISP] * ( 0.5 * pow( v_ncons[V_DISP], 2.0 ) + v_ncons[P_DISP] / v_ncons[R_DISP] /
        ( paramsc->g1 - 1.0 ) + paramsc->g1 * paramsc->p01 / ( paramsc->g1 - 1.0 ) / v_ncons[R_DISP] ); // b1 * r1 * E1 
    v_cons[R_GAS] = b2 * v_ncons[R_GAS]; // b2 * r2
    v_cons[V_GAS] = v_cons[R_GAS] * v_ncons[V_GAS]; // b2 * r2* v2
    v_cons[P_GAS] = v_cons[R_GAS] * ( 0.5 * pow( v_ncons[V_GAS], 2.0 ) + ( v_ncons[P_GAS] + paramsc->g2 * paramsc->p02 )/ v_ncons[R_GAS] /
		( paramsc->g2 - 1.0 ) / (1.0 + paramsc->b_virial * v_ncons[R_GAS])) ; // b2 * r2 * E2

    for (int i = 0; i < number_of_scalars; i++)
        v_cons[Z0 + i] = v_ncons[Z0 + i] * v_ncons[R_DISP] * v_ncons[B_DISP];
}


/* Преобразование однофазного вектора примитивных переменных в вектор "консервативных"

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons[M_REDUCTION] - вектор примитивных переменных (in)
   phase - идентификатор фазы, для которой осуществляется преобразование - газовая или дисперсная (in)

   v_cons[M_REDUCTION] - вектор "консервативных" переменных (out) */
void convert_noncons_to_cons_reduction( struct ParametersCommon *paramsc, double v_ncons[M_REDUCTION], Phase phase, double v_cons[M_REDUCTION] ) {

    /* масса */
    v_cons[R] = v_ncons[R];
    
    /* импульс */
    v_cons[V] = v_ncons[R] * v_ncons[V];

    /* энергия */
    switch ( phase ) {
        case GAS_PHASE:
            v_cons[P] = v_ncons[R] * ( 0.5 * pow( v_ncons[V], 2.0 ) + v_ncons[P] / v_ncons[R] / ( paramsc->g2 - 1.0 ) );
            break;
        case DISPERSED_PHASE:
            v_cons[P] = v_ncons[R] * ( 0.5 * pow( v_ncons[V], 2.0 ) + v_ncons[P] / v_ncons[R] /
                ( paramsc->g1 - 1.0 ) + paramsc->g1 * paramsc->p01 / ( paramsc->g1 - 1.0 ) / v_ncons[R] );
            break;
        default:
            printf( "\nconvert_noncons_to_cons_reduction -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

// Формирование части полного двухфазного вектора по сокращенному однофазному
// reduced_vector[M_REDUCTION] - редуцированный однофазный вектор (in)
// phase - идентификатор фазы, для которой по редуцированному вектору строится часть полного - газовая или дисперсная (in)
// full_vector[M] - полный двухфазный вектор (out)
void convert_reduced_to_full( double reduced_vector[M], Phase phase, double full_vector[M_REDUCTION] ) {

    switch ( phase ) {
        case GAS_PHASE:
            full_vector[R_GAS] = reduced_vector[R];
            full_vector[V_GAS] = reduced_vector[V];
            full_vector[P_GAS] = reduced_vector[P];
            break;
        case DISPERSED_PHASE:
            full_vector[R_DISP] = reduced_vector[R];
            full_vector[V_DISP] = reduced_vector[V];
            full_vector[P_DISP] = reduced_vector[P];
            break;
        default:
            printf( "\nconvert_reduced_to_full -> wrong value of phase variable.\n" );
            exit( EXIT_FAILURE );
    }

}

// Расчет вектора дифференциального "потока" по вектору "консервативных" переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// flux - рассчитываемый вектор дифференциального "потока"
void diff_flux_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* flux, int number_of_scalars ) {

    double v_ncons[M]; // вектор примитивных переменных
        
    convert_cons_to_noncons( paramsc, v_cons, v_ncons, number_of_scalars );

    double b2 = 1 - v_ncons[B_DISP]; // объемная доля газовой фазы

    (*flux)[B_DISP] = 0.0; // поток объемной доли дисперсной фазы
    (*flux)[R_DISP] = v_cons[V_DISP]; // поток массы в дисперсной фазе
    (*flux)[V_DISP] = v_ncons[B_DISP] * ( v_ncons[R_DISP] * pow( v_ncons[V_DISP], 2.0 ) + v_ncons[P_DISP] ); // поток импульса в дисперсной фазе
    (*flux)[P_DISP] = v_ncons[V_DISP] * ( v_cons[P_DISP] + v_ncons[B_DISP] * v_ncons[P_DISP] ); // поток энергии в дисперсной фазе
    (*flux)[R_GAS] = v_cons[V_GAS]; // поток массы в газовой фазе
    (*flux)[V_GAS] = b2 * ( v_ncons[R_GAS] * pow( v_ncons[V_GAS], 2.0 ) + v_ncons[P_GAS] ); // поток импульса в газовой фазе
    (*flux)[P_GAS] = v_ncons[V_GAS] * ( v_cons[P_GAS] + b2 * v_ncons[P_GAS] ); // поток энергии в газовой фазе

    for ( int i = 0; i < number_of_scalars; i++ )
        (*flux)[Z0 + i] = v_cons[V_DISP] * v_ncons[Z0 + i];
}

// Расчет вектора дифференциального "потока" по вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_ncons - вектор примитивных переменных
// flux - рассчитываемый вектор дифференциального "потока"
void diff_flux_ncons( const struct ParametersCommon* paramsc, const double v_ncons[M], array1D* flux, int number_of_scalars ) {

    double v_cons[M]; // вектор "консервативных" переменных

    convert_noncons_to_cons( paramsc, v_ncons, v_cons, number_of_scalars );

    double b2 = 1 - v_ncons[B_DISP]; // объемная доля газовой фазы
    
    (*flux)[B_DISP] = 0.0; // поток объемной доли дисперсной фазы
    (*flux)[R_DISP] = v_cons[V_DISP]; // поток массы в дисперсной фазе
    (*flux)[V_DISP] = v_ncons[B_DISP] * ( v_ncons[R_DISP] * pow( v_ncons[V_DISP], 2.0 ) + v_ncons[P_DISP] ); // поток импульса в дисперсной фазе
    (*flux)[P_DISP] = v_ncons[V_DISP] * ( v_cons[P_DISP] + v_ncons[B_DISP] * v_ncons[P_DISP] ); // поток энергии в дисперсной фазе
    (*flux)[R_GAS] = v_cons[V_GAS]; // поток массы в газовой фазе
    (*flux)[V_GAS] = b2 * ( v_ncons[R_GAS] * pow( v_ncons[V_GAS], 2.0 ) + v_ncons[P_GAS] ); // поток импульса в газовой фазе
    (*flux)[P_GAS] = v_ncons[V_GAS] * ( v_cons[P_GAS] + b2 * v_ncons[P_GAS] ); // поток энергии в газовой фазе

    for ( int i = 0; i < number_of_scalars; i++ )
        (*flux)[Z0 + i] = v_cons[V_DISP] * v_ncons[Z0 + i];
}

/* Расчет однофазного вектора дифференциального "потока" по сокращенному однофазному вектору примитивных переменных

   paramsc - структура с основными параметрами вычислительного эксперимента (in)
   v_ncons - сокращенный однофазный вектор примитивных переменных (in)
   phase - фаза - газовая или дисперсная, для которой рассчитывается "поток" (in)

   flux - рассчитываемый сокращенный однофазный вектор дифференциального "потока" (out) */
void diff_flux_ncons_reduced( struct ParametersCommon *paramsc, double *v_ncons, Phase phase, double *flux ) {

    double v_cons[M_REDUCTION];     /* сокращенный однофазный вектор "консервативных" переменных */
    
    /* преобразование вектора неконсервативных переменных в вектор "консервативных" */
    convert_noncons_to_cons_reduction( paramsc, v_ncons, phase, v_cons );
    
    flux[R] = v_cons[V];                                        /* поток массы */
    flux[V] = v_ncons[R] * pow( v_ncons[V], 2.0 ) + v_ncons[P]; /* поток импульса */
    flux[P] = v_ncons[V] * ( v_cons[P] + v_ncons[P] );          /* поток энергии */

}

// Расчет неконсервативного вектора правых частей по вектору "консервативных" переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// rhst - неконсервативный вектор правых частей
void rhst_cons( const struct ParametersCommon* paramsc, const double v_cons[M], array1D* rhst, int number_of_scalars ) {

    double v_ncons[M]; // вектор примитивных переменных
        
    convert_cons_to_noncons( paramsc, v_cons, v_ncons, number_of_scalars );

    (*rhst)[B_DISP] = - v_ncons[V_DISP]; // уравнение компактирования
    (*rhst)[R_DISP] = 0.0; // ЗСМ в дисперсной фазе
    (*rhst)[V_DISP] = v_ncons[P_GAS]; // ЗСИ в дисперсной фазе
    (*rhst)[P_DISP] = v_ncons[P_GAS] * v_ncons[V_DISP]; // ЗСЭ в дисперсной фазе
    (*rhst)[R_GAS] = 0.0; // ЗСМ в газовой фазе
    (*rhst)[V_GAS] = - (*rhst)[V_DISP]; // ЗСИ в газовой фазе
    (*rhst)[P_GAS] = - (*rhst)[P_DISP]; // ЗСЭ в газовой фазе

    for ( int i = 0; i < number_of_scalars; i++ )
        (*rhst)[Z0 + i] = 0.0;
}

// Расчет неконсервативного вектора правых частей по вектору примитивных переменных
// paramsc - структура с основными параметрами вычислительного эксперимента
// v_cons - вектор "консервативных" переменных
// rhst - неконсервативный вектор правых частей
void rhst_ncons( const struct ParametersCommon* paramsс, const double v_ncons[M], array1D* rhst, int number_of_scalars  ) {

    double p_i = calc_p_i( paramsс, v_ncons ); // усредненное межфазное давление
    double u_i = calc_u_i( paramsс, v_ncons ); // усредненная межфазная скорость
    
    (*rhst)[B_DISP] = - u_i; // уравнение компактирования
    (*rhst)[R_DISP] = 0.0; // ЗСМ в дисперсной фазе
    (*rhst)[V_DISP] = p_i; // ЗСИ в дисперсной фазе
    (*rhst)[P_DISP] = p_i * u_i; // ЗСЭ в дисперсной фазе
    (*rhst)[R_GAS] = 0.0; // ЗСМ в газовой фазе
    (*rhst)[V_GAS] = - p_i; // ЗСИ в газовой фазе
    (*rhst)[P_GAS] = - (*rhst)[P_DISP]; // ЗСЭ в газовой фазе

    for ( int i = 0; i < number_of_scalars; i++ )
        (*rhst)[Z0 + i] = 0.0;

}

// Отладочная печать компонент вектора в ячейке
// debug_info - структура с отладочной информацией
void vector_debug_print( const struct DebugInfo *debug_info ) {

    printf( "Vector of primitive variables in cell:\n" );
    printf( "b1 = %.3e\n", debug_info->current_cell_vncons[B_DISP] );
    printf( "r1 = %.3e\n", debug_info->current_cell_vncons[R_DISP] );
    printf( "v1 = %.3e\n", debug_info->current_cell_vncons[V_DISP] );
    printf( "p1 = %.3e\n", debug_info->current_cell_vncons[P_DISP] );
    printf( "r2 = %.3e\n", debug_info->current_cell_vncons[R_GAS] );
    printf( "v2 = %.3e\n", debug_info->current_cell_vncons[V_GAS] );
    printf( "p2 = %.3e\n", debug_info->current_cell_vncons[P_GAS] );

}

// Печать полной информации о проблемном месте в случае аварийной остановки в методе Годунова
// debug_info - структура с отладочной информацией (in)
// index - индекс проблемного элемента в массиве давлений (in)
void debug_print( struct DebugInfo *debug_info, int index ) {

    // проблема с давлением в какой фазе?
    switch ( index ) {
        case GAS_LEFT:
        case GAS_RIGHT:
            // компоненты, отвечающие давлениям в газовой фазе
            printf( "for the gas phase\n" );
            break;
        case DISP_LEFT:
        case DISP_RIGHT:
            // компоненты, отвечающие давлениям в дисперсной фазе
            printf( "for the dispersed phase\n" );
            break;
        }

    printf( "\nDebug report: current cell %d (x = %f), neighbour cell %d\n", debug_info->current_cell, debug_info->current_cell_x,
        debug_info->neighbour_cell );

    printf( "\nCurrent cell primitive vector:\n\n" );
    printf( "b1 = %.3e\n", debug_info->current_cell_vncons[B_DISP] );
    printf( "r1 = %.3e\n", debug_info->current_cell_vncons[R_DISP] );
    printf( "v1 = %.3e\n", debug_info->current_cell_vncons[V_DISP] );
    printf( "p1 = %.3e\n", debug_info->current_cell_vncons[P_DISP] );
    printf( "r2 = %.3e\n", debug_info->current_cell_vncons[R_GAS] );
    printf( "v2 = %.3e\n", debug_info->current_cell_vncons[V_GAS] );
    printf( "p2 = %.3e\n", debug_info->current_cell_vncons[P_GAS] );

    printf( "\nNeighbour cell primitive vector:\n\n" );
    printf( "b1 = %.3e\n", debug_info->neighbour_cell_vncons[B_DISP] );
    printf( "r1 = %.3e\n", debug_info->neighbour_cell_vncons[R_DISP] );
    printf( "v1 = %.3e\n", debug_info->neighbour_cell_vncons[V_DISP] );
    printf( "p1 = %.3e\n", debug_info->neighbour_cell_vncons[P_DISP] );
    printf( "r2 = %.3e\n", debug_info->neighbour_cell_vncons[R_GAS] );
    printf( "v2 = %.3e\n", debug_info->neighbour_cell_vncons[V_GAS] );
    printf( "p2 = %.3e\n\n", debug_info->neighbour_cell_vncons[P_GAS] );

}

// Приведение параметров задачи к безразмерному виду
// paramsc - структура с основными параметрами вычислительного эксперимента, часть полей проходит процедуру приведения к безразмерному виду (in/out)
// params1d - структура с параметрами одномерной задачи, часть полей проходит процедуру приведения к безразмерному виду (in/out)
void dimensionalization( struct ParametersCommon *paramsc, struct Parameters1d *params1d, struct Parameters2d *params2d ) {
	
    double density_scale = paramsc->mass_scale / pow( paramsc->length_scale, 3.0 ); // характерный масштаб плотности
    double velocity_scale = paramsc->length_scale / paramsc->time_scale; // характерный масштаб скорости
    double pressure_scale = paramsc->mass_scale / paramsc->length_scale / pow( paramsc->time_scale, 2.0 ); // характерный масштаб давления

    double energy_scale = pressure_scale / density_scale;

    // шаг по времени
    if ( paramsc->constant_time_step )
        paramsc->dt /= paramsc->time_scale;

    // УРС
    paramsc->p01 /= pressure_scale; // константа в уравнении состояния дисперсной фазы
    paramsc->p02 /= pressure_scale; // константа в уравнении состояния газовой фазы
    // коэффициент динамической вязкости в газе
    paramsc->mu2 /= ( paramsc->mass_scale / ( paramsc->length_scale * paramsc->time_scale ) );

    // вывод результатов
    paramsc->stop_time /= paramsc->time_scale;

    // фоновые параметры для особого случая отсутствия дисперсной фазы по одну из сторон от разрыва
    paramsc->background_density /= density_scale; // плотность дисперсной фазы
    paramsc->background_velocity /= velocity_scale; // скорость дисперсной фазы
    paramsc->background_pressure /= pressure_scale; // давление дисперсной фазы

    // "физика"
    paramsc->particle_diameter /= paramsc->length_scale; // диаметр частицы дисперсной фазы
    paramsc->compaction_viscosity /= (pressure_scale * paramsc->time_scale);
    paramsc->interface_drag_coef /= (density_scale / paramsc->time_scale);
	
    paramsc->tau_parameter /= (energy_scale / paramsc->mass_scale);

    if (paramsc->program_name == ONED2PHC || paramsc->program_name == ONED3PHC){

        // границы расчетной области
        params1d->left_boundary_x /= paramsc->length_scale;
        params1d->right_boundary_x /= paramsc->length_scale;

        // начальные условия
        for ( int i_block = 0; i_block < params1d->ic_blocks_number; i_block++ ) {
            params1d->block_values[i_block][R_DISP] /= density_scale; // плотность дисперсной фазы
            params1d->block_values[i_block][V_DISP] /= velocity_scale; // скорость дисперсной фазы
            params1d->block_values[i_block][P_DISP] /= pressure_scale; // давление дисперсной фазы
            params1d->block_values[i_block][R_GAS] /= density_scale; // плотность газовой фазы
            params1d->block_values[i_block][V_GAS] /= velocity_scale; // скорость газовой фазы
            params1d->block_values[i_block][P_GAS] /= pressure_scale; // давление газовой фазы
        }
    
        // параметры вдува
        if ( params1d->left_bc == INFLOW || params1d->right_bc == INFLOW || params1d->left_bc == COMPLEX_INFLOW || params1d->right_bc == COMPLEX_INFLOW ) {
            params1d->inflow_values[R_DISP] /= density_scale; // плотность дисперсной фазы
            params1d->inflow_values[V_DISP] /= velocity_scale; // скорость дисперсной фазы
            params1d->inflow_values[P_DISP] /= pressure_scale; // давление дисперсной фазы
            params1d->inflow_values[R_GAS] /= density_scale; // плотность газовой фазы
            params1d->inflow_values[V_GAS] /= velocity_scale; // скорость газовой фазы
            params1d->inflow_values[P_GAS] /= pressure_scale; // давление газовой фазы
        }

        // параметры сложного вдува
        if ( params1d->left_bc == COMPLEX_INFLOW || params1d->right_bc == COMPLEX_INFLOW ) {
            params1d->change_inflow_time /= paramsc->time_scale;
            params1d->changed_inflow_values[R_DISP] /= density_scale; // плотность дисперсной фазы
            params1d->changed_inflow_values[V_DISP] /= velocity_scale; // скорость дисперсной фазы
            params1d->changed_inflow_values[P_DISP] /= pressure_scale; // давление дисперсной фазы
            params1d->changed_inflow_values[R_GAS] /= density_scale; // плотность газовой фазы
            params1d->changed_inflow_values[V_GAS] /= velocity_scale; // скорость газовой фазы
            params1d->changed_inflow_values[P_GAS] /= pressure_scale; // давление газовой фазы
        }
    }
    else if (paramsc->program_name == TWOD2PHC){

        double velocity_scale1 = paramsc->length_scale / paramsc->time_scale;    /* характерный масштаб скорости */
        double velocity_scale2 = paramsc->length_scale / paramsc->time_scale;   // (S) Ввёл 2-ую функцию для скорости просто для различия (которое возможно и не нужно)

        /* границы расчетной области */
        params2d->left_boundary_x /= paramsc->length_scale;
        params2d->right_boundary_x /= paramsc->length_scale;

        /* начальные условия */
        for ( int i_block = 0; i_block < params2d->x_blocks_number; i_block++ ) 
        {
            for ( int j_block = 0; j_block < params2d->y_blocks_number; j_block++ ) 
            {
                params2d->block_values[i_block][j_block][R_DISP_2D] /= density_scale;     /* плотность дисперсной фазы */
                params2d->block_values[i_block][j_block][V_DISP_2D] /= velocity_scale1;    /* скорость дисперсной фазы */
                params2d->block_values[i_block][j_block][U_DISP_2D] /= velocity_scale2;    /* скорость дисперсной фазы */
                params2d->block_values[i_block][j_block][P_DISP_2D] /= pressure_scale;    /* давление дисперсной фазы */
                params2d->block_values[i_block][j_block][R_GAS_2D] /= density_scale;     /* плотность газовой фазы */
                params2d->block_values[i_block][j_block][V_GAS_2D] /= velocity_scale1;    /* скорость газовой фазы */
                params2d->block_values[i_block][j_block][U_GAS_2D] /= velocity_scale2;    /* скорость газовой фазы */
                params2d->block_values[i_block][j_block][P_GAS_2D] /= pressure_scale;    /* давление газовой фазы */
            }
        }

        if ( params2d->left_bc == INFLOW ) 
        {
            params2d->inflow_values_left_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_left_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->right_bc == INFLOW ) 
        {
            params2d->inflow_values_right_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_right_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->down_bc == INFLOW ) 
        {
            params2d->inflow_values_down_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_down_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_down_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_down_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_down_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_down_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_down_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_down_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->up_bc == INFLOW ) 
        {
            params2d->inflow_values_up_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_up_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_up_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_up_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_up_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_up_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_up_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_up_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->left_bc == COMPLEX_INFLOW ) 
        {
            params2d->inflow_values_left_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_left_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_left_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_left_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_left_bc[P_GAS_2D] /= pressure_scale;    
        }
        if ( params2d->right_bc == COMPLEX_INFLOW ) 
        {
            params2d->inflow_values_right_bc[R_DISP_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_DISP_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_DISP_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_DISP_2D] /= pressure_scale;    
            params2d->inflow_values_right_bc[R_GAS_2D] /= density_scale;     
            params2d->inflow_values_right_bc[V_GAS_2D] /= velocity_scale1;   
            params2d->inflow_values_right_bc[U_GAS_2D] /= velocity_scale2;    
            params2d->inflow_values_right_bc[P_GAS_2D] /= pressure_scale;    
        }

    }
}

void current_block_number(struct ParametersCommon *params, struct Parameters1d *params1d, struct Parameters2d *params2d, int step_number_x, int step_number_y, int *number_of_block){

    number_of_block[0] = 0;

    if (params->program_name == ONED2PHC || params->program_name == ONED3PHC){
        for (int i = 0; i < params1d->ic_blocks_number; i++){
	    if (step_number_x <= params1d->cell_end[i] && step_number_x >= params1d->cell_begin[i])
	        number_of_block[0] = i;
        }
    }
    else if (params->program_name == TWOD2PHC){
        number_of_block[1] = 0;
        for (int i = 0; i < params2d->x_blocks_number; i++){
            for (int j = 0; j < params2d->y_blocks_number; j++){
                if (step_number_x <= params2d->x_end[i][j] && step_number_x >= params2d->x_begin[i][j] && step_number_y <= params2d->y_end[i][j] && step_number_y >= params2d->y_begin[i][j]){
                    number_of_block[0] = i;
                    number_of_block[1] = j;
                }
            }
        }
    }
    else{
        printf("Wrong program name is used in function current_block_number_2d\n");
        exit(EXIT_FAILURE);
    }
   
}

