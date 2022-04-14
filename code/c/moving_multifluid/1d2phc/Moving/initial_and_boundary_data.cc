#include "initial_and_boundary_data.h"

//    Функция заполняет массив начальных данных
//    *conservative - массив значений (out)
void initiate_data ( struct Parameters *params, struct Conservative_vector *conservative,
                    int *status, struct TimeMoment *time_mom, int *zones_count )
{
    FILE *restart_info_file = NULL; // дескриптор файла с информацией о последнем доступном файле с промежуточными результатами
    int i_block;  // индекс блоков для задания начальных условий
    struct Conservative_vector null_conservative = {0.0, 0.0, 0.0};

    // попытка считать файл с информацией о рестарте, определяется значение params->is_continue
    //try_to_restart( params, time_mom, restart_info_file );
    if ( params->is_сontinue ) // запускаемый расчет является продолжением предыдущего
    {
        //read_last_time_moment( params, time_mom ); // считывание последнего момента времени
        *zones_count = (int)(params->number_of_output_files * time_mom->curr_t / params->t_fin);
        //read_solution( params, conservative, time_mom );
    }
    else
    {
        for ( i_block = 0; i_block < 2; i_block++ ) // 2 блока начальных данных - слева и справа от поршня
        {
            for ( int i = params->blocks[i_block].number_of_left_cell ; 
                i <= params->blocks[i_block].number_of_right_cell ; i++ )
                if ( INNER == status[i] || BOUNDARY == status[i] )
                    calc_conservative_variables( params, params->blocks[i_block].initial_primitive_vector,
                    &(conservative[i]) );
                else
                    conservative[i] = null_conservative;
        }
    }
}

void initiate_status( struct Parameters params, double left, double right, int *status )
{
    double grid_step = ( params.coordinate_of_right_boundary - params.coordinate_of_left_boundary ) /
        params.number_of_cells; // шаг сетки
    for ( int i = 0 ; i < params.number_of_cells ; i++ )
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
                    if ( right >= ( i + 1 ) * grid_step ) // если правая граница тела правее правой границы ячейки
                        status[i] = GHOST;
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
                    else // если правая граница тела левее правой границы ячейки
                        status[i] = GHOST;
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
    for ( int i = 0 ; i < params.number_of_cells ; i++ )
    {
        if ( GHOST == status[i] )
        {
            if ( i != 0 && INNER == status[i - 1] && OUTER == status[i + 1] )
                status[i - 1] = BOUNDARY;
            else if ( i != params.number_of_cells - 1 && INNER == status[i + 1] 
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

void boundary_conditions ( struct Parameters *params, struct Conservative_vector *conservative, 
struct Conservative_vector *behind_left_boundary_conservative_vector, struct Conservative_vector *behind_right_boundary_conservative_vector )
{
	if ( WALL == params->type_of_left_border )
	{
		*behind_left_boundary_conservative_vector = conservative[0];
		behind_left_boundary_conservative_vector->specific_momentum *= -1;
	}
	else if ( PARAMETERS_BEHIND_SHOCK_WAVE == params->type_of_left_border )
	{}
	else if ( ZERO_ORDER == params->type_of_left_border )
		*behind_left_boundary_conservative_vector = conservative[0];
	else
	{
        printf( "\nincorrect boundary type\n" );
        system ( "Pause" );
    }

	if ( WALL == params->type_of_right_border )
	{
		*behind_right_boundary_conservative_vector = conservative[params->number_of_cells - 1];
		behind_right_boundary_conservative_vector->specific_momentum *= -1;
	}
	else if ( PARAMETERS_BEHIND_SHOCK_WAVE == params->type_of_right_border )
	{}
	else if ( ZERO_ORDER == params->type_of_right_border )
		*behind_right_boundary_conservative_vector = conservative[params->number_of_cells - 1];
	else
	{
        printf( "\nincorrect boundary type\n" );
        system ( "Pause" );
    }
}