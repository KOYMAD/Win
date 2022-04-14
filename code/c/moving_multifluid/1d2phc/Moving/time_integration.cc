#include "time_integration.h"

// Определяет, продолжать расчет или достигнуто требуемое время/количество шагов
// params - структура с параметрами вычислительного эксперимента
// time_iter_num - текущее количество шагов по времени
// curr_time - текущий момент времени
// Возвращает true, если нужно продолжать расчет; false - иначе
void is_time_evolution( const struct Parameters* params, struct TimeMoment time_mom, bool *condition_to_enter_the_cycle )
{
    switch( params->exit_time_cycle )
    {
        case ITERATIONS:
            // критерий останова - заданное количество шагов по времени
            if ( time_mom.steps_num < params->time_cycle_iterations )
                *condition_to_enter_the_cycle = true;
            else
                *condition_to_enter_the_cycle = false;
            break;
        case FINAL_TIME:
            // критерий останова - заданный момент времени
            if ( time_mom.curr_t < params->t_fin )
                *condition_to_enter_the_cycle = true;
            else
                *condition_to_enter_the_cycle = false;
            break;
        default:
            printf( "\nis_time_evolution -> wrong exit_time_cycle value\n\n" );
            system ( "Pause" );
    }
}

//Функция выбирает шаг по времени и возвращает его
//params - структура с параметрами вычислительного эксперимента (in)
//*conservative - массив консервативных переменных (out)
// status - массив статусов расчетных ячеек
void calc_time_step ( struct Parameters *params, int *status, struct Conservative_vector *conservative, double body_velocity, double *time_step )
{
    double char_time; // время за которое возмущение не успеет распространиться за край ячейки
    switch ( params->time_step_method )
    {
        case DYNAMIC_TIME_STEP:
        {
            double delta = ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary ) / params->number_of_cells;
            double c; // скорость звука в текущей ячейке
            for ( int i = 0 ; i < params->number_of_cells; i++ )
            {
                struct Primitive_vector primitive;
                if ( INNER == status[i] || BOUNDARY == status[i] )
                {
                    calc_primitive_variables( params, conservative[i], &primitive );
                    calc_sound_velocity ( params, primitive, &c ); // расчёт скорости звука в текущей ячейке
                    char_time = ( delta / ( abs ( primitive.velosity ) + c ) );	
                    if ( char_time < *time_step ) *time_step = char_time;
                }
            }
            if ( (*time_step) * params->cfl > delta / abs ( body_velocity ) )
            {
	        printf( "\ncalc_time_step -> too big time step or body velocity\n" );
	        system ( "Pause" );
	    }
            break;
        }
        case CONSTANT_TIME_STEP:
            *time_step = params->dt;
            break;
        default:
            system ( "Pause" );
    }
    *time_step *= params->cfl;
}

void time_integration ( struct Parameters *params, struct Conservative_vector *conservative,
                       double time_step, int *status, double body_velocity, double *coordinate_of_left_boundary_of_body,
                      double *coordinate_of_right_boundary_of_body )
{
    struct Flux_vector *array_of_fluxes; // массив потоков через рёбра
    get_memory_for_1D_flux_vector_array ( params->number_of_cells + 1, &array_of_fluxes );

    // консервативные переменные слева и справа от границ расчётной области
    struct Conservative_vector behind_left_boundary_conservative_vector;
    struct Conservative_vector behind_right_boundary_conservative_vector;
    // учёт граничных условий - расчёт консервативных переменных слева и справа от расчётной области
    boundary_conditions ( params, conservative, &behind_left_boundary_conservative_vector,
        &behind_right_boundary_conservative_vector );
    // расчёт потоков на левой и правой границах расчётной области
    godunov_flux ( params, behind_left_boundary_conservative_vector, conservative[0], &(array_of_fluxes[0]) );
    godunov_flux ( params, conservative[params->number_of_cells - 1], behind_right_boundary_conservative_vector,
        &(array_of_fluxes[params->number_of_cells]) );
    // вектор консервативных переменных на контактном разрыве на левой границе подвижного тела
    struct Conservative_vector contact_discontinuity_conservative_on_left_border;
    // вектор консервативных переменных на контактном разрыве на правой границе подвижного тела
    struct Conservative_vector contact_discontinuity_conservative_on_right_border;
    int index_of_left_ghost_cell, index_of_right_ghost_cell; // индексы левой и правой виртуальных ячеек

    // расчёт потоков
    for ( int i = 0 ; i < params->number_of_cells - 1 ; i++ ) // цикл по правым рёбрам ячеек вплоть до предпоследней
    {
        if ( INNER == status[i] ) // если ячейка внутренняя
            // аргументы функции Римана - параметры в левой и правой ячейках соответственно
            godunov_flux ( params, conservative[i], conservative[i+1], &array_of_fluxes[i+1] );
        else if ( BOUNDARY == status[i] ) // если ячейка граничная
        {
            if ( INNER == status[i+1] ) // если ячейка справа внутренняя
                // аргументы функции Римана - параметры в левой и правой ячейках соответственно
                godunov_flux ( params, conservative[i], conservative[i+1], &array_of_fluxes[i+1] );
            else if ( GHOST == status[i+1] ) // если ячейка справа виртуальная
            {
                index_of_left_ghost_cell = i + 1;
                // в качестве параметров справа выбирается величина, зависящая от скорости границы тела
	        calc_flux_through_edge_of_ghost_cell ( params, conservative[i], body_velocity, FROM_LEFT_TO_RIGHT,
                    &array_of_fluxes[i+1], &contact_discontinuity_conservative_on_left_border );
            }
        }
        else if ( GHOST == status[i] ) // если ячейка виртуальная
        {
            if ( BOUNDARY == status[i+1] ) // если ячейка справа граничная
            {
                index_of_right_ghost_cell = i;
                // в качестве параметров слева выбирается величина, зависящая от скорости границы тела
                calc_flux_through_edge_of_ghost_cell ( params, conservative[i+1], body_velocity, FROM_RIGHT_TO_LEFT,
                    &array_of_fluxes[i+1], &contact_discontinuity_conservative_on_right_border );
            }
            else continue;
        }
	else if ( OUTER == status[i] ) {} // если ячейка внешняя, потоки на её рёбрах не рассчитываются
	else
	{
	    printf( "\ntime_integration -> wrong status value\n\n" );
	    system ( "Pause" );
	}
    }

    // нахождение консервативных переменных на новом шаге по времени
    for ( int i = 0 ; i < params->number_of_cells ; i++ ) // цикл по ячейкам
    {
        if ( INNER == status[i] || BOUNDARY == status[i] ) // если ячейка внутренняя или граничная
            calc_conservative_on_next_time_moment ( params, &(conservative[i]), array_of_fluxes[i],
            array_of_fluxes[i+1], time_step );
        else continue;
    }

    // передвижение тела и определение новых координат его левой границы
    moving_the_body ( body_velocity, time_step, coordinate_of_left_boundary_of_body,
        coordinate_of_right_boundary_of_body );

    // пересчёт статусов ячеек и инициализация ячеек с новым статусом
    calc_status ( params, *coordinate_of_left_boundary_of_body, *coordinate_of_right_boundary_of_body,
        contact_discontinuity_conservative_on_left_border, contact_discontinuity_conservative_on_right_border,
        status, conservative, index_of_left_ghost_cell, index_of_right_ghost_cell );
}

// расчёт потока через ребро, соединяющее граничную ячейку и виртуальную, содержащую в себе границу подвижного тела
void calc_flux_through_edge_of_ghost_cell ( struct Parameters *params, struct Conservative_vector conservative,
                                           double body_velocity, enum Direction dir, struct Flux_vector *flux, 
                                           struct Conservative_vector *contact_discontinuity_conservative )
{
    struct Primitive_vector Riemann_primitive_solution; // вектор примитивных переменных решения задачи Римана
    struct Conservative_vector ghost_conservative; // вектор консервативных переменных в виртуальной ячейке
    struct Primitive_vector ghost_primitive; // вектор примитивных переменных в виртуальной ячейке
    struct Primitive_vector current_primitive; // вектор примитивных переменных в текущей граничной ячейке
    calc_primitive_variables ( params, conservative, &current_primitive );
    ghost_primitive = current_primitive;
    ghost_primitive.velosity = 2 * body_velocity - current_primitive.velosity;
    calc_conservative_variables ( params, ghost_primitive, &ghost_conservative );
    struct Primitive_vector contact_discontinuity_primitive; // примитивные переменные на контактном разрыве
    // предварительное решение задачи Римана и определение давления на контактном разрыве
    if ( FROM_LEFT_TO_RIGHT == dir )
    {
        solve_Riemann_problem ( params, conservative, ghost_conservative, &Riemann_primitive_solution,
            &contact_discontinuity_primitive );
    }
    else if ( FROM_RIGHT_TO_LEFT == dir )
    {
        solve_Riemann_problem ( params, ghost_conservative, conservative, &Riemann_primitive_solution,
            &contact_discontinuity_primitive );
    }
    // расчёт консервативных переменных на контактном разрыве
    calc_conservative_variables ( params, contact_discontinuity_primitive, contact_discontinuity_conservative );
    //calc_conservative_variables ( params, Riemann_primitive_solution, Riemann_conservative_solution );
    // скорость звука на контактном разрыве
    double contact_discontinuity_sound_velocity;
    calc_sound_velocity ( params, contact_discontinuity_primitive, &contact_discontinuity_sound_velocity );
    // если тело движется влево быстрее скорости звука на контактном разрыве
    if ( body_velocity < -contact_discontinuity_sound_velocity )
    {
        if ( FROM_LEFT_TO_RIGHT == dir ) diff_flux_ncons( params, Riemann_primitive_solution, flux );
        else if ( FROM_RIGHT_TO_LEFT == dir ) diff_flux_ncons( params, current_primitive, flux );
    }
    // если тело движется медленнее скорости звука на контактном разрыве
    else if ( body_velocity < contact_discontinuity_sound_velocity )
    {
        if ( FROM_LEFT_TO_RIGHT == dir ) godunov_flux ( params, conservative, *contact_discontinuity_conservative, flux );
        else if ( FROM_RIGHT_TO_LEFT == dir ) godunov_flux ( params, *contact_discontinuity_conservative, conservative, flux );
        else
        {
            printf( "\ncalc_flux_through_edge_of_ghost_cell -> wrong Direction value\n\n" );
            system ( "Pause" );
        }
    }
    // если тело движется вправо быстрее скорости звука на контактном разрыве
    else
    {
        if ( FROM_LEFT_TO_RIGHT == dir ) diff_flux_ncons( params, current_primitive, flux );
        else if ( FROM_RIGHT_TO_LEFT == dir ) diff_flux_ncons( params, Riemann_primitive_solution, flux );
    }
}

void calc_conservative_on_next_time_moment ( struct Parameters *params, struct Conservative_vector *conservative,
                                    struct Flux_vector left_flux, struct Flux_vector right_flux, double time_step )
{
    double grid_step = ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary ) /
        params->number_of_cells; // шаг сетки

    conservative->specific_mass -= ( right_flux.mass_flux - left_flux.mass_flux ) * time_step / grid_step;
    conservative->specific_momentum -= ( right_flux.momentum_flux - left_flux.momentum_flux ) * time_step / grid_step;
    conservative->specific_energy -= ( right_flux.energy_flux - left_flux.energy_flux ) * time_step / grid_step;
}

// int i (in) - индекс виртуальной ячейки (содержащей границу тела)
//void calc_status ( Parameters *params, double coordinate_of_left_boundary_of_body, double coordinate_of_right_boundary_of_body,
//    struct Conservative_vector contact_discontinuity_conservative_on_left_border,
//    struct Conservative_vector contact_discontinuity_conservative_on_right_border,
//    int *status, struct Conservative_vector *conservative, int index_of_left_ghost_cell, int index_of_right_ghost_cell )
//{
//    struct Conservative_vector null_conservative = {0.0, 0.0, 0.0};
//    double grid_step = ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary )
//        / params->number_of_cells;
//
//    // левая граница тела сдвинулась направо в ячейку i + 1
//    if ( coordinate_of_left_boundary_of_body >= ( index_of_left_ghost_cell + 1 ) * grid_step )
//    {
//        status[index_of_left_ghost_cell + 1] = GHOST;
//        status[index_of_left_ghost_cell] = BOUNDARY;
//        conservative[index_of_left_ghost_cell] = contact_discontinuity_conservative_on_left_border;
//        status[index_of_left_ghost_cell - 1] = INNER;
//    }
//    else if ( coordinate_of_left_boundary_of_body >= index_of_left_ghost_cell * grid_step ) // осталась в той же ячейке
//    {}
//    else // сдвинулась налево в ячейку i - 1
//    {
//        status[index_of_left_ghost_cell - 2] = BOUNDARY;
//        status[index_of_left_ghost_cell - 1] = GHOST;
//        status[index_of_left_ghost_cell] = OUTER;
//        conservative[index_of_left_ghost_cell] = null_conservative;
//    }
//
//    // правая граница тела сдвинулась направо в ячейку i + 1
//    if ( coordinate_of_right_boundary_of_body >= ( index_of_right_ghost_cell + 1 ) * grid_step ) // сдвинулась направо в ячейку i + 1
//    {
//        status[index_of_right_ghost_cell + 2] = BOUNDARY;
//        status[index_of_right_ghost_cell + 1] = GHOST;
//        status[index_of_right_ghost_cell] = OUTER;
//        conservative[index_of_right_ghost_cell] = null_conservative;
//    }
//    else if ( coordinate_of_right_boundary_of_body >= index_of_right_ghost_cell * grid_step ) // осталась в той же ячейке
//    {}
//    else // сдвинулась налево в ячейку i - 1
//    {
//        status[index_of_right_ghost_cell - 1] = GHOST;
//        status[index_of_right_ghost_cell] = BOUNDARY;
//        conservative[index_of_right_ghost_cell] = contact_discontinuity_conservative_on_right_border;
//        status[index_of_right_ghost_cell + 1] = INNER;
//    }
//}