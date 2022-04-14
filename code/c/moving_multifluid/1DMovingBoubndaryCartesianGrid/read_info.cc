#include "read_info.h"

// Считывание файла с параметрами задачи, заполнение полей соответствующей структуры и проверка корректности задания параметров
// params - структура с параметрами вычислительного эксперимента
void read_parameters( struct Parameters *params )
{
    FILE *parameters_file;
    char string[MAX_STRING_SIZE]; // для считывания строковой информации из файла
    int i_block; // для считывания блоков с начальными условиями
        
    // все параметры задачи находятся в файле parameters2D.txt
    if ( ( fopen_s( &parameters_file, "parameters.dat", "rt" ) ) != 0 )
    {
        printf( "\nread_parameters -> can't open file parameters.dat for reading\n\n" );
        system("Pause");
    }

    // считывание заголовка файла
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );

    // расчетная сетка
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->number_of_cells) ); // количество ячеек
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->coordinate_of_left_boundary) ); // координата левой границы расчетной области
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->coordinate_of_right_boundary) ); // координата правой границы расчетной области

    // подвижное тело
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->initial_coordinate_of_left_boundary_of_body) ); // координата левой границы тела
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->initial_coordinate_of_right_boundary_of_body) ); // координата правой границы тела
    if ( params->initial_coordinate_of_right_boundary_of_body - params->initial_coordinate_of_left_boundary_of_body < 
        ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary ) / params->number_of_cells )
    {
        printf( "\nread_parameters -> too small body\n" );
        system("Pause");
    }

    // движение тела
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->body_velosity) ); // скорость тела
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->body_cross_section_devided_to_mass) );

    // шаг по времени 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );              // заголовок раздела 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->exit_time_cycle) );    // условие выхода из цикла по времени
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->time_cycle_iterations) );    // число заходов в цикл
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->t_fin) );    // конечный момент времени
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->time_step_method) );   // метод расчёта шага по вемени
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->cfl) );   // число cfl постоянный шаг по времени
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->dt) );  // величина шага

    // параметры обезразмеривания
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );              // заголовок раздела 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->length_diml) ); 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->time_diml) ); 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->mass_diml) ); 

    // запись выходных данных
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );              // заголовок раздела 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->number_of_output_files) );
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->step_of_output_files) );

    // численный метод 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );                                  // заголовок раздела 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->flux_calculation_method) );  // номер разностной схемы 
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->max_iter_num) );      // максимальное количество итераций 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->p_max_ratio) );  // максимальный перепад по давлению слева
                                                                                           //и справа от разрыва, при котором еще
                                                                                           //используется начальное приближение из
                                                                                           //решения линеаризованной задачи
                                                                                           //в методе Годунова
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->approximation_order) ); // порядок разностной схемы по пространству
    if ( params->approximation_order != 1 && params->approximation_order != 2 )
    {
        printf( "\nread_parameters -> approximation_order should be 1 or 2.%d\n\n", params->approximation_order );
        system("Pause");
    }
	
    // начальные условия 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );                                  // заголовок раздела 
    for ( i_block = 0; i_block < 2; i_block++ ) // цикл по блокам начальных условий - левому и правому
    {
        fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // заголовок - номер блока
        fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).number_of_left_cell ) );
        fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).number_of_right_cell ) );
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).initial_primitive_vector.density ) );
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).initial_primitive_vector.velosity ) );
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &( (params->blocks[i_block]).initial_primitive_vector.pressure ) );
    }

    // малые константы
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );  // заголовок раздела 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->eps_general) );  // регулярное малое число для сравнения вещественных чисел,
                                                                                           //использования в математических функциях, контроля
                                                                                           //сходимости итерационных процессов (если иное не оговорено особо) 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->eps_spatial) ); // малое число для проверки условия совпадения двух точек
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->big) ); // большое число

    // показатель адиабаты 
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );  // заголовок раздела 
    fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->g) );  // заголовок раздела 

    // тип каждого из 4 рёбер изначально прямоугольной границы
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE );  // заголовок раздела - тип границы
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->type_of_left_border ) );
    if ( params->type_of_left_border == PARAMETERS_BEHIND_SHOCK_WAVE )
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->Mach_number ) );
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(params->type_of_right_border ) );
    if ( params->type_of_right_border == PARAMETERS_BEHIND_SHOCK_WAVE )
        fscanf_s( parameters_file, "%s %lf", string, MAX_STRING_SIZE, &(params->Mach_number ) );
    

    // флаг отладочного режима
    fscanf_s( parameters_file, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    int dtmp; // вспомогательная переменная для считывания целых чисел
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // считывание флага отладочного режима
    switch ( dtmp )
    {
        case 0:
            params->is_debug_regime = false;
            break;
        case 1:
            params->is_debug_regime = true;
            break;
        default:
            printf( "\nread_parameters -> debug should be 0 or 1.\n\n" );
            system("Pause");
    }
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &(dtmp) ); // считывание флага запуска расчёта с начала или с некоторого момента
    switch ( dtmp )
    {
        case 0:
            params->is_сontinue = false;
            break;
        case 1:
            params->is_сontinue = true;
            break;
        default:
            printf( "\nread_parameters -> is_сontinue should be 0 or 1.\n\n" );
            system("Pause");
    }
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // считывание флага отладочного режима
    switch ( dtmp )
    {
        case 0:
            params->is_calculate_drag_force = false;
            break;
        case 1:
            params->is_calculate_drag_force = true;
            break;
        default:
            printf( "\nread_parameters -> is_calculate_drag_force should be 0 or 1.\n\n" );
            system("Pause");
    }
    fscanf_s( parameters_file, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // считывание флага инициализации во внешних ячейках
    if ( 0 == dtmp ) params->is_initiate_parameters_in_outer_cells = false;
    else if ( 1 == dtmp ) params->is_initiate_parameters_in_outer_cells = true;
    else
    {
        printf( "\nread_parameters -> is_initiate_parameters_in_outer_cells should be 0 or 1.\n\n" );
        system("Pause");
    }

    fclose( parameters_file );
}

void dimless_parameters( struct Parameters *params )
{
    double length_diml = params->length_diml;
    double time_diml = params->time_diml;
    double dens_diml = params->mass_diml / pow(length_diml, 3);
    double speed_diml = length_diml / time_diml;
    double press_diml = params->mass_diml / ( length_diml * pow(time_diml, 2));

    params->coordinate_of_left_boundary /= length_diml;
    params->coordinate_of_right_boundary /= length_diml;

    params->body_velosity /= speed_diml;
    params->body_cross_section_devided_to_mass /= pow( length_diml, 2 ) / params->mass_diml;

    params->initial_coordinate_of_left_boundary_of_body /= length_diml;
    params->initial_coordinate_of_right_boundary_of_body /= length_diml;

    params->dt /= time_diml;
    params->t_fin /= time_diml;

    for ( int i_block = 0; i_block < 2; i_block++ ) // цикл по блокам начальных условий 
    {
        (params->blocks[i_block]).initial_primitive_vector.density /= dens_diml;
        (params->blocks[i_block]).initial_primitive_vector.velosity /= speed_diml;
        (params->blocks[i_block]).initial_primitive_vector.pressure /= press_diml;
    }
}

// файл с последним моментом времени, на который сохранены данные расчёта, содержит размерное время
void read_last_time_moment( struct Parameters *params, struct TimeMoment *time_mom )
{
    FILE *in;
    if ( ( fopen_s( &in, "last_time_moment.dat", "rt" ) ) != 0 )
    {
        printf( "\nread_last_time_moment -> can't open file read_last_time_moment.dat for reading%d\n\n", params->is_сontinue );
        system ( "Pause" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( in, "%06d", &(time_mom->steps_num) );
    fscanf_s( in, "%lf", &(time_mom->curr_t) );
    time_mom->curr_t /= params->time_diml; // обезразмеривание считанного момента времени
    fclose( in );
}

void read_solution ( struct Parameters *params, struct Conservative_vector *conservative, struct TimeMoment *time_mom )
{
}