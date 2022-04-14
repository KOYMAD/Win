#ifndef __STRUCTURES_H_
#define __STRUCTURES_H_

#include "constants.h"

struct Primitive_vector // вектор примитивных переменных
{
    double density; // плотность
    double velosity; // скорость
    double pressure; // давление
};

struct Conservative_vector // вектор консервативных переменных
{
    double specific_mass; // удельная масса
    double specific_momentum; // удельный импульс
    double specific_energy; // удельная энергия
};

struct Flux_vector // вектор потоков
{
    double mass_flux; // поток массы
    double momentum_flux; // поток импульса
    double energy_flux; // поток энергии
};

// блок для задания начальных данных
struct Block
{
    int number_of_left_cell; // левая ячейка блока
    int number_of_right_cell; // правая ячейка блока
    struct Primitive_vector initial_primitive_vector; // начальные данные в блоке
};

struct Parameters // параметры расчёта
{
    // сетка
    int number_of_cells; // число ячеек
    double coordinate_of_left_boundary; // координата левой границы
    double coordinate_of_right_boundary; // координата правой границы

    // координаты подвижного тела
    double initial_coordinate_of_left_boundary_of_body; // координата левой границы тела
    double initial_coordinate_of_right_boundary_of_body; // координата правой границы тела

    // движение тела
    double body_velosity; // скорость тела
    // площадь поперечного сечения, отесённая к массе тела
    double body_cross_section_devided_to_mass; 

    // условие выхода из цикла по времени: 
    int exit_time_cycle; // 1 - по прошествии time_cycle_iterations итераций, 
                         // 2 - по достижении величиной текущего момента времени значения t_fin
    int time_cycle_iterations; // число заходов в цикл
    double t_fin; // конечный момент времени

    // метод расчёта шага по времени:
    int time_step_method; // 1 - динамический, из условия устойчивости Куранта-Фридрихса-Леви
                          // 2 - статический, dt
    double cfl; // число Куранта-Фридрихса-Леви
    double dt; // шаг по времени
    
    // параметры обезразмеривания
    double length_diml;
    double time_diml;
    double mass_diml;

    // запись выходных данных
    int number_of_output_files; // число записываемых за расчёт файлов распределений параметров газа при exit_time_cycle == 2
    int step_of_output_files; // шаг записи в файл распределений параметров газа при exit_time_cycle == 1

    // численный метод
    int flux_calculation_method;
    int max_iter_num;
    double p_max_ratio;
    int approximation_order; // порядок аппроксимации разностной схемы по пространству
    
    // начальные условия
    struct Block blocks[2]; // массив блоков

    // малые константы
    double eps_general;
    double eps_spatial;
    double big;

    // уравнение состояния газа
    double g; // показатель адиабаты

    // типы граничных условий 1 - стенка, 2 - параметры за УВ с числом Маха M
    // 3 - условие экстраполяции нулевого порядка
    int type_of_left_border;
    int type_of_right_border;
    double Mach_number; // число Маха УВ

    // флаг отладочного режима
    bool is_debug_regime; // true - включен режим отладочной записи файлов
    bool is_сontinue; // true - включена возможность старта расчёта не с начального момента времени
    bool is_calculate_drag_force; // true - включен режим расчёта и вывода сил, действующих на границы
    bool is_initiate_parameters_in_outer_cells; // true - внешние ячейки инициализируются
    // параметрами в соответствии в блоками начальных данных, иначе числом -big
};

// структура, определяющая текущий момент расчетамомент времени и номер шага
struct TimeMoment
{
    int steps_num;  // текущее количество шагов
    double curr_t;  // текущее время
};

#endif /* __STRUCTURES_H_ */