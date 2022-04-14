﻿// struct.h
// Общие структуры данных и константы для двухфазных сжимаемых солверов
// вне завимисоти от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 28 августа 2018 г.

#ifndef __STRUCT_H_
#define __STRUCT_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>


using namespace std;

// неизменные константы не для вынесения в конфигурационный файл

#define M 12 // большое число - размерность вектора-решения полной системы 
#define N 5 // количество лагранжевых скаляров с запасом
#define M1D 7 // размерность вектора-решения полной системы 1d задачи
#define M2D 9 // размерность вектора-решения полной системы 2d задачи

#define SMALL_NUMBER 1.e-8  // маленькое число
#define	MAX_STRING_SIZE 1000 // максимальное число символов в строковых переменных
#define INFINITY1 1.e10 // большое число для процедуры выбора шага по времени
#define M_REDUCTION 3 // размерность вектора-решения редуцированной системы для случая отсутствия градиента объемной доли дисперсной фазы
#define K_GENERAL_CASE 4 // размерность вектора для описания условия на контактном разрыве в дисперсной фазе в общем случае наличия дисперсной фазы по обе стороны                          // от первоначального разрыва
#define K_SPECIAL_CASE 3 // размерность вектора для описания условия на контактном разрыве в дисперсной фазе в специальном случае отсутствия дисперсной фазы по одну из сторон от первоначального разрыва
#define MAX_IC_BLOCKS_NUMBER 100 // максимальное количество блоков, на которые можно разбить расчетную область для задания начальных условий
#define MAX_SENSORS_NUM 100 // максимальное количество датчиков для записи величин в определенных точках
#define ADIABATIC_CURVE_PNUM 101 // количество точек на адиабате - множестве состояний, куда можно перевести среду, пустив по ней ударную волну или волну разрежения, для отладки метода Годунова
#define ADAIABATIC_CURVE_P_B 1.e-5 // нижняя граница графика адиабаты на оси ординат (давлений), для отладки метода Годунова
#define ADAIABATIC_CURVE_P_T 0.1 // верхняя граница графика адиабаты на оси ординат (давлений), для отладки метода Годунова
#define RELAX_CASES_NUM 10 // количество веток алгоритма в процедуре расчета релаксации давлений
#define GODUNOV_CASES_NUM 10 // количество веток алгоритма в методе Годунова

// новые типы для работы с многомерными матрицами с использованием контейнера vector

typedef vector<double> array1D; // одномерный динамический вещественный массив
typedef vector<array1D> array2D; // двумерный динамический вещественный массив

// новые типы из перечислений

enum ReturnCodes { SUCCESS, NEGATIVE_PRESSURE, GODUNOV_FAILS }; // коды, возвращаемые расчетными функциями
enum Direction { LEFT, RIGHT }; // идентификатор направления
enum Phase { GAS_PHASE = 1, DISPERSED_PHASE = 2 }; // идентификатор фазы
enum Model { BAER_NUNZIATO = 1, SAUREL_ABGRALL = 2 }; // модель двухфазной среды
enum Program { ONED2PHC = 1, ONED3PHC = 2, TWOD2PHC = 3 }; // запускаемая программа
// ROGUE - формула для засыпки сферических частиц из
// Rogue et al. Experimental and numerical investigation of the shock-induced fluidization of a particles bed // Shock Waves. - 1998. - V. 8. - P. 29 - 45.
// HOUIM - формула для засыпки сферических частиц из
// Houim R.W., Oran E.S. A multiphase model for compressible granular-gaseous flows: formulation and initial tests // J. Fluid Mech. - 2016. - V. 789. - P. 166 - 220.
// TANINO - формула для системы цилиндров из
// Tanino Y., Nepf H.M. Laboratory investigation of mean drag in a random array of rigid, emergent cylinders // J. Hydraul. Eng. - 2008. - V.134, No. 1. - P. 34 - 41.
// SCHWENDEMAN - формула из
//  D. W. Schwendeman , C. W. Wahle & A. K. Kapila. A study of detonation evolution and structure for a model of compressible two-phase reactive flow // Combustion Theory and Modelling. - 2008. - 12:1. - P. 159 - 204.
enum FrictionForce { ROGUE = 1, HOUIM = 2, TANINO = 3, SCHWENDEMAN = 4, HOMENKO = 5 }; // формулы для расчета силы трения

// SEREBRYAKOV - формула из
// Серебряков М.Е. Внутренняя баллистика ствольных систем и пороховых ракет. М.: Оборонгиз, 1962. 
// И. В. Семенов, П. С. Уткин, И. Ф. Ахмедьянов, И. С. Меньшов, Применение многопроцессорной вычислительной техники для решения задач внутренней баллистики, Выч. мет. программирование, 2011, том 12, выпуск 1, 183–193
// SCHWENDEMAN_BURNING - формула из 
// D. W. Schwendeman , C. W. Wahle & A. K. Kapila. A study of detonation evolution and structure for a model of compressible two-phase reactive flow // Combustion Theory and Modelling. - 2008. - 12:1. - P. 159 - 204.
enum BurningLaw { SEREBRYAKOV_BURNING = 1, SCHWENDEMAN_BURNING = 2 }; // формулы для расчета газоприхода 

// SCHWENDEMAN_COMPACTION - формула из 
// D. W. Schwendeman , C. W. Wahle & A. K. Kapila. A study of detonation evolution and structure for a model of compressible two-phase reactive flow // Combustion Theory and Modelling. - 2008. - 12:1. - P. 159 - 204.
// SAUREL_COMPACTION - формула из
// SAUREL, R., FAVRIE, N., PETITPAS, F., LALLEMAND, M., & GAVRILYUK, S. (2010). Modelling dynamic and irreversible powder compaction. Journal of Fluid Mechanics, 664, 348-396. doi:10.1017/S0022112010003794
// Sebastien Bodard, Olivier Jalbaud, Richard Saurel, Yves Burtschell, and Emmanuel Lapebie. (2016). Phase diagram of crushed powders. Physics of Fluids, 28, 123301. doi.org/10.1063/1.4968195
enum CompactionPressureRelaxation { SCHWENDEMAN_COMPACTION = 1, SAUREL_COMPACTION = 2}; // формулы для расчета конфигурационного давления 

// разностные схемы
enum Scheme { CIR1 = 0, // первая модификация схемы Куранта-Изаксона-Рис, только для БН
              CIR2 = 1, // вторая модификация схемы Куранта-Изаксона-Рис, только для БН
              CIR3 = 2, // третья модификация схемы Куранта-Изаксона-Рис, только для БН
              CIR4 = 3, // четвертая модификация схемы Куранта-Изаксона-Рис, только для БН
              GODUNOV = 4, // метод Годунова, только для БН
              HLL = 5, // метод Хартена-Лакса-Ван Лира (HLL), только для СА 
              RUSANOV = 6, // метод Русанова, только для СА 
              HLLC = 7 }; // метод Хартена-Лакса-Ван Лира-контактный (HLLC), только для СА 

enum NumericalMethodPhysics { EXPLICIT_EULER = 1, NEWTON = 2 }; // численный метод для решения физики

// индексы компонент вектора примитивных переменных для одномерной полной системы уравнений Баера-Нунциато
enum Full_vector_index_1D { B_DISP = 0, // объемная доля дисперсной фазы
                         R_DISP = 1, // плотность дисперсной фазы
                         V_DISP = 2, // скорость дисперсной фазы
                         P_DISP = 3, // давление дисперсной фазы
                         R_GAS = 4, // плотность газовой фазы
                         V_GAS = 5, // скорость газовой фазы
                         P_GAS = 6, // давление газовой фазы
                         Z0    = 7, // лагранжев скаляр №0
                         Z1    = 8};// лагранжев скаляр №1

// индексы компонент вектора примитивных переменных для двумерной системы уравнений
enum Full_vector_index_2D { B_DISP_2D = 0, // объемная доля дисперсной фазы
                            R_DISP_2D = 1, // плотность дисперсной фазы
                            V_DISP_2D = 2, // скорость дисперсной фазы по x
                            U_DISP_2D = 3, // скорость дисперсной фазы по y
                            P_DISP_2D = 4, // давление дисперсной фазы
                            R_GAS_2D = 5, // плотность газовой фазы
                            V_GAS_2D = 6, // скорость газовой фазы по x
                            U_GAS_2D = 7, // скорость газовой фазы по y
                            P_GAS_2D = 8 }; // давление газовой фазы

// индексы компонент вектора примитивных переменных для системы уравнений для отдельной фазы
enum Reduced_vector_index_1D { R = 0, // плотность
                            V = 1, // скорость
                            P = 2 }; // давление

/* Вектор "консервативных" переменных для полной системы уравнений
   (0) объемная доля дисперсной фазы ( b1 )
   (1) b1 * плотность дисперсной фазы ( b1 * r1 )
   (2) b1 * r1 * скорость дисперсной фазы ( b1 * r1 * v1 )
   (3) b1 * r1 * полная удельная энергия дисперсной фазы ( b1 * r1 * E1 )
   (4) b2 * r2
   (5) b2 * r2 * v2
   (6) b2 * r2 * E2 */

// тип граничного условия
enum Boundary_condition { WALL = 1, // стенка
                          FREE = 2, // свободное втекание/истечение
                          INFLOW = 3, // втекание с заданными параметрами
                          COMPLEX_INFLOW = 4, // втекание с одними заданными параметрами до заданного момента времени и с другими заданными после него
                          WALL_DISP = 5  }; // стенка для дисперсной фазы

enum Boundary_direction { LEFT_BOUNDARY = 1, // левая граница
                          RIGHT_BOUNDARY = 2, // правая граница
                          UP_BOUNDARY = 3, // верхняя граница
                          DOWN_BOUNDARY = 4 }; // нижняя граница

enum Block_types { NONE = 0, INNER = 1, OUTER = 2 }; // внутренний или внешний блок. 0 - нет ячеек

// возможные соотношения между значениями объемной доли дисперсной фазы слева и справа от разрыва,
// отличающиеся с точки зрения численного метода
enum Disp_phase_cases { BOTH_GRAD, // дисперсная фаза в обеих ячейках, причем перепад значительный
                        LEFT_ONLY_GRAD, // дисперсная фаза только слева от разрыва, причем перепад значительный
                        RIGHT_ONLY_GRAD, // дисперсная фаза только справа от разрыва, причем перепад значительный
                        NO_GRAD }; //перепад между значениями дисперсной фазы слева и справа от разрыва незначительный

enum Pressures_idexes { GAS_LEFT, GAS_RIGHT, DISP_LEFT, DISP_RIGHT }; // индексы в массиве давлений в методе Годунова

// структура с параметрами вычислительного эксперимента, не зависящих от размерности задачи
struct ParametersCommon {
    
    // поля, которые считываются из конфигурационного файла
    int program_name; // название запускаемой программы

    // модель среды
    Model media_model; // для двухфазной среды - модель БН или СА
    bool porous_body; // дисперсная фаза - неподвижный пористый каркас или нет?

    // шаг по времени
    bool constant_time_step; // является ли шаг по времени постоянным или выбирается динамически
    double cfl; // число Куранта-Фридрихса-Леви (имеет смысл, если шаг выбирается динамически)
    double dt; // шаг по времени (имеет смысл, если расчет ведется с постоянным шагом)

    // численный метод
    int numerical_method; // номер разностной схемы, см. enum Scheme
    int max_iter_num; // максимальное количество итераций для поиска давления и скорости на разрыве, только для GODUNOV
    double p_max_ratio; // максимальный перепад по давлению слева и справа от разрыва, при котором
                        // еще используется начальное приближение из решения линеаризованной задачи, только для GODUNOV
    
    bool pressure_relaxation; // использовать релаксацию давления между фазами?
    bool pressure_relaxation_compaction;
    int pressure_relaxation_compaction_formula;
    double tau_parameter;
    double n_parameter;
    double volume_fraction_compaction;

    bool velocity_relaxation; // использовать релаксацию скорости между фазами?
    int approximation_order; // порядок аппроксимации метода - 1 или 2

    // начальные условия
    bool isContinue; // является ли запускаемый расчет продолжанием старого?

    // уравнения состояния фаз
    bool use_real_eos; // использовать реальное УРС для вычисления параметров двучленного УРС?
    double g1, g2; // показатели адибаты в дисперсной и газовой фазах
    double p01, p02; // константы в двучленных уравнения состояния фаз
    double b_virial; // константа в вириальном уравнении состояния газовой фазы
    double mu2; // коэффициент вязкости газа
    int material_number1, material_number2; // номер материала. если нет, то равен 0
    double specific_heat_cv1, specific_heat_cv2; // удельные теплоемости дисперсной и газовой фазы
    
    // режим расчета
    bool is_debug; // требуется ли выполнять проверку корректности матричных операций и решения систем уравнений

    // вывод результатов
    bool is_output_on_time; // является ли критерием вывода результатов достижение заданного момента времени
    double stop_time; // момент времени окончания расчета (если is_output_on_time == true)
    int stop_steps_number; // требуемое количество шагов по времени (если is_output_on_time == false)
    int output_number; // количество файлов с промежуточными результатами
    bool is_tecplot_header; // добавить ли заголовок для визуализации в Tecplot?

    // датчики
    bool are_sensors; // требуется ли использовать датчики в заданных точках?
    int sensors_num; // фактическое количество датчиков

    // "малые" константы в коде
    double eps_general; // регулярное малое число для сравнения вещественных чисел, использования в математических функциях,
                        // контроля сходимости итерационных процессов (если иное не оговорено особо)
    double eps_ludcmp; // малое число в функции ludcmp (math_utilc.cc), оставлено как было в оригинальной работе
    double eps_decouple; // максимально допустимый перепад по объемной доле слева и справа, при котором учитываются неконсервативные члены,
                         // иначе - расщепление конвективной части уравнений для фаз
    double eps_thin_layer; // точность решения системы нелинейных алгебраических уравнений "тонкого слоя" для проверки
    double eps_contact; // малое значение для задании скорости в случае стационарного контактного разрыва
    double eps_disp_abs; // малое значение объемной доли дисперсной фазы, начиная с которой считается, что она отсутствует
    double eps_cut_out; // малое значение объемной доли дисперсной фазы, начиная с которой при выводе результатов параметры
                        // дисперсной фазы кладуся равными фоновым параметрам
    double eps_relax_compaction; // малое число для обнуления переменных или выражений в релаксации давлений с компактированием,
                                 // равных машинному нулю

    // константы обезразмеривания
    bool use_dimensions; // задача будет определяться размерными параметрами?
    
    double mass_scale; // характерный масштаб массы
    double time_scale; // характерный масштаб времени
    double length_scale; // характерный масштаб длины
    double temperature_scale; // характерный масштаб температуры

    // фоновые параметры для ячеек с малым содержанием дисперсной фазы
    double background_density; // фоновое значение плотности дисперсной фазы
    double background_velocity; // фоновое значение скорости дисперсной фазы
    double background_pressure; // фоновое значение давления дисперсной фазы

    bool nice_output; // записывать в ячейки с малым содержанием дисперсной фазы фоновые параметры при выводе результатов?

    // "физика" задачи
    bool is_physics; // будут ли учитываться члены в правых частях, описывающие "физику" задачи?

    int substeps_num; // количество подшагов, на которое разбивается газодинамический шаг, для интегрирования СОДУ,
                      // описывающей межфазное взаимодействие

    double particle_diameter; // диаметр частиц дисперсной фазы

    int numerical_method_physics; // численный метод решения физики

    // сила трения
    bool is_friction_force; // будет ли учитываться сила межфазного трения?
    FrictionForce fric_force_formula_num; // формула для расчета силы трения
    double interface_drag_coef; // коэффициент сопротивления
    bool is_Cd_const; 
    double Cd;
    double Cd_formula;

    // компактирование
    bool is_compaction; // будет ли учитываться компактирование?
    double compaction_viscosity; // вязкость компактирования

    // горение
    bool is_burning; // будет ли учитываться горение?
	int ignition_condition; // условие воспламенения
    BurningLaw burning_formula_num; // формула для расчета горения
    // свойства пороха - формула SEREBRYAKOV:
    double e_coef; // половина толщины свода горения порохового элемента
    double U_coef; // множитель, зависящий в общем случае от начальной температуры пороха, параметров обдувающего потока и состава пороха
    double nu_coef; // постоянный показатель степени в законе горения
    double X_coef1; // безразмерный коэффициент формы порохового элемента на первой стадии горения
    double lambda_coef1; // безразмерный коэффициент формы порохового элемента на первой стадии горения
    double mu_coef1; // безразмерный коэффициент формы порохового элемента на первой стадии горения
    double X_coef2; // безразмерный коэффициент формы порохового элемента на второй стадии горения
    double lambda_coef2; // безразмерный коэффициент формы порохового элемента на второй стадии горения
    double mu_coef2; // безразмерный коэффициент формы порохового элемента на второй стадии горения
    double zk; // критерий окончания горения - максимальная относительная толщина горящего свода
    // Формула SCHWENDEMAN_BURNING
    double reaction_rate_prefactor; // предэкспоненциальный множитель в уравнении химической реакции
    double p_ignition; // давление "возгорания". при давлениие больше p_ignition начинается горение
    double heat_release; // тепловыделение q

    // теплоперенос
    bool is_heat_transfer; // будет ли учитываться передача тепла между фазами
    double heat_transfer_coefficient; // коэффициент межфазного теплопереноса
    double T_ignition; //Температура воспламенения   
    // поля, которые заполняются в программе
    char input_file_directory[MAX_STRING_SIZE]; // директория, в которой находится конфигурационный файл задачи
    char output_file_directory[MAX_STRING_SIZE]; // директория, куда будут записаны результаты решения задачи

};

// структура с параметрами вычислительного эксперимента в одномерной постановке
struct Parameters1d {
    
    // поля, которые считываются из конфигурационного файла
       
    // информация о лагранжевых скалярах
    int number_of_scalars; // число скаляров 
    double scalar_values[MAX_IC_BLOCKS_NUMBER][N]; // массив характеристик скаляров

    // расчетная сетка (равномерная)
    double left_boundary_x; // координата левой границы расчетной области
    double right_boundary_x; // координата правой границы расчетной области
    int cells_number; // число ячеек в расчетной области

    int ic_blocks_number; // количество блоков, на которые разбита область для задания начальных условий
    int cell_begin[MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек начала блока (нумерация ячеек начинается с 0)
    int cell_end[MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек конца блока (нумерация ячеек начинается с 0)
    double block_values[MAX_IC_BLOCKS_NUMBER][M1D]; //  векторы примитивных переменных в блоках
    int initial_burning[MAX_IC_BLOCKS_NUMBER]; // наличие в блоке начального горения
    double ignition_time; // время поджига
    double additional_heat; // энергия от поджига
    double mass_input; // внесение массы при поджиге
    double S[MAX_IC_BLOCKS_NUMBER];// Площадь канала в блоке
    // граничные условия
    Boundary_condition left_bc; // граничное условие на левой границе
    Boundary_condition right_bc; // граничное условие на правой границе

    double inflow_values[M1D+N]; // параметры втекающего потока
    double change_inflow_time; // момент времени переключения на другие параметры втекания
    double changed_inflow_values[M1D+N]; // другие параметры втекающего потока после заданного момента времени

    // точное решение, только для БН
    bool build_exact_solution; // требуется ли строить точное решение
    int cells_number_for_exact_solution; // число ячеек для построения точного решения,
                                         // только для начальных условий из двух блоков
    
    double sensors_location[MAX_SENSORS_NUM]; // координаты датчиков

};

// структура с параметрами вычислительного эксперимента для 2d задачи
struct Parameters2d {    
    
    // расчетная сетка (равномерная)
    double left_boundary_x; // абсцисса левой границы расчетной области
    double right_boundary_x; // абсцисса правой границы расчетной области
    double down_boundary_y; // ордината нижней границы расчетной области
    double up_boundary_y; // ордината верхней границы расчетной области
    int cells_number_x; // число ячеек по x
    int cells_number_y; // число ячеек по y
       
    int x_blocks_number; // количество блоков, на которые разбита область для задания начальных условий по x
    int y_blocks_number; // количество блоков, на которые разбита область для задания начальных условий по y

    int x_begin[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек начала блока по x (нумерация ячеек начинается с 0) 
    int x_end[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек конца блока по x (нумерация ячеек начинается с 0)

    int y_begin[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек начала блока по y (нумерация ячеек начинается с 0)
    int y_end[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив номеров ячеек конца блока по y (нумерация ячеек начинается с 0)

    int block_type[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER]; // массив типов блоков

    double block_values[MAX_IC_BLOCKS_NUMBER][MAX_IC_BLOCKS_NUMBER][M2D]; // векторы примитивных переменных в блоках
    
    // граничные условия
    Boundary_condition left_bc;  // граничное условие на левой границе
    Boundary_condition right_bc; // граничное условие на правой границе
    Boundary_condition down_bc;  // граничное условие на нижней границе
    Boundary_condition up_bc;    // граничное условие на верхней границе

    double inflow_values_left_bc[M2D];  // параметры втекающего потока через левую ганицу
    double inflow_values_right_bc[M2D]; // параметры втекающего потока через правую границу
    double inflow_values_down_bc[M2D];  // параметры втекающего потока через нижнюю границу
    double inflow_values_up_bc[M2D];    // параметры втекающего потока через верхнюю границу
        
    // координаты датчиков
    double sensors_location_x[MAX_SENSORS_NUM];
    double sensors_location_y[MAX_SENSORS_NUM];
    
};

// для расщепления по пространственным направлениям
enum Direction2d { X_DIRECTION, Y_DIRECTION };

// структура, определяющий текущий момент расчета - либо момент времени, либо количество шагов
struct TimeMoment {
    int steps_num; // текущее количество шагов
    double curr_t; // текущее время
};

// структура с данными для отладки кода и исследования используемых численных методов
struct DebugInfo {

    int current_cell; // номер текущей ячейки
    double current_cell_x; // координата центра текущей ячейки
    double current_cell_vncons[M]; // вектор неконсервативных переменных в текущей ячейке
    
    // для отладки функций расчета потока
    int neighbour_cell; // номер соседней ячейки
    double neighbour_cell_vncons[M]; // вектор неконсервативных переменных в соседней ячейке

    // для отладки процедуры релаксации давлений
    FILE *relaxation_out; // дескриптор файла для записи информации о работе алгоритма релаксации давлений
    int relaxation_cases[RELAX_CASES_NUM]; // счетчики попаданий в различные ветки алгоритма при расчете релаксации давлений

    // для статистики частоты работы различных частей алгоритма метода Годунова
    FILE *godunov_out; // дескриптор файла для записи информации о работе алгоритма метода Годунова
    int godunov_cases[GODUNOV_CASES_NUM];
};

#endif // __STRUCT_H_