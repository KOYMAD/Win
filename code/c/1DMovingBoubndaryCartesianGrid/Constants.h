#ifndef __CONSTANTS_H_
#define __CONSTANTS_H_

#define MAX_BLOCS_NUM       100

#define MAX_STRING_SIZE     250 // максимальное количество символов в строковой переменной

#define D                   3 // размерность основного вектора

#define _PI_                3.1415926536

enum Type_of_border { WALL = 1,
                      PARAMETERS_BEHIND_SHOCK_WAVE = 2,
                      ZERO_ORDER = 3 };

enum Time_cycle_stop_choice { ITERATIONS = 1,       // цикл продолжается заданное число повторений
                              FINAL_TIME = 2 };     // цикл продолжается до достижения временем финального значения

enum Time_step_choice { DYNAMIC_TIME_STEP = 1,      // динамический выбор шага по времени
                        CONSTANT_TIME_STEP = 2 };   // выбор постоянного шага

// статусы ячеек
enum Cell_status { INNER = 1,      // внутренняя 
                   BOUNDARY = 2,   // граничная
                   GHOST = 3,      // виртуальная
                   OUTER = 4 };    // внешняя

// направление расчёта потока
enum Direction { FROM_LEFT_TO_RIGHT = 1,    // слева направо
				 FROM_RIGHT_TO_LEFT = 2 };  // справа налево

enum Law_of_motion { CONSTANT_VELOCITY = 1,     // постоянная скорость
                     NEWTON = 2 };              // закон Ньютона

#endif /* __CONSTANTS_H_ */