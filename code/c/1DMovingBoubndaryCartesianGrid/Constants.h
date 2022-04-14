#ifndef __CONSTANTS_H_
#define __CONSTANTS_H_

#define MAX_BLOCS_NUM       100

#define MAX_STRING_SIZE     250 // ������������ ���������� �������� � ��������� ����������

#define D                   3 // ����������� ��������� �������

#define _PI_                3.1415926536

enum Type_of_border { WALL = 1,
                      PARAMETERS_BEHIND_SHOCK_WAVE = 2,
                      ZERO_ORDER = 3 };

enum Time_cycle_stop_choice { ITERATIONS = 1,       // ���� ������������ �������� ����� ����������
                              FINAL_TIME = 2 };     // ���� ������������ �� ���������� �������� ���������� ��������

enum Time_step_choice { DYNAMIC_TIME_STEP = 1,      // ������������ ����� ���� �� �������
                        CONSTANT_TIME_STEP = 2 };   // ����� ����������� ����

// ������� �����
enum Cell_status { INNER = 1,      // ���������� 
                   BOUNDARY = 2,   // ���������
                   GHOST = 3,      // �����������
                   OUTER = 4 };    // �������

// ����������� ������� ������
enum Direction { FROM_LEFT_TO_RIGHT = 1,    // ����� �������
				 FROM_RIGHT_TO_LEFT = 2 };  // ������ ������

enum Law_of_motion { CONSTANT_VELOCITY = 1,     // ���������� ��������
                     NEWTON = 2 };              // ����� �������

#endif /* __CONSTANTS_H_ */