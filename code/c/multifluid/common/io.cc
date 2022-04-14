// io.cc
// Общие функции ввода/вывода, не зависящие от размерности задачи
// (c) Уткин Павел, 2018
// Создан: 26 августа 2018 г.

#include "io.h"

// Заполнение структуры с основными параметрами вычислительного эксперимента
// params - дескриптор конфигурационного файла задачи
// paramsc - структура с основными параметрами вычислительного эксперимента
void fill_parameters_common( FILE *params, struct ParametersCommon* paramsc ) {

    int dtmp; // для считывания целочисленных переменных
    char string[MAX_STRING_SIZE]; // для считывания строковой информации из файла

    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->program_name) );
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок, обозначающий, что считывается общая часть, не зависящая от размерности задачи

    // модель среды
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->media_model) ); // модель двухфазной среды
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // дисперсная фаза - неподвижный пористый каркас или нет?
    switch ( dtmp ) {
        case 0:
            paramsc->porous_body = false;
            break;
        case 1:
            paramsc->porous_body = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> porous_body should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }

    // шаг по времени
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // шаг - постоянный или нет
    switch ( dtmp ) {
        case 0:
            paramsc->constant_time_step = false;
            break;
        case 1:
            paramsc->constant_time_step = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> constant_time_step should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->cfl) ); // число cfl
    if ( !paramsc->constant_time_step ) {
        if ( paramsc->cfl <= 0.0 ) {
            printf( "\nfill_parameters_struct -> cfl number should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->dt) ); // постоянный шаг по времени
    if ( paramsc->constant_time_step ) {
        if ( paramsc->dt <= 0.0 ) {
            printf( "\nfill_parameters_struct -> dt should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }

    // численный метод
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->numerical_method) ); // номер разностной схемы
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->max_iter_num) ); // максимальное количество итераций для метода Голунова
    if ( paramsc->numerical_method == GODUNOV ) {
        if ( paramsc->max_iter_num <= 0 ) {
            printf( "\nfill_parameters_struct -> max_iter_num should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p_max_ratio) );  // максимальный перепад по давлению слева и справа от разрыва, при котором еще
                                                                                     // используется начальное приближение из решения линеаризованной задачи
                                                                                     // в методе Годунова
    if ( paramsc->numerical_method == GODUNOV ) {
        if ( paramsc->p_max_ratio <= 0.0 ) {
            printf( "\nfill_parameters_struct -> p_max_ratio should be positive.\n\n" );
            exit( EXIT_FAILURE );
        }
    }

    // релаксация давлений
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // подзаголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // релаксацию давлений фаз делать или нет? Только для SA и BN с компактированием
    switch ( dtmp ) {
        case 0:
            paramsc->pressure_relaxation = false;
            break;
        case 1:
            paramsc->pressure_relaxation = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> pressure_relaxation should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // учитывать ли компактирование в операторе релаксации давлений?
    switch ( dtmp ) {
        case 0:
            paramsc->pressure_relaxation_compaction = false;
            break;
        case 1:
            paramsc->pressure_relaxation_compaction = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> pressure_relaxation_compaction should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // номер формулы для расчета конфигурационного давления в релаксации давлений с компактированием
    if ( paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction ) {
        switch ( dtmp ) {
            case 1:
                paramsc->pressure_relaxation_compaction_formula = SCHWENDEMAN_COMPACTION;
                break;
            case 2:
                paramsc->pressure_relaxation_compaction_formula = SAUREL_COMPACTION;
                break;
            default:
                printf( "\nfill_parameters_struct -> pressure_relaxation_compaction_formula should be 1 or 2.\n\n" );
                exit( EXIT_FAILURE );
        }
    } 
    if (paramsc->program_name == TWOD2PHC && paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction && paramsc->pressure_relaxation_compaction_formula == SCHWENDEMAN_COMPACTION){
        printf( "\nfill_parameters_struct -> pressure_relaxation_compaction_formula should be 2. 1 has not been implemented for 2d2phc yet.\n\n" );
        exit( EXIT_FAILURE );
    }

    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->tau_parameter) ); 
    printf( "\nfill_parameters_struct -> tau_parameter = %lf.\n\n", paramsc->tau_parameter );
    if (paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction && paramsc->pressure_relaxation_compaction_formula == SAUREL_COMPACTION){
        if (paramsc->tau_parameter <= 0.0){
            printf( "\nfill_parameters_struct -> tau_parameter should be positive\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->n_parameter) ); 
    if (paramsc->pressure_relaxation && paramsc->pressure_relaxation_compaction && paramsc->pressure_relaxation_compaction_formula == SAUREL_COMPACTION){
        if (paramsc->n_parameter <= 1.0){
            printf( "\nfill_parameters_struct -> n_parameter should be greater than 1.0\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->volume_fraction_compaction) ); 
    if (paramsc->pressure_relaxation_compaction){
        if (paramsc->volume_fraction_compaction < 0.0 || paramsc->volume_fraction_compaction > 1.0){
            printf( "\nfill_parameters_struct -> volume_fraction_compaction should be equal to [0;1]\n\n" );
            exit( EXIT_FAILURE );
        }
    }

    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // релаксацию скоростей фаз делать или нет?
    switch ( dtmp ) {
        case 0:
            paramsc->velocity_relaxation = false;
            break;
        case 1:
            paramsc->velocity_relaxation = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> velocity_relaxation should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->approximation_order) ); // порядок аппроксимации метода
    if ( paramsc->approximation_order != 1 && paramsc->approximation_order != 2 ) {
        printf( "\nfill_parameters_struct -> approximation_order should be 1 or 2\n\n" );
        exit( EXIT_FAILURE );
    }
            
    // уравнения состояния фаз
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // применять реальное УРС для вычисления параметров двучленного УРС или нет?
    switch ( dtmp ) {
        case 0:
            paramsc->use_real_eos = false;
            break;
        case 1:
            paramsc->use_real_eos = true;
            break;
        default:
            printf( "\nfill_parameters_struct -> use_real_eos should be 0 or 1\n\n" );
            exit( EXIT_FAILURE );
    }

    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок подраздела
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->g1) ); // показатель адиабаты в дисперсной фазе
    if ( paramsc->g1 <= 1.0 ) {
        printf( "\nfill_parameters_struct -> g1 should be grater than 1.0.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p01) ); // параметр p01 в УРС для дисперсной фазы
    if ( paramsc->p01 < 0.0 ) {
        printf( "\nfill_parameters_struct -> p01 should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->material_number1) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->specific_heat_cv1) );
	
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок подраздела
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->g2) ); // показатель адиабаты в газовой фазе
    if ( paramsc->g2 <= 1.0 ) {
        printf( "\nfill_parameters_struct -> g2 should be grater than 1.0.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p02) ); // параметр p02 в УРС для газовой фазы
    if ( paramsc->p02 < 0.0 ) {
        printf( "\nfill_parameters_struct -> p02 should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->b_virial) ); // параметр b в вириальном УРС для газовой фазы
    if ( paramsc->b_virial < 0.0 ) {
        printf( "\nfill_parameters_struct -> b should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mu2) ); // коэффициент динамической вязкости в газе
    if ( paramsc->mu2 < 0.0 ) {
        printf( "\nfill_parameters_struct -> mu_g should be non-negative.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->material_number2) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->specific_heat_cv2) );

    // режим расчета
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp );
    switch ( dtmp ) {
        case 1:
            paramsc->is_debug = true; // расчет с проверками матричных операций и корректности решения систем уравнений
            break;
        case 0:
            paramsc->is_debug = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_debug should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }

    // вывод
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp );
    switch ( dtmp ) {
        case 1:
            paramsc->is_output_on_time = true; // критерий вывода результатов - достижение заданного момента времени
            break;
        case 0:
            paramsc->is_output_on_time = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_output_on_time should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->stop_time) ); // момент времени окончания расчета
    if ( paramsc->is_output_on_time ) {
        if ( paramsc->stop_time <= 0.0 ) {
            printf( "\nfill_parameters_struct -> stop_time should be a positive value.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->stop_steps_number) ); // требуемое количество шагов по времени
    if ( !paramsc->is_output_on_time ) {
        if ( paramsc->stop_steps_number <= 0.0 ) {
            printf( "\nfill_parameters_struct -> stop_steps_number should be a positive value.\n\n" );
            exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->output_number) ); // количество файлов с промежуточными результатами
    if ( paramsc->output_number <= 0 ) {
        printf( "\nfill_parameters_struct -> output_number should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // добавить ли заголовок для визуализации в Tecplot?
    switch ( dtmp ) {
        case 1:
            paramsc->is_tecplot_header = true;
            break;
        case 0:
            paramsc->is_tecplot_header = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_tecplot_header should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    // датчики
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок подраздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // требуется ли использовать датчики в заданных точках?
    switch ( dtmp ) {
        case 1:
            paramsc->are_sensors = true;
            break;
        case 0:
            paramsc->are_sensors = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> are_sensors should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    if ( paramsc->are_sensors ) { 
        fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->sensors_num) ); // фактическое количество датчиков
    }
    if ( paramsc->are_sensors && paramsc->sensors_num <= 0 ) {
        printf( "\nfill_parameters_struct -> sensors_num should be positive.\n\n" );
        exit( EXIT_FAILURE );
    }
   
    // малые константы
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_general) ); // регулярное малое число для сравнения вещественных чисел,
                                                                                    // использования в математических функциях, контроля
                                                                                    // сходимости итерационных процессов (если иное не оговорено особо)
    if ( paramsc->eps_general <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_general should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_ludcmp) ); // малое число в функции ludcmp (math_utilc.cc),
                                                                                   // оставлено как было в оригинальной работе
    if ( paramsc->eps_ludcmp <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_ludcmp should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_decouple) ); // максимально допустимый перепад по объемной доле слева и справа,
                                                                                     // при котором учитываются неконсервативные члены, иначе - расщепление
                                                                                     // конвективной части уравнений для фаз
    if ( paramsc->eps_decouple <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_decouple should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_thin_layer) ); // точность решения системы нелинейных алгебраических уравнений
                                                                                       // "тонкого слоя" для проверки
    if ( paramsc->eps_thin_layer <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_thin_layer should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_contact) ); // малое значение для задании скорости в случае стационарного
                                                                                    // контактного разрыва
    if ( paramsc->eps_contact <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_contact should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_disp_abs) ); // малое значение объемной доли дисперсной
                                                                                     // фазы, начиная с которой считается, что она
                                                                                     // отсутствует
    if ( paramsc->eps_disp_abs <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_disp_abs should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_cut_out) ); // малое значение объемной доли дисперсной
                                                                                    // фазы, начиная с которой при выводе
                                                                                    // результатов параметры дисперсной фазы
                                                                                    // кладутся равными фоновым параметрам
    if ( paramsc->eps_cut_out <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_cut_out should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }

    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->eps_relax_compaction) ); // малое число для обнуления переменных 
                                                                                             // или выражений в релаксации давлений 
                                                                                             // с компактированием, равных машинному нулю
    if ( paramsc->eps_relax_compaction <= 0.0 ) {
        printf( "\nfill_parameters_struct -> eps_relax_compaction should be a positive value.\n\n" ); 
        exit( EXIT_FAILURE );
    }

    // константы обезразмеривания
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // задача будет определяться размерными параметрами или нет?
    switch ( dtmp ) {
        case 1:
            paramsc->use_dimensions = true;
            break;
        case 0:
            paramsc->use_dimensions = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> use_dimensions should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mass_scale) ); // характерный масштаб массы
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->time_scale) ); // характерный масштаб времени
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->length_scale) ); // характерный масштаб длины
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->temperature_scale) ); // характерный масштаб температуры

    // фоновые значения параметров дисперсной фазы для ячеек с малым значением объемной доли дисперсной фазы
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->background_density) ); // фоновое значение плотности дисперсной фазы
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->background_velocity) ); // фоновое значение скорости дисперсной фазы
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->background_pressure) ); // фоновое значение давления дисперсной фазы
   
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // записывать в ячейки с малым содержанием дисперсной фазы
                                                                 // фоновые параметры при выводе результатов?
    switch ( dtmp ) {
        case 1:
            paramsc->nice_output = true;
            break;
        case 0:
            paramsc->nice_output = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> nice_output should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }

    // "физика"
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // будут ли учитываться члены, описывающие физику?
    switch ( dtmp ) {
        case 1:
            paramsc->is_physics = true;
            break;
        case 0:
            paramsc->is_physics = false;
            break;
        default:
            printf( "\nfill_parameters_struct -> is_physics should be 0 or 1.\n\n" );
            exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->substeps_num) ); // количество подшагов, на которое
                                                                                    // разбивается газодинамический шаг, для
                                                                                    // интегрирования СОДУ, описывающей
                                                                                    // межфазное взаимодействие
    if ( paramsc->is_physics && paramsc->substeps_num <= 0.0 ) {
        printf( "\nfill_parameters_struct -> substeps_num should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->particle_diameter) ); // диаметр частиц дисперсной фазы
    if ( paramsc->is_physics && paramsc->particle_diameter <= 0.0 ) {
        printf( "\nfill_parameters_struct -> particle_diameter should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
    
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->numerical_method_physics) );
    if ( paramsc->is_physics && ( paramsc->numerical_method_physics != EXPLICIT_EULER && paramsc->numerical_method_physics != NEWTON ) ){
        printf( "\nfill_parameters_struct -> numerical_method_physics should be 1 or 2.\n\n" );
        exit( EXIT_FAILURE );
    }

    // сила межфазного трения
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // будет ли учитываться сила межфазного трения?
    if ( paramsc->is_physics ) {
        switch ( dtmp ) {
            case 1:
                paramsc->is_friction_force = true;
                break;
            case 0:
                paramsc->is_friction_force = false;
                break;
            default:
                printf( "\nfill_parameters_struct -> is_friction_force should be 0 or 1.\n\n" );
                exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // номер формулы для расчеты силы межфазного трения
    if ( paramsc->is_physics && paramsc->is_friction_force ) {
        switch ( dtmp ) {
            case 1:
                paramsc->fric_force_formula_num = ROGUE;
                break;
            case 2:
                paramsc->fric_force_formula_num = HOUIM;
                break;
            case 3:
                paramsc->fric_force_formula_num = TANINO;
                break;
	    case 4:
		paramsc->fric_force_formula_num = SCHWENDEMAN;
                break;
            case 5:
                paramsc->fric_force_formula_num = HOMENKO;
                break;

            default:
                printf( "\nfill_parameters_struct -> fric_force_formula_num should be 1, 2, 3, 4 or 5.\n\n" );
                exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->interface_drag_coef) );
    if ( paramsc->program_name == ONED3PHC || TWOD2PHC){
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // является ли коэффициент Cd постоянным?
	if ( paramsc->is_physics ) {
	    switch ( dtmp ) {
		case 1:
		    paramsc->is_Cd_const = true;
		    break;
		case 0:
		    paramsc->is_Cd_const = false;
		    break;
		default:
		    printf( "\nfill_parameters_struct -> is_Cd_const should be 0 or 1.\n\n" );
		    exit( EXIT_FAILURE );
	    }
	}
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->Cd) );
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->Cd_formula) );
    }

    // компактирование
    fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
    fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // будет ли учитываться компактирование?
    if ( paramsc->is_physics ) {
        switch ( dtmp ) {
            case 1:
		paramsc->is_compaction = true;
                break;
            case 0:
                paramsc->is_compaction = false;
                break;
            default:
                printf( "\nfill_parameters_struct -> is_compaction should be 0 or 1.\n\n" );
                exit( EXIT_FAILURE );
        }
    }
    fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->compaction_viscosity) );

    // горение и теплоперенос в двухфазной программе
    if (paramsc->program_name == ONED2PHC || paramsc->program_name == TWOD2PHC){
	// горение
	fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // будет ли учитываться горение?

	if ( paramsc->is_physics ) {
	    switch ( dtmp ) {
		case 1:
		    paramsc->is_burning = true;
		    break;
		case 0:
		    paramsc->is_burning = false;
		    break;
		default:
		    printf( "\nfill_parameters_struct -> is_burning should be 0 or 1.\n\n" );
		    exit( EXIT_FAILURE );
	    }
	}
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &(paramsc->ignition_condition)); // выбор условия воспламенения

        fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // номер формулы для расчета горения
        if ( paramsc->is_physics && paramsc->is_burning ) {
            switch ( dtmp ) {
                case 1:
                    paramsc->burning_formula_num = SEREBRYAKOV_BURNING;
                    break;
                case 2:
                    paramsc->burning_formula_num = SCHWENDEMAN_BURNING;
                    break;
                default:
                    printf( "\nfill_parameters_struct -> burning_formula_num should be 1 or 2.\n\n" );
                    exit( EXIT_FAILURE );
            }

        }
        // параметры для формулы SEREBRYAKOV_BURNING
        fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->e_coef) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->U_coef) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->nu_coef) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->X_coef1) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->lambda_coef1) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mu_coef1) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->X_coef2) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->lambda_coef2) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->mu_coef2) );
        fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->zk) );
        // параметры для формулы SCHWENDEMAN_BURNING
        fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->reaction_rate_prefactor) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->p_ignition) );
        // продолжение общих параметров горения
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->heat_release) );
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->T_ignition) );
	// теплоперенос
	fscanf_s( params, "%s", string, MAX_STRING_SIZE ); // заголовок раздела
	fscanf_s( params, "%s %d", string, MAX_STRING_SIZE, &dtmp ); // будет ли учитываться теплоперенос?
	if ( paramsc->is_physics ) {
	    switch ( dtmp ) {
		case 1:
		    paramsc->is_heat_transfer = true;
		    break;
		case 0:
		    paramsc->is_heat_transfer = false;
		    break;
		default:
		    printf( "\nfill_parameters_struct -> is_heat_transfer should be 0 or 1.\n\n" );
		    exit( EXIT_FAILURE );
	    }
	}
	fscanf_s( params, "%s %lf", string, MAX_STRING_SIZE, &(paramsc->heat_transfer_coefficient) );
    }
}

// Запись файла с информацией для возможности продолжения расчета
// files_directory - директория, в которой происходит расчет (in)
// time_moment - строка с информацией о моменте времени или числе шагов, для которых
// в последний раз был записан файл с промежуточными результатами (in)
// filename - имя последнего записанного файла (in)
void write_restart_info( char *files_directory, char *time_moment, char *filename ) {

    FILE *restart_info_file; // дескриптор файла с информацией о последнем доступном файле с промежуточными результатами
    char string[MAX_STRING_SIZE]; // для формирования имени файла

    // информация, необходимая для рестарта, находится в файле restart_info.dat
    strcpy_s( string, files_directory );
    strcat_s( string, "\\restart_info.dat" );
    if ( ( fopen_s( &restart_info_file, string, "wt" ) ) != 0 ) {
        printf( "\nwrite_restart_info -> can't open file %s for writing\n\n", string );
        exit( EXIT_FAILURE );
    }

    // запись информации о моменте времени и об имени последнего файла с промежуточными результатами
    fprintf( restart_info_file, "%s %s", time_moment, filename );

    fclose( restart_info_file );

}	



/* Считывание файла с решением для реализации возможности продолжения прерванного расчета 

   params - структура с параметрами вычислительного эксперимента (in)
   restart_info_file - дескриптор файла с информацией о последнем доступном файле с промежуточными результатами (in)

   **initial_solution - массив векторов в примитивных переменных - начальных условий в каждой ячейке (out) */
void read_solution( struct ParametersCommon *paramsc, struct Parameters1d *params1d, FILE *restart_info_file, double **initial_solution ) {

    // char string[MAX_STRING_SIZE];   /* для считывания строковой информации из файла */

}

// Попытка считать файл с информацией о рестарте и, в случае успеха, получение дескриптора
// файла с последними записанными распределениями
// paramsc - структура с основными параметрами вычислительного эксперимента (in) 
// files_directory - директория, в которой происходит расчет (in)
// time_mom - структура, определяющая текущий момент времени (out)
// file_to_read - дескриптор файла с распределениями для рестарта (out)
void try_to_restart( struct ParametersCommon *paramsc, char *files_directory, struct TimeMoment *time_mom, FILE *file_to_read ) {

    FILE *restart_info_file; // дескриптор файла с информацией о последнем доступном файле с промежуточными результатами
    char string[MAX_STRING_SIZE]; // для формирования имени файла и считывания строковой информации
    char mom_string[MAX_STRING_SIZE]; // для считывания информации о моменте времени, которому соответствуют последние записанные распределения

    // информация, необходимая для рестарта, находится в файле restart_info.dat */
    strcpy_s( string, files_directory );
    strcat_s( string, "\\restart_info.dat" );
    if ( ( fopen_s( &restart_info_file, string, "rt" ) ) != 0 ) {
        // нет файла с информацией о рестарте, начинаем расчет с начальных условий
        paramsc->isContinue = false;
        return;
    }
    else {
        // есть файл с информацией о рестарте, продолжаем расчет
        paramsc->isContinue = true;
        printf( "\n> File with restart data restart_info.dat is detected. " );
    }

    // обработка содержимого файла restart_info.dat
    fscanf_s( restart_info_file, "%s %s", mom_string, MAX_STRING_SIZE, string, MAX_STRING_SIZE );

    // запись текущего момента времени
    if ( paramsc->is_output_on_time )
        time_mom->curr_t = atof( mom_string );
    else
        time_mom->steps_num = atoi( mom_string );

    // открытие файла с записанными распределениями на считывание
    if ( ( fopen_s( &file_to_read, string, "rt" ) ) != 0 ) {
        printf( "\ntry_to_restart -> can't open file %s for writing\n\n", string );
        exit( EXIT_FAILURE );    
    }
    printf( "Data from %s is loading.\n", string );

}
