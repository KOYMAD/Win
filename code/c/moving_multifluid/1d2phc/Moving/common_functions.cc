#include "common_functions.h"

//������� ��� ���������� ������� �������������� ���������� �� ������� ���������������� ����������
//U[D] - ������� �������������� ���������� (out)
//primitive[D] - ������ ���������������� ���������� (in)

void calc_conservative_variables ( struct Parameters *params, struct Primitive_vector primitive, struct Conservative_vector *conservative )
{
    if ( primitive.density <= 0.0 )
    {
        printf( "\ncalc_conservative_variables: negative density\n" );
        printf( "\t primitive parameters \n%e\n%e\n%e\n)", primitive.density, primitive.velosity, primitive.pressure );
        system ( "Pause" );
    }
    if ( primitive.pressure <= 0.0 )
    {
        printf( "\ncalc_conservative_variables: negative pressure\n" );
        printf( "\t primitive parameters \n%e\n%e\n%e\n)", primitive.density, primitive.velosity, primitive.pressure );
        system ( "Pause" );
    }
    conservative->specific_mass = primitive.density;
    conservative->specific_momentum = primitive.density * primitive.velosity;
    conservative->specific_energy = 0.5 * primitive.density * pow( primitive.velosity, 2 ) + ( primitive.pressure / ( params->g - 1 ) );
}

//������� ����������� ������ �������������� ���������� � ������ ���������������� ����������
//urp[D] - (out)
//U[D] - (in)

void calc_primitive_variables ( struct Parameters *params, struct Conservative_vector conservative, struct Primitive_vector *primitive )
{
    primitive->density = conservative.specific_mass;
    if ( conservative.specific_mass > 0.0 )
    {
        primitive->velosity = conservative.specific_momentum / conservative.specific_mass;
    }
    else
    {
        printf( "\ncalc_primitive_variables: negative specific mass\n" );
        printf( "\t conservative parameters \n%e\n%e\n%e\n)", conservative.specific_mass, conservative.specific_momentum, conservative.specific_energy );
        system ( "Pause" );
    }
    primitive->pressure = ( params->g - 1 ) * ( conservative.specific_energy - 0.5 * conservative.specific_mass * pow(primitive->velosity,2) );
    if ( primitive->pressure <= 0.0 )
    {
        printf( "\ncalc_primitive_variables: negative specific energy\n" );
        printf( "\t conservative parameters \n%e\n%e\n%e\n)", conservative.specific_mass, conservative.specific_momentum, conservative.specific_energy );
        system ( "Pause" );
    }
}

// ������ �������� �����
// params - ��������� � ����������� ��������������� ������������
// v_ncons - ������ ����������� ����������
// c - �������� �����
void calc_sound_velocity ( const struct Parameters* params, struct Primitive_vector primitive, double* c )
{
    *c = sqrt( params->g * primitive.pressure / primitive.density );
    if ( primitive.pressure < 0.0 || primitive.density < 0.0 )
    {
        printf( "\ncalc_sound_velocity -> pressure or density are smaller then zero\n" );
        system ( "Pause" );
    }
}

/* ������� max, ������������ ���������� �� ���� �����
   a - ������ ����� (in)
   b - ������ ����� (in) */
double max( double a, double b ) {

    return ( a < b ) ? b : a;

}

/* ������� min, ������������ ���������� �� ���� �����
   a - ������ ����� (in)
   b - ������ ����� (in) */
double min( double a, double b ) {

    return ( a < b ) ? a : b;

}