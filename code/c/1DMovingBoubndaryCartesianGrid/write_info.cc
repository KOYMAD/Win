#include "write_info.h"

void write_results ( struct Parameters *params, int *status, struct Conservative_vector *conservative, struct TimeMoment time_mom )
{
    FILE *out;
    char a[50] = "tecplot_solution_";
    char b[15];
    sprintf(b, "%06d%s", time_mom.steps_num, "_");
    strcat(a,b);
    sprintf(b, "%e", params->time_diml * time_mom.curr_t );
    strcat(a,b);
    char c[5] = ".dat";
    strcat(a,c);
 
    // открытие файла на запись
    out = fopen( a, "wt" );
    if ( NULL == out )
    {
        printf( "\nwrite_results -> can't open file for results binary writing\n" );
        system ( "Pause" );
        exit( EXIT_FAILURE );
    }

    fprintf( out, "TITLE = \"Example: Simple XY Plot\"\nVARIABLES = \"Distance\", \"Density\", \"Velocity\", \"Pressure\"\nZONE T=\"%lf\", I=%i, F=POINT\n", params->time_diml * time_mom.curr_t, params->number_of_cells );

    double length_diml = params->length_diml;
    double time_diml = params->time_diml;
    double dens_diml = params->mass_diml / pow(length_diml, 3);
    double speed_diml = length_diml / time_diml;
    double press_diml = params->mass_diml / ( length_diml * pow(time_diml, 2));
    struct Primitive_vector primitive;
    double grid_step = ( params->coordinate_of_right_boundary - params->coordinate_of_left_boundary )
        / params->number_of_cells;
    for ( int i = 0 ; i < params->number_of_cells ; i++ )
    {
        if ( INNER == status[i] || BOUNDARY == status[i] )
        {
            calc_primitive_variables ( params, conservative[i], &primitive );
            fprintf( out, "%e\t%e\t%e\t%e\n", grid_step * ( i + 0.5 ) * length_diml, primitive.density * dens_diml,
                primitive.velosity * speed_diml, primitive.pressure * press_diml );
        }
        else fprintf( out, "%e\t%e\t%e\t%e\n", grid_step * ( i + 0.5 ) * length_diml, 0.0, 0.0, 0.0 );
    }
    fclose(out);
}

void write_piston_trajectory ( struct Parameters *params, double coordinate_of_left_boundary_of_body,
struct TimeMoment time_mom )
{
    FILE *out;
    //if ( ( out = fopen( "piston_trajectory.dat", "r" ) ) == NULL )
    //    puts ( "piston_trajectory.dat is not exist or be used by another program\n" );
    //else
    //    puts ( "piston_trajectory.dat is open\n" );
    //fclose ( out );
    if ( 0.0 == time_mom.curr_t )
        if ( ( out = fopen ( "piston_trajectory.dat", "wt" ) ) != NULL )
            fprintf ( out, "%e\t%e\n", coordinate_of_left_boundary_of_body, time_mom.curr_t );
        else system ( "Pause" );
    else
        if ( ( out = fopen ( "piston_trajectory.dat", "at" ) ) != NULL )
            fprintf ( out, "%e\t%e\n", coordinate_of_left_boundary_of_body, time_mom.curr_t );
        else system ( "Pause" );

    fclose ( out );
}