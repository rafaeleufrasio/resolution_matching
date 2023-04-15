#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "stinytim.h"

#define GET_INPUT    1
#define USE_DEFAULT  2


/*---------------------------------------------------------------------------------------------*/
void Process_command_line( int argc, char *argv[], char *filename, paramstruct *paraminfo )
{
	int	i;
	float	value;


	if ( argc < 2 )
	{
		printf( "stiny1 parameter_file [jitter=x] [wmag=x] [scan=x] [pfile=file]\n" );
		exit(0);
	}

	strcpy( filename, argv[1] );

	for ( i = 2; i < argc; ++i )
	{ 
		if ( Find_string(argv[i], "JITTER") != NULL )
		{
			sscanf( &argv[i][7], "%f", &value );
			paraminfo->jitter = value;
			if ( paraminfo->jitter < 0.0 )
			{
				printf( "ERROR : JITTER < 0\n" );
				exit(0);
			} 
		}
		else if ( Find_string(argv[i], "WMAG") != NULL )
		{
			sscanf( &argv[i][5], "%f", &value );
			paraminfo->wavelength_mag = value;
			if ( paraminfo->wavelength_mag <= 0.0 )
			{
				printf( "ERROR : WMAG <= 0\n" );
				exit(0);
			}
		}
		else if ( Find_string(argv[i], "PFILE") != NULL )
		{
			sscanf( &argv[i][6], "%s", paraminfo->pupilfile );
		}
		else if ( Find_string(argv[i], "SCAN") != NULL )
		{
			/* SCAN parameter used only for MIPS 24, 70 WF, 70 superres modes */

			sscanf( &argv[i][5], "%f", &value );
			paraminfo->scan_offset = value;
		}
		else
		{
			printf( "Unrecognized parameter %s ignored\n", argv[i] );
		}
	}

} /* Process_command_line */

/*---------------------------------------------------------------------------------------------*/
void Ask_wavelength( paramstruct *paraminfo )
{
	paraminfo->wavelength[0] = 0.0;
	while ( paraminfo->wavelength[0] < paraminfo->min_wavelength ||
		paraminfo->wavelength[0] > paraminfo->max_wavelength )
	{
		printf( "\nEnter wavelength in microns (%.1f-%.1f) : ",
		   paraminfo->min_wavelength, paraminfo->max_wavelength );
		scanf( "%f", &paraminfo->wavelength[0] );
	}
	paraminfo->num_waves = 1;
	strcpy( paraminfo->spectrum_file, "(none)" );
}

/*---------------------------------------------------------------------------------------------*/
void Ask_bandpass( paramstruct *paraminfo )
{
	int	choice;
	FILE	*file;

	printf( SEPARATOR );

	choice = 0;
	while ( choice < 1 || choice > 3 )
	{
		printf( "\nSelect one of the following :\n");
		printf( "   1) Use default throughput curve\n" );
		printf( "   2) Specify a single wavelength\n" );
		printf( "   3) Use a user-supplied throughput curve\n\n" );
		printf( "Choice (1-3) : " );
		scanf( "%d", &choice );
	}

	/* do nothing if using default throughput curve; *
	 * things are already set up for it              */
				
	if ( choice == 2 )
	{
		Ask_wavelength( paraminfo );
	}
	else if ( choice == 3 )
	{ 
		printf( "\nEnter name of throughput curve file : " );
		scanf( "%s", paraminfo->bandpass_file );

		/* see if file actually exists */

		file = fopen( paraminfo->bandpass_file, "r" );
		if ( file == NULL )
		{
			printf( "ERROR : Could not open bandpass file %s\n",
					paraminfo->bandpass_file );
			exit(0);
		}
		fclose( file );
	}

} /* Ask_bandpass */

/*---------------------------------------------------------------------------------------------*/
void Init_settings( paramstruct *paraminfo )
{
	paraminfo->adjust_field_aberrations = 1;
	paraminfo->bandpass_file[0] = '\0';
	paraminfo->pupilfile[0] = '\0';
	paraminfo->is_subsampled = 0;
	paraminfo->jitter = 0.0;
	paraminfo->num_positions = 1;
	paraminfo->x[0] = 0;
	paraminfo->y[0] = 0;
	paraminfo->num_waves = 0;
	paraminfo->scan_offset = 0.0;
	paraminfo->smart_skip = 1;
	paraminfo->use_kernel = 1;
	paraminfo->use_map = 1;
	paraminfo->use_new_secondary = 1;
	paraminfo->wavelength_mag = 1.0;
	paraminfo->weight_limit = 0.01;
	paraminfo->write_pupil = 0;
	paraminfo->write_wave = 0;
	paraminfo->write_nyquist_psf = 0;
	paraminfo->x_pixels = 0;
	paraminfo->y_pixels = 0;

	paraminfo->focus = 0.0;
	paraminfo->xcoma = 0.0;
	paraminfo->ycoma = 0.0;
	paraminfo->xastig = 0.0;
	paraminfo->yastig = 0.0;
	paraminfo->spherical = 0.0;

	Init_entry_list();

} /* Init_settings */

/*---------------------------------------------------------------------------------------------*/
void Read_position_file( paramstruct *paraminfo, char *filename )
{
	FILE	*file;
	int	i, num_read;

	file = fopen( filename, "r" );
	if ( file == NULL )
	{
		printf( "ERROR : Could not open position list file %s\n", filename );
		exit(0);
	}

	i = 0;
	while ( !feof(file) )
	{
		if ( i == MAX_POSITIONS )
		{
			printf( "ERROR : Only %d positions allowed.\n", MAX_POSITIONS );
			fclose( file );
			exit(0);
		}

				
		num_read = fscanf( file, "%f %f", &paraminfo->x[i], &paraminfo->y[i] );
		if ( num_read == EOF )
			break;
		else if ( num_read != 2 )
		{
			printf( "ERROR : Problems reading position list file\n" );
			fclose( file );
			exit(0);
		}

		++i;
	}

	fclose( file ); 
	paraminfo->num_positions = i;

} /* Read_position_file */

/*---------------------------------------------------------------------------------------------*/
void Ask_field_position( paramstruct *paraminfo )
{
	float	x, y;
	char	line[MAX_STRING];

	x = paraminfo->min_xfield - 1;
	y = paraminfo->min_yfield - 1;

	printf( SEPARATOR );

	while ( x < paraminfo->min_xfield || x > paraminfo->max_xfield || 
		y < paraminfo->min_yfield || y > paraminfo->max_yfield )
	{
		printf( "\nEnter X & Y offset (arcmin) from the center of the instrument field.\n" );
		printf( "X range is %g to %g arcmin; Y range is %g to %g arcmin.\n", 
			paraminfo->min_xfield, paraminfo->max_xfield,
			paraminfo->min_yfield, paraminfo->max_yfield );
		printf( "Or, enter filename, preceded by @, of list containing positions.\n\n" );
		printf( "Field offset : " );
		scanf( "%s", line );

		if ( line[0] != '@' )
		{
			sscanf( line, "%f", &x );
			scanf( "%f", &y );
			paraminfo->num_positions = 1;
			paraminfo->x[0] = x;
			paraminfo->y[0] = y;
		}
		else
		{
			Read_position_file( paraminfo, &line[1] );
			break;
		}
	}

} /* Ask_field_position */

/*---------------------------------------------------------------------------------------------*/
void Ask_pixel_position( paramstruct *paraminfo )
{
	float	x, y;
	int	xc, yc;
	char	line[MAX_STRING];

	x = -1;
	y = -1;

	printf( SEPARATOR );

	xc = paraminfo->x_pixels/2;
	yc = paraminfo->y_pixels/2;

	while ( x < 0 || x >= paraminfo->x_pixels || y < 0 || y >= paraminfo->y_pixels )
	{

		printf( "\nEnter X & Y position on detector in integer pixels.\n" );
		printf( "X range is 0-%d, Y range is 0-%d pixels.\n", 
			paraminfo->x_pixels-1, paraminfo->y_pixels-1 );
		printf( "Center is (X,Y) = (%d,%d)\n", xc, yc );
		printf( "Or, enter filename, preceded by @, of list containing positions.\n\n" );
		printf( "Position : " ); 

		scanf( "%s", line );
		
		if ( line[0] != '@' )
		{
			sscanf( line, "%f", &x );
			scanf( "%f", &y );
			paraminfo->num_positions = 1;
			paraminfo->x[0] = x;
			paraminfo->y[0] = y;
		}
		else
		{
			Read_position_file( paraminfo, &line[1] );
			break;
		}
	}

} /* Ask_pixel_position */

/*------------------------------------------------------------------------*/
void Ask_IRAC( paramstruct *paraminfo )
{
	Ask_pixel_position( paraminfo );
	Ask_bandpass( paraminfo );
}

/*------------------------------------------------------------------------*/
void Ask_IRS( paramstruct *paraminfo )
{
	if ( paraminfo->camera == IRS_PEAKUP_BLUE || paraminfo->camera == IRS_PEAKUP_RED )
	{
		Ask_bandpass( paraminfo );
	}
	else
	{
		printf( SEPARATOR );
		Ask_wavelength( paraminfo );
	}
}

/*------------------------------------------------------------------------*/
void Ask_MIPS( paramstruct *paraminfo )
{
	int 	camera;

	camera = paraminfo->camera;

	if ( camera == MIPS_24 || camera == MIPS_70_WF || camera == MIPS_70_SUPER )
		Ask_pixel_position( paraminfo );
	else if ( camera == MIPS_160 )
		Ask_field_position( paraminfo );

	if ( camera != MIPS_70_SED )
		Ask_bandpass( paraminfo );
	else
	{
		printf( SEPARATOR );
		Ask_wavelength( paraminfo );
	}
}

/*------------------------------------------------------------------------*/
void Ask_psf_size( paramstruct *paraminfo )
{
	int	subfactor, camera;
	float	diam, max_size, min_size, c_wavelength;
	char	choice[9];


	camera = paraminfo->camera;

	/* SIRTF PSFs have about 95% EE at r=45.69 arcsec at 24.0 microns */
	/* if spectrographic mode, use specified wavelength, else use mean wavelength */

	if ( paraminfo->num_waves == 1 )
		c_wavelength = paraminfo->wavelength[0];
	else
		c_wavelength = (paraminfo->wavelength[0] + paraminfo->wavelength[paraminfo->num_waves-1]) / 2;

	printf( SEPARATOR );
	printf( "At wavelength of %.1f microns, 95%% encircled energy \n", c_wavelength );
	printf( "at a radius of about %.3f arcminutes.\n\n", c_wavelength/24.0*45.69/60.0);
 
	max_size = floor(Max_psf_diameter(paraminfo)*10.0) / 10.0;

	diam = 0.0;
	min_size = (0.2 / 3.15) * paraminfo->wavelength[0];  /* arcmin */

	if ( max_size <= 125.0 )
	{
		printf( "Maximum computable PSF diameter = %.1f arcseconds\n", max_size );
		printf( "Minimum recommended PSF diameter = %.2f arcseconds\n", min_size*60.0 ); 
		printf( "Nominal pixel size = %.2f x %.2f arcseconds\n", 
			paraminfo->det_pixel_size_arcsec_x, paraminfo->det_pixel_size_arcsec_y );
		if ( paraminfo->x_pixels != 0 )
			printf( "Detector field size = %.2f by %.2f arcseconds\n\n",
				paraminfo->det_pixel_size_arcsec_x * paraminfo->x_pixels, 
				paraminfo->det_pixel_size_arcsec_y * paraminfo->y_pixels );

		while ( diam < 3*paraminfo->det_pixel_size_arcsec_x || diam > max_size )
		{
			printf( "Enter approximate diameter of PSF model in arcseconds : " );
			scanf( "%f", &diam );
		}
	}
	else
	{
		printf( "Maximum computable PSF diameter = %.3f arcminutes\n", max_size/60.0 );
		printf( "Minimum recommended PSF diameter = %.2f arcminutes\n", min_size ); 
		printf( "Nominal pixel size = %.4f x %.4f arcminutes\n", 
			paraminfo->det_pixel_size_arcsec_x/60.0, paraminfo->det_pixel_size_arcsec_y/60.0 );
		if ( paraminfo->x_pixels != 0 )
			printf( "Detector field size = %.2f by %.2f arcminutes\n\n",
				paraminfo->det_pixel_size_arcsec_x * paraminfo->x_pixels / 60.0, 
				paraminfo->det_pixel_size_arcsec_y * paraminfo->y_pixels / 60.0 );

		while ( diam < 3*paraminfo->det_pixel_size_arcsec_x || diam > max_size )
		{
			printf( "Enter approximate diameter of PSF model in arcminutes : " );
			scanf( "%f", &diam );
			diam = diam * 60.0;
		}
	}

	Set_grid_size( diam, paraminfo );

	printf( "\nDo you want to generate an subsampled PSF? (Y/N) : " );
	scanf( "%s", choice );

	if ( choice[0] == 'Y' || choice[0] == 'y' )
	{
		subfactor = 0;
		while ( subfactor < 2 || subfactor > 10 )
		{
			printf( "\nEnter an INTEGER subsampling factor (2-10) : " );
			scanf( "%d", &subfactor );
  		}

		paraminfo->psf_grid_size *= subfactor;
		paraminfo->psf_pixel_size_arcsec_x = paraminfo->det_pixel_size_arcsec_x / subfactor;
		paraminfo->psf_pixel_size_arcsec_y = paraminfo->det_pixel_size_arcsec_y / subfactor;
		if ( max_size <= 125.0 )
			printf("\nPSF model pixel size = %g arcsec\n", paraminfo->psf_pixel_size_arcsec_x );
		else
			printf("\nPSF model pixel size = %g arcmin\n", paraminfo->psf_pixel_size_arcsec_x/60 );

		paraminfo->is_subsampled = 1;
	}
	else
	{
		paraminfo->psf_pixel_size_arcsec_x = paraminfo->det_pixel_size_arcsec_x;
		paraminfo->psf_pixel_size_arcsec_y = paraminfo->det_pixel_size_arcsec_y;
	}

	printf( "\nPSF diameter is %d pixels (%g arcsec, %g arcmin)\n",
		paraminfo->psf_grid_size, 
		paraminfo->psf_grid_size*paraminfo->psf_pixel_size_arcsec_x,
		paraminfo->psf_grid_size*paraminfo->psf_pixel_size_arcsec_x/60.0 );

} /* Ask_psf_size */
 
/*------------------------------------------------------------------------*/
void Ask_spectrum( paramstruct *paraminfo )
{
	int	choice, i;
	float	tot;

	printf( SEPARATOR );

	choice = 0;
	while ( choice < 1 || choice > 4 )
	{
		printf( "\nSelect object spectrum :\n" );
		printf( "  1) Blackbody \n" );
		printf( "  2) Power law : F(nu) = nu^i \n" );
		printf( "  3) Power law : F(lambda) = lambda^i \n" );
		printf( "  4) User-provided spectrum \n\n" );
		printf( "Choice (1-4) : " );
		scanf( "%d", &choice );
	}

	switch ( choice )
	{
		case 1 : Blackbody( paraminfo );
			 break;
		case 2 : Power_law_nu( paraminfo );
			 break;
		case 3 : Power_law_lambda( paraminfo );
			 break;
		case 4 : printf( "\nNOTE : Spectrum file must contain wavelength (microns)\n" );
			 printf( "       and flux pairs (see manual).\n\n" );
			 printf( "Enter name of spectrum file : " );
			 scanf( "%s", paraminfo->spectrum_file );
			 Spectrum_file( paraminfo );
			 break;
	       default : printf( "ERROR : Invalid selection\n" );
			 exit(0);
	}

	/* normalize weights to a total of 1.0 */

        tot = 0.0;
        for ( i = 0; i < paraminfo->num_waves; ++i )
                tot = tot + paraminfo->weight[i];
        for ( i = 0; i < paraminfo->num_waves; ++i )
                paraminfo->weight[i] /= tot;

} /* Ask_spectrum */

/*------------------------------------------------------------------------*/
void Ask_general( paramstruct *paraminfo )
{
	if ( paraminfo->num_waves > 1 )
		Ask_spectrum( paraminfo );
	Ask_psf_size( paraminfo );

	printf( SEPARATOR );
	printf( "\nEnter rootname of output image file(s) : " );
	scanf( "%s", paraminfo->rootname );
}

/*------------------------------------------------------------------------*/
int main( int argc, char *argv[] )
{
	paramstruct paraminfo;
	char	parameter_file[MAX_STRING];

	if ( argc < 2 )
	{
		printf( "Call is : stiny1 paramfile\n" );
		exit(0);
	}

	Init_settings( &paraminfo );
	Process_command_line( argc, argv, parameter_file, &paraminfo );
 
	printf("\n\n                    Tiny Tim / Spitzer \n");
	printf("                       Version %s \n", VERSION );
	printf("                          June 2006 \n");
	printf("                  Developed by John Krist \n");
	printf("               For the Spitzer Science Center \n\n" );
	paraminfo.camera = 0;

	while ( paraminfo.camera < 1 || paraminfo.camera > IRS_PEAKUP_RED )
	{
		printf( "---------- IRAC ----------     ------------ MIPS ---------- \n" );
		printf( "1) Channel 1 (3.5 microns)     5) 24 microns \n" ); 
		printf( "2) Channel 2 (4.5 microns)     6) 70 microns (Wide Field) \n" ); 
		printf( "3) Channel 3 (6.3 microns)     7) 70 microns (Narrow Field) \n" ); 
		printf( "4) Channel 4 (8.0 microns)     8) 70 microns (S.E.D.) \n" );
		printf( "                               9) 160 microns \n" );
		printf( "--------- IRS ------------ \n" );
		printf( "10) Long-Low \n");
		printf( "11) Long-High \n");
		printf( "12) Short-Low \n");
		printf( "13) Short-High \n");
		printf( "14) Peakup (Blue) \n");
		printf( "15) Peakup (Red) \n\n");
		printf( "Selection : " );
		scanf( "%d", &paraminfo.camera );

		if ( paraminfo.camera < 1 || paraminfo.camera > IRS_PEAKUP_RED )
			printf( "ERROR : Invalid selection\n\n" );
	}

	Read_pupil_file( &paraminfo );
 
	if ( paraminfo.camera >= IRAC_35 && paraminfo.camera <= IRAC_80 )
		Ask_IRAC( &paraminfo );
	else if ( paraminfo.camera >= MIPS_24 && paraminfo.camera <= MIPS_160 )
		Ask_MIPS( &paraminfo );
	else if ( paraminfo.camera >= IRS_LONG_LOW && paraminfo.camera <= IRS_PEAKUP_RED )
		Ask_IRS( &paraminfo );

	Ask_general( &paraminfo );

	Write_parameters( parameter_file, &paraminfo );

	printf( SEPARATOR );
	printf( "\nTo generate PSF, enter the following command :\n" );
	printf( "           stiny2 %s\n\n", argv[1] );

	return( 0 );

} /* main */

