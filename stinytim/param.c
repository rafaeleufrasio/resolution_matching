#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include "stinytim.h"

/*-------------------------------------------------------------------------------------*/
void Write_parameters( char *filename, paramstruct *paraminfo )
{
	FILE 	*file;
	int	camera, i;

	camera = paraminfo->camera;

	file = fopen( filename, "w" );
	if ( file == NULL )
	{
		printf( "ERROR : Could not open %s for output.\n", filename );
		exit(0);
	}

	fprintf( file, "%s    = Version of Tiny Tim that generated this file\n", VERSION ); 
	fprintf( file, "%s      = rootname of output image file(s)\n", paraminfo->rootname );
	fprintf( file, "%d         = instrument code\n", paraminfo->camera );
	fprintf( file, "%s        = instrument name\n", paraminfo->camera_name );
	fprintf( file, "%f  = PSF model X-axis pixel size in arcsec (might be subsampled)\n",
		paraminfo->psf_pixel_size_arcsec_x );
	fprintf( file, "%f  = PSF model X-axis pixel size in arcsec (might be subsampled)\n",
		paraminfo->psf_pixel_size_arcsec_y );
	fprintf( file, "%d  = Is PSF subsampled? (0 = no, 1 = yes)\n", paraminfo->is_subsampled );
	fprintf( file, "%d    = model grid size (n by n pixels)\n", 
		paraminfo->psf_grid_size );
	fprintf( file, "%d    = Nyquist PSF grid size (n by n pixels)\n",
		paraminfo->nyquist_grid_size );

	fprintf( file, "%f  = Jitter (arcsec) (ignored if < 0.001)\n", paraminfo->jitter );

	if ( paraminfo->camera == IRAC_35 || paraminfo->camera == IRAC_45 )
	{
		if ( paraminfo->is_subsampled )
			fprintf( file, "0    = Convolve PSF with IRAC charge diffusion kernel? (1=yes)\n" );
		else
			fprintf( file, "0    = Convolve PSF with IRAC charge diffusion kernel? (1=yes)\n" );
	}

	fprintf( file, "%s   = instrument throughput table file used\n",
		paraminfo->bandpass_file );

	fprintf( file, "%s   = spectrum type or spectrum file used\n", 
		paraminfo->spectrum_file );

	fprintf( file, "%d        = Skip wavelengths with low weights? (1=yes)\n",
		paraminfo->smart_skip );
	fprintf( file, "%f  = Weight skip low weight limit/max weight\n",
		paraminfo->weight_limit );

	fprintf( file, "%d    = number of wavelengths in model\n", paraminfo->num_waves );
	for ( i = 0; i < paraminfo->num_waves; ++i )
		fprintf( file, "%f  %g     = Wavelength %d (microns) & weight\n",
			paraminfo->wavelength[i], paraminfo->weight[i], i+1 );

	if ( camera != MIPS_70_SED && camera < IRS_LONG_LOW )
	{ 
		fprintf( file, "%d        = number of detector positions\n", paraminfo->num_positions );
		for ( i = 0; i < paraminfo->num_positions; ++i )
			fprintf( file, "%g %g     = X & Y position %d\n",
				paraminfo->x[i], paraminfo->y[i], i+1 );
		fprintf( file, "%d    = Adjust for field dependent aberrations? (1=yes)\n",
			paraminfo->adjust_field_aberrations );
	}

	fprintf( file, "%d    = Use zonal error maps? (1=yes)\n",
		paraminfo->use_map );
	fprintf( file, "%d    = Use new secondary map? (1=yes)\n",
		paraminfo->use_new_secondary );
	fprintf( file, "%d    = Write pupil map? (1=yes)\n",
		paraminfo->write_pupil );
	fprintf( file, "%d    = Write OPD map? (1=yes)\n",
		paraminfo->write_wave );
	fprintf( file, "%d    = Write Nyquist PSF? (1=yes)\n",
		paraminfo->write_nyquist_psf );

	if ( camera == MIPS_24 || camera == MIPS_70_WF || camera == MIPS_70_SUPER )
		fprintf( file, "%g  = MIPS scan mirror offset (arcsec)\n", paraminfo->scan_offset );
 
	Write_entry_list( file );
	Delete_entry_list();

	fclose( file );

} /* Write_parameters */

/*---------------------------------------------------------------------------------------*/
void Read_parameters( char *filename, paramstruct *paraminfo )
{
	FILE	*file;
	int	camera, i;
	char	file_version[20];


	Init_entry_list();

	file = fopen( filename, "r" );
	if ( file == NULL )
	{
		printf( "ERROR : Could not open %s for input.\n", filename );
		exit(0);
	}

	
	sscanf( Get_entry(file), "%s", file_version );
	if ( strstr(VERSION, file_version) == NULL )
	{
		printf("ERROR : Parameter file not produced by this version of Tiny Tim\n" );
		fclose( file );
		exit(0);
	}
 
	sscanf( Get_entry(file), "%s", paraminfo->rootname );
	sscanf( Get_entry(file), "%d", &paraminfo->camera );
	camera = paraminfo->camera;
	sscanf( Get_entry(file), "%s", paraminfo->camera_name );
	sscanf( Get_entry(file), "%f", &paraminfo->psf_pixel_size_arcsec_x );
	sscanf( Get_entry(file), "%f", &paraminfo->psf_pixel_size_arcsec_y );
	sscanf( Get_entry(file), "%d", &paraminfo->is_subsampled );
	sscanf( Get_entry(file), "%d", &paraminfo->psf_grid_size );
	sscanf( Get_entry(file), "%d", &paraminfo->nyquist_grid_size );
	sscanf( Get_entry(file), "%f", &paraminfo->jitter );

	if ( paraminfo->camera == IRAC_35 || paraminfo->camera == IRAC_45 )
		sscanf( Get_entry(file), "%d", &paraminfo->use_kernel );

	sscanf( Get_entry(file), "%s", paraminfo->bandpass_file );
	sscanf( Get_entry(file), "%s", paraminfo->spectrum_file );
	sscanf( Get_entry(file), "%d", &paraminfo->smart_skip );
	sscanf( Get_entry(file), "%f", &paraminfo->weight_limit );

	sscanf( Get_entry(file), "%d", &paraminfo->num_waves );
	for ( i = 0; i < paraminfo->num_waves; ++i )
		sscanf( Get_entry(file), "%f %f", &paraminfo->wavelength[i], &paraminfo->weight[i] );

	if ( camera != MIPS_70_SED && camera < IRS_LONG_LOW )
	{
		sscanf( Get_entry(file), "%d", &paraminfo->num_positions );
		for ( i = 0; i < paraminfo->num_positions; ++i )
			sscanf( Get_entry(file), "%f %f", &paraminfo->x[i], &paraminfo->y[i] );
		sscanf( Get_entry(file), "%d", &paraminfo->adjust_field_aberrations );
	}
	else
	{
		paraminfo->num_positions = 1;
	}

	sscanf( Get_entry(file), "%d", &paraminfo->use_map );
	sscanf( Get_entry(file), "%d", &paraminfo->use_new_secondary );
	sscanf( Get_entry(file), "%d", &paraminfo->write_pupil );
	sscanf( Get_entry(file), "%d", &paraminfo->write_wave );
	sscanf( Get_entry(file), "%d", &paraminfo->write_nyquist_psf );

	if ( camera == MIPS_24 || camera == MIPS_70_WF || camera == MIPS_70_SUPER )
		sscanf( Get_entry(file), "%f", &paraminfo->scan_offset );
	else
		paraminfo->scan_offset = 0.0;

	Read_aberrations( file, paraminfo );
	Read_ota( file, paraminfo );
	Read_general( file, paraminfo );
	Read_instrument( file, paraminfo );

	fclose( file );

	Delete_entry_list();
}
 
