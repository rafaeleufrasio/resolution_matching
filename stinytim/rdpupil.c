#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stinytim.h"

/*---------------------------------------------------------------------------
*  Read_ota :
*     Read parameters from the pupil or parameter file pertaining to the
*     telescope.
*---------------------------------------------------------------------------*/
void Read_ota( FILE *file, paramstruct *paraminfo )
{
        sscanf( Get_entry(file), "%f", &paraminfo->diam_telescope_mm );
        sscanf( Get_entry(file), "%f", &paraminfo->secondary_radius );
	sscanf( Get_entry(file), "%f", &paraminfo->ota_rotation );
}

/*---------------------------------------------------------------------------
*  Read_general :
*	Read parameters (focal ratio, detector size) that are common to all
*	cameras.
*---------------------------------------------------------------------------*/
void Read_general( FILE *file, paramstruct *paraminfo )
{
	sscanf( Get_entry(file), "%f %f", &paraminfo->min_wavelength,
		&paraminfo->max_wavelength );

	sscanf( Get_entry(file), "%f", &paraminfo->det_pixel_size_arcsec_x );
	sscanf( Get_entry(file), "%f", &paraminfo->det_pixel_size_arcsec_y );

	/* x_pixels, y_pixels are the number of detector pixels in x and y */

	if ( (paraminfo->camera != MIPS_70_SED && paraminfo->camera < IRS_LONG_LOW) ||
             (paraminfo->camera == IRS_PEAKUP_BLUE || paraminfo->camera == IRS_PEAKUP_RED) )
		sscanf( Get_entry(file), "%d %d", &paraminfo->x_pixels, &paraminfo->y_pixels );
 
	/* x_offset, y_offset are the field offsets of the detector center *
	 * in arcminutes in the telescope focal plane.			   */

	sscanf( Get_entry(file), "%f %f", &paraminfo->x_offset, &paraminfo->y_offset );
	paraminfo->xfield = paraminfo->x_offset;
	paraminfo->yfield = paraminfo->y_offset;

	sscanf( Get_entry(file), "%f", &paraminfo->detector_rotation );
	sscanf( Get_entry(file), "%d", &paraminfo->detector_x_flip );
	sscanf( Get_entry(file), "%d", &paraminfo->detector_y_flip );

} /* Read_general */

/*-----------------------------------------------------------------------------
*  Read_field_aberrations :
*	Read the coefficient matrices for the field-dependent aberration
*	relations.  MIPS has its own set that are read with a different 
*	routine.
*----------------------------------------------------------------------------*/
void Read_field_aberrations( FILE *file, paramstruct *paraminfo )
{
        float   *c;
        int     j, z;


        for ( z = 0; z <= 5; ++z )
        {
                for ( j = 0; j <= 3; ++j )
                {
                        c = &paraminfo->caber[z][j][0]; 
                        sscanf( Get_entry(file), "%f %f %f %f",
                                &c[0], &c[1], &c[2], &c[3] );
                }
        }

} /* Read_field_aberrations */

/*---------------------------------------------------------------------------
*  Read_IRAC :
*	Read parameters specific to the IRAC cameras
*---------------------------------------------------------------------------*/
void Read_IRAC( FILE *file, paramstruct *paraminfo )
{
	int	j;

	/* read field-dependent aberration equations */

	Read_field_aberrations( file, paraminfo );

	/* read charge diffusion kernel (IRAC 1 & 2 only) */

	if ( paraminfo->camera == IRAC_35 || paraminfo->camera == IRAC_45 )
	{
		for ( j = 0; j < 3; ++j )
	    		sscanf( Get_entry(file), "%f %f %f", &paraminfo->kernel[j][0], 
				&paraminfo->kernel[j][1],&paraminfo->kernel[j][2] );
	}

} /* Read_IRAC */

/*---------------------------------------------------------------------------
*  Read_MIPS :
*	Read parameters pertaining to the MIPS cameras; specifically,
*	the field-dependent aberration matrices.
*---------------------------------------------------------------------------*/
void Read_MIPS( FILE *file, paramstruct *paraminfo )
{
	float 	*c;
	int	camera, iscan, iaber, j;


	camera = paraminfo->camera;

	if ( camera == MIPS_70_SED )
		return;
	else if ( camera == MIPS_160 )
	{
		sscanf( Get_entry(file), "%f %f", 
			&paraminfo->min_xfield, &paraminfo->max_xfield );
		sscanf( Get_entry(file), "%f %f", 
			&paraminfo->min_yfield, &paraminfo->max_yfield );
		Read_field_aberrations( file, paraminfo );
		return;
	} 

	sscanf( Get_entry(file), "%f", &paraminfo->scan_angle );

	/* There is a separate matrix for each field-dependent aberration *
	 * for each reference scan position.				  */
  
	for ( iscan = 0; iscan <= 2; ++iscan )   /* three reference positions */
	{
	   /* mipsscan[] is the reference scan position (field offset   *
	    * in arcsec							*/

	   sscanf( Get_entry(file), "%f", &paraminfo->mipsscan[iscan] );

	   /* read aberrations at the center of field field at each  *
	    * reference scan position.				     */

	   for ( iaber = 0; iaber <= 5; ++iaber )
		sscanf( Get_entry(file), "%f", &paraminfo->mipsaber[iscan][iaber] );

	   /* read field-dependent aberration matrices for each  *
	    * reference scan position.				 */

	   for ( iaber = 0; iaber <= 5; ++iaber )
	   {
		for ( j = 0; j <= 3; ++j )
		{
                        c = &paraminfo->mipscaber[iscan][iaber][j][0]; 
                        sscanf( Get_entry(file), "%f %f %f %f",
                                &c[0], &c[1], &c[2], &c[3] );
		}
	   }
	}

} /* Read_MIPS */

/*-------------------------------------------------------------------------*/
void Read_IRS( FILE *file, paramstruct *paraminfo )
{
	return;
}

/*---------------------------------------------------------------------------
*  Interpolate_bandpass :
*	Interpolate the throughput curve given in paraminfo->wavelength,weight
*	to provide finer or coarser sampling.  The subsampling factor is
*	given in paraminfo->wavelength_mag, which is >1 for subsampling and
*	<1 for undersampling.
*----------------------------------------------------------------------------*/
void Interpolate_bandpass( paramstruct *paraminfo )
{
	float	old_wave[MAX_WAVELENGTHS], old_weight[MAX_WAVELENGTHS];
	float	minw, maxw, w, dw, slope, yint;
	int	i, j1, j2, new_num_waves, old_num_waves;


	old_num_waves = paraminfo->num_waves;

	for ( i = 0; i < old_num_waves; ++i )
	{ 
		old_wave[i] = paraminfo->wavelength[i];
		old_weight[i] = paraminfo->weight[i];
	}

	minw = old_wave[0];
	maxw = old_wave[old_num_waves-1];

	new_num_waves = (int)(paraminfo->wavelength_mag * old_num_waves);
	if ( new_num_waves < 1 )
	{
		new_num_waves = 1;
		paraminfo->wavelength[0] = (minw + maxw) / 2.0;
		paraminfo->weight[0] = 1.0;
		return;
	}
	else if ( new_num_waves > MAX_WAVELENGTHS )
		new_num_waves = MAX_WAVELENGTHS;

	dw = (maxw - minw) / (new_num_waves - 1);

	j1 = 0;
	j2 = 0;

	for ( i = 0; i < new_num_waves; ++i )
	{
		w = i * dw + minw;
		while ( old_wave[j2] < w && j2 < old_num_waves-1 )
			++j2;
		if ( j2 == 0 )
			j2 = 1;
		else
			j1 = j2 - 1;

		slope = (old_weight[j2] - old_weight[j1]) / (old_wave[j2] - old_wave[j1]);
		yint = old_weight[j2] - slope * old_wave[j2];
		paraminfo->wavelength[i] = w;
		paraminfo->weight[i] = slope * w + yint;
	}

	paraminfo->num_waves = new_num_waves;

} /* Interpolate_bandpass */ 
 
/*---------------------------------------------------------------------------
*  Read_throughput :
*	Read the bandpass file.
*---------------------------------------------------------------------------*/
void Read_throughput( paramstruct *paraminfo )
{
	FILE	*file;
	int	i, status;

	file = fopen( paraminfo->bandpass_file, "r" );
	if ( file == NULL )
	{
		printf( "ERROR : Could not open throughput file %s\n", 
			paraminfo->bandpass_file );
		exit(0);
	}

	status = fscanf( file, "%d", &paraminfo->num_waves );
	if ( status != 1 )
	{
		printf( "ERROR : First entry in throughput file is not an integer\n" );
		fclose( file );
		exit(0);
	}

	for ( i = 0; i < paraminfo->num_waves; ++i )
	{
		status = fscanf( file, "%f %f", &paraminfo->wavelength[i], &paraminfo->weight[i] );
		if ( status != 2 )
		{
			printf( "ERROR : Bad wavelength & weight pair in throughput file\n" );
			fclose( file );
			exit(0);
		}
	}

	fclose(file);

	if ( paraminfo->wavelength_mag != 1.0 )
		Interpolate_bandpass( paraminfo );

} /* Read_throughput */

/*---------------------------------------------------------------------------
*  Read_aberrations :
*	Read the Zernike coefficients for the center of the detector
*	field.  Aberrations are in microns RMS for 33% obscured
*	Zernikes.  These values are added to the default field-center
*	aberrations.
*---------------------------------------------------------------------------*/
void Read_aberrations( FILE *file, paramstruct *paraminfo )
{
	int	i;

	for ( i = 1; i <= LAST_ZERNIKE; ++i )
		sscanf( Get_entry(file), "%f", &paraminfo->zval[i] );

} /* Read_aberrations */

/*---------------------------------------------------------------------------
*  Read_instrument :
*	Call the proper routine to read the instrument parameters.
*--------------------------------------------------------------------------*/
void Read_instrument( FILE *file, paramstruct *paraminfo )
{
	switch ( paraminfo->camera )
	{
		case IRAC_35 :
		case IRAC_45 :
		case IRAC_63 :
		case IRAC_80 :	
				Read_IRAC( file, paraminfo );
				break;

		case MIPS_24 :
		case MIPS_70_SUPER :
		case MIPS_70_WF : 
		case MIPS_70_SED :
		case MIPS_160 :
				Read_MIPS( file, paraminfo );
				break;

		case IRS_LONG_LOW :
		case IRS_LONG_HIGH :
		case IRS_SHORT_LOW :
		case IRS_SHORT_HIGH :
		case IRS_PEAKUP_BLUE :
		case IRS_PEAKUP_RED :
				Read_IRS( file, paraminfo );
				break;

	        default :  printf( "ERROR (Read_instrument) : unknown instrument\n" );
			   exit(0);
	}

} /* Read_instrument */

/*-----------------------------------------------------------------------------
*  Read_pupil_file :
* 	Read instrument-specific parameters.
*-----------------------------------------------------------------------------*/
void Read_pupil_file( paramstruct *paraminfo )
{
	FILE	*file;
	char	filename[MAX_STRING];
	char	*rootname[] = { " ", 
			"irac1", "irac2", "irac3", "irac4", 
			"mips24", "mips70wf", "mips70super", "mips70sed", "mips160",
			"irslonglow", "irslonghigh", "irsshortlow", "irsshorthigh",
			"irspeakupblue", "irspeakupred" };
	char	*camera_name[] = { " ",
			"IRAC_3.5_micron", "IRAC_4.5_micron", "IRAC_6.3_micron", "IRAC_8.0_micron",
			"MIPS_24_micron", "MIPS_70_micron_wf", "MIPS_70_micron_super",
			"MIPS_70_micron_sed", "MIPS_160_micron",
			"IRS_long_low", "IRS_long_high", "IRS_short_low", "IRS_short_high",
			"IRS_peakup_blue", "IRS_peakup_red" };


	strcpy( paraminfo->camera_name, camera_name[paraminfo->camera] );

	/* read aberration file (list of Zernike coefficients for the *
	 * center of the instrument's field).			      */

	Default_dir( filename );
	strcat( filename, rootname[paraminfo->camera] );
	strcat( filename, ".zer" );
	file = fopen( filename, "r" );
	if ( file == NULL )
	{
		printf( "ERROR : Could not open Zernike file %s\n", filename );
		exit(0);
	}
	Read_aberrations( file, paraminfo );
	fclose( file );

	/* Read the OTA configuration file (ota.pup) */

        Default_dir( filename );
        strcat( filename, "ota.pup" );
        file = fopen( filename, "r" );
        if ( file == NULL )
        {
                printf( "ERROR : Could not open pupil file %s\n", filename );
                exit(0);
        }
	Read_ota( file, paraminfo );
	fclose( file );

	/* Read the instrument configuration file (ie. .pup file) */

	if ( paraminfo->pupilfile[0] == '\0' )
	{
		Default_dir( paraminfo->pupilfile );
		strcat( paraminfo->pupilfile, rootname[paraminfo->camera] );
		strcat( paraminfo->pupilfile, ".pup" );
	}

	file = fopen( paraminfo->pupilfile, "r" );
	if ( file == NULL )
	{
		printf( "ERROR : Could not open pupil file %s\n", paraminfo->pupilfile );
		exit(0);
	}
	Read_general( file, paraminfo );
	Read_instrument( file, paraminfo );
	fclose( file );

	/* read system throughput curve */

	if ( paraminfo->camera != MIPS_70_SED && 
	     (paraminfo->camera < IRS_LONG_LOW || paraminfo->camera > IRS_SHORT_HIGH) )
	{
		if ( paraminfo->bandpass_file[0] == '\0' )
		{
			Default_dir( paraminfo->bandpass_file );
			strcat( paraminfo->bandpass_file, rootname[paraminfo->camera] );
			strcat( paraminfo->bandpass_file, ".band" );
		}
		Read_throughput( paraminfo );
	}
	else
	{
		strcat( paraminfo->bandpass_file, "(none)" );
		paraminfo->weight[0] = 1.0;
	}

} /* Read_pupil_file */

