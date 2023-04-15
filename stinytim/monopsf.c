/* File      :  monopsf.c
 *
 * Contents  :
 *      Convert_OPD      : Convert OPD from microns to waves of error
 *      Write_pupil_map  : Write out aperture or wavefront error function to a file
 *      Compute_mono_psf : Compute the monochromatic PSF.
 *
 * Author    :  John Krist
 * Date      :  January 2000 (modified for SIRTF)
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "stinytim.h"

/*-------------------------------------------------------------------------
*  Convert_OPD :
*       Convert OPD from microns to waves at "wavelength" microns.
*
*  Inputs :
*           dim : Dimension of "pupil" function array (dim by dim).
*    wavelength : Wavelength in microns of desired OPD.
*         pupil : Complex valued pupil function defined in microns
*
*  Returns :
*       Returns the wavelength converted OPD in the "pupil" array.
*-----------------------------------------------------------------------------*/
void Convert_OPD( int dim, float wavelength, complex **pupil )
{
	int     x, y;

	for ( y = 0; y < dim; ++y )
		for ( x = 0; x < dim; ++x )
			pupil[y][x].i /= wavelength;

} /* Convert_OPD */

/*-----------------------------------------------------------------------------
*  Write_pupil_map :
*       Write the aperture or opd array to a FITS files.
*
*  Inputs :
*	pupil : complex pupil array
*	  dim : dimension of "pupil" (dim by dim)
*   data_type : flag indicating pupil (1) or opd (2) data to be written.
*
*----------------------------------------------------------------------------*/
void Write_pupil_map( complex **pupil, int dim, int data_type )
{
	float	**temp;
	int	i, j;
	char    name[MAX_STRING];


	/* Create a temporary array for the shifted image */

	temp = Alloc_image( dim, dim );

	if ( data_type == PUPIL_FILE )
	{
		for ( j = 0; j < dim; ++j )
			for ( i = 0; i < dim; ++i )
				temp[j][i] = pupil[j][i].r;
		strcpy( name, "tinytim_pupil.fits" );
		printf( " ** Writing aperture pattern to tinytim_pupil.fits **\n" );
	}
	else
	{
		for ( j = 0; j < dim; ++j )
			for ( i = 0; i < dim; ++i )
				temp[j][i] = pupil[j][i].i;
		strcpy( name, "tinytim_opd.fits" );
		printf( " ** Writing OPD to tinytim_opd.fits **\n" );
	}
	fflush( stdout );

	Shift_to_center( temp, dim );

        Write_FITS( name, temp[0], dim, dim, NULL, data_type );

	Free_image( temp );

} /* Write_pupil_map */

/*-----------------------------------------------------------------------------
*  Compute_mono_psf :
*       Compute the monochromatic PSF for a specified wavelength.
*
*  Inputs :
*        paraminfo : parameter information structure
*       wave_index : wavelength index
*          monopsf : 2D array to contain monochromatic PSF; must be allocated
*		 	by calling routine
*            pupil : pupil function array (complex structure);  assumed to be
*                       shifted so that the aperture center is at (0,0)
*
*  Returns :
*       Returns the monochromatic, detector-integrated PSF in "monopsf".
*----------------------------------------------------------------------------*/
void Compute_mono_psf( paramstruct *paraminfo, int wave_index, 
		float **monopsf, complex **pupil )
{
	int	dim, i, j;
	float	**temp, nyquist_scale_arcsec;


	dim = paraminfo->nyquist_grid_size;

	/* Convert OPD to current wavelength */

	Convert_OPD( dim, paraminfo->wavelength[wave_index], pupil );

	/* If the "write_pupil" or "write_wave" flags are set, write the  *
	 * arrays to disk (do this just once).				  */

	if ( paraminfo->write_pupil == 1 )
	{
		Write_pupil_map( pupil, dim, PUPIL_FILE );
		paraminfo->write_pupil = 0;
	}

	if ( paraminfo->write_wave == 1 )
	{
		Write_pupil_map( pupil, dim, OPD_FILE );
		paraminfo->write_wave = 0;
	}

	Compute_nyquist_psf( pupil, dim );   /* "pupil" now contains Nyquist sampled PSF */

	/* If "write_nyquist_psf" is set, write the Nyquist sampled      *
	 * PSF to the file "tinytim_nyquist_psf.fits".  The PSF is first *
	 * shifted to the center.  The program then stops (hoping that   *
	 * the exit() routine will take care of cleaning up the arrays). */

	if ( paraminfo->write_nyquist_psf == 1 )
	{
		temp = Alloc_image( dim, dim );

		for ( j = 0; j < dim; ++j )
			for ( i = 0; i < dim; ++i )
				temp[j][i] = pupil[j][i].r;

		Shift_to_center( temp, dim );

	     	printf("*** Writing Nyquist sampled PSF to tinytim_nyquist_psf.fits ***\n");
	     	Write_FITS( "tinytim_nyquist_psf.fits", temp[0], dim, dim, NULL, NYQUISTPSF_FILE );

		Free_image( temp );
		paraminfo->write_nyquist_psf = 0;
	}

	/* Integrate Nyquist-sampled PSF onto detector pixels */

	nyquist_scale_arcsec = 3600.0 * 360.0 * paraminfo->wavelength[wave_index] / 
				(4 * M_PI * (paraminfo->diam_telescope_mm * 1.0e3));

	Integrate_psf( pupil, nyquist_scale_arcsec, dim, monopsf, 
		       paraminfo->psf_pixel_size_arcsec_x, paraminfo->psf_pixel_size_arcsec_y, 
		       paraminfo->psf_grid_size, paraminfo->jitter );

} /* Compute_mono_psf */

