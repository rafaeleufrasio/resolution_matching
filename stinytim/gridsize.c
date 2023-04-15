#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stinytim.h"

#define NUM_GRID_SIZES  3
static int  GRID_SIZES[] = { 512, 1024, 1280 };
 
/*-------------------------------------------------------------------------------------------*/
float Max_psf_diameter( paramstruct *paraminfo )
{
	float	nyquist_arcsec;

	/* Compute the maximum PSF diameter in arcsec that Tiny Tim can compute *
	 * for a given wavelength and instrument configuration.			*/
	/* Assume shortest wavelength is in paraminfo->wavelength[0]; *
	 * nyquist_arcsec = Nyquist sampling spacing in arcseconds    */

	nyquist_arcsec = 3600.0 * 360.0 * paraminfo->wavelength[0] / 
				(4 * M_PI * (paraminfo->diam_telescope_mm * 1.0e3));

	/* Nyquist-sampled PSF can get grotty beyond central 1/2,  * 
	 * but let it go a bit past that.			   */

	return( GRID_SIZES[NUM_GRID_SIZES-1] * nyquist_arcsec * 0.75 );

} /* Max_psf_diameter */ 

/*-------------------------------------------------------------------------------------------*/
void Set_grid_size( float psf_diam_arcsec, paramstruct *paraminfo )
{
	float	nyquist_arcsec, nyquist_diam_arcsec;
	int	i;

	/* assume shortest wavelength is in paraminfo->wavelength[0] */

	nyquist_arcsec = 3600.0 * 360.0 * paraminfo->wavelength[0] / 
			   (4 * M_PI * (paraminfo->diam_telescope_mm * 1.0e3));

	/* look for smallest grid size that provides a "psf_diam_arcsec" psf */

	for ( i = 0; i < NUM_GRID_SIZES; ++i )
	{
		nyquist_diam_arcsec = GRID_SIZES[i] * nyquist_arcsec * 0.75;
		if ( psf_diam_arcsec < nyquist_diam_arcsec )
		{
			paraminfo->nyquist_grid_size = GRID_SIZES[i];
			break;
		}
	}

	paraminfo->psf_grid_size = (int)(psf_diam_arcsec / paraminfo->det_pixel_size_arcsec_x);

} /* Get_grid_size */

