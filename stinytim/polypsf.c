/*  File     :  polypsf.c
 *
 *  Contents :
 *      Flux_normalize : Normalize PSF to total flux of 1.0
 *      Add_mono_to_poly  : Weighted addition of a monochromatic PSF to the
 *                              current polychromatic PSF array.
 *      Compute_poly_psf : Compute polychromatic PSF for a given position.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/* #include <malloc.h> */
#include "stinytim.h"

/*-------------------------------------------------------------------------*/
void Flux_normalize( float **psf, int n )
{
	int	x, y;
	float 	sum;

	sum = 0.0;

	for ( y = 0; y < n; ++y )
		for ( x = 0; x < n; ++x )
			sum += psf[y][x];

	if ( sum != 0.0 )
	{
		for ( y = 0; y < n; ++y )
			for ( x = 0; x < n; ++x )
				psf[y][x] /= sum;
	}
}

/*-------------------------------------------------------------------------
*  Add_mono_to_poly :
*       Weighted addition of monochromatic PSF to the "polypsf" array.
*
*  Inputs :
*       paraminfo :  parameter information structure
*      wave_index :  current wavelength index
*         monopsf :  monochromatic PSF image
*        polypsf  :  polychromatic PSF image
*
*-------------------------------------------------------------------------*/
void Add_mono_to_poly( paramstruct *paraminfo, int wave_index, float **monopsf, float **polypsf )
{
	int     x, y;
	float   weight;

	/* Flux normalize "monopsf", and then add it to "polypsf"  *
	 * with appropriate weighting.				   */

	Flux_normalize( monopsf, paraminfo->psf_grid_size );

	if ( paraminfo->num_waves < 2 )      /* Only one wavelength */
	{

		for ( y = 0; y < paraminfo->psf_grid_size; ++y )
			for ( x = 0; x < paraminfo->psf_grid_size; ++x )
				polypsf[y][x] = monopsf[y][x];
	}
	else
	{
		weight = paraminfo->weight[wave_index];
		for ( y = 0; y < paraminfo->psf_grid_size; ++y )
			for ( x = 0; x < paraminfo->psf_grid_size; ++x )
				polypsf[y][x] += weight * monopsf[y][x];
	}

} /* Add_mono_to_poly */

/*-------------------------------------------------------------------------
*  Convolve_kernel :
*	Convolve IRAC 1 or 2 PSF with 3 by 3 pixel charge diffusion kernel.
*-------------------------------------------------------------------------*/
void Convolve_kernel( paramstruct *paraminfo, float **psf )
{
	int	first_row, i, j, ix, iy, nx, ny, x, y;
	float	sum, *temp_ptr, *temp[3];


	nx = ny = paraminfo->psf_grid_size;

	for ( j = 0; j < 3; ++j )
	{
		temp[j] = (float *)malloc( nx * sizeof(float) );
		memcpy( &temp[j][0], &psf[j][0], nx*sizeof(float) );
 	}
	first_row = 0;
 
	for ( y = 0; y < ny; ++y )
	{
		if ( y > 1 && y < ny-1 )
		{
			temp_ptr = temp[0];
			temp[0] = temp[1];
			temp[1] = temp[2];
			temp[2] = temp_ptr;
			memcpy( &temp[2][0], &psf[y+1][0], paraminfo->psf_grid_size*sizeof(float) );
			++first_row;
		}

		for ( x = 0; x < nx; ++x )
		{
			sum = 0.0;

			for ( j = 0; j < 3; ++j )
			{
				iy = y - first_row - 1 + j;
				if ( iy < 0 )
					iy = 0;
				else if ( iy > 2 )
					iy = 2;

				for ( i = 0; i < 3; ++i )
				{
					ix = x - 1 + i;
					if ( ix < 0 )
						ix = 0;
					else if ( ix >= nx )
						ix = nx - 1;

					sum += temp[iy][ix] * paraminfo->kernel[j][i];
				}
			}

			psf[y][x] = sum;
		}
	}

	for ( j = 0; j < 3; ++j )
		free( temp[j] );

} /* Convolve_kernel */

/*-------------------------------------------------------------------------
*  Compute_poly_psf :
*       Compute the polychromatic PSF for the current position.
*
*  Inputs :
*	paraminfo : parameter information structure
*	    pupil : complex array of dimension paraminfo->crit_dim^2
*	  polypsf : 2D array to contain final PSF; must be allocated by
*                     calling routine
*         monopsf : 2D array to contain monochromatic PSF; must be
*		      allocated by calling routine
*  Outputs :
*       Returns the polychromatic PSF in "polypsf".
*
*-------------------------------------------------------------------------*/
void Compute_poly_psf( paramstruct *paraminfo, complex **pupil, float **polypsf, float **monopsf )
{
	int     i, iwave, crit_dim, pupil_drawn, x, y;
	int     file_written, skip_psf;
	float   **temp, max_weight, weight_limit;
	char	pupil_name[MAX_STRING];


	/* Define name for the temporary pupil/OPD file */

	strcpy( pupil_name, paraminfo->rootname );
	strcat( pupil_name, ".pupil" );

	for ( y = 0; y < paraminfo->psf_grid_size; ++y )
		for ( x = 0; x < paraminfo->psf_grid_size; ++x )
			polypsf[y][x] = 0.0;

	/* if using "Smart Skip" feature, determine max weight */

	if ( paraminfo->smart_skip != 0 )
	{
		max_weight = paraminfo->weight[0];
		for ( i = 1; i < paraminfo->num_waves; ++i )
		    	if ( paraminfo->weight[i] > max_weight )
				max_weight = paraminfo->weight[i];
		weight_limit = paraminfo->weight_limit * max_weight;
	}

	file_written = 0;
	crit_dim = paraminfo->nyquist_grid_size;
	pupil_drawn = 0;

	/* Compute monochromatic PSF for each wavelength */

	for ( iwave = 0; iwave < paraminfo->num_waves; ++iwave )
	{
                if ( paraminfo->num_waves > 1 )
		{
		   	if ( (paraminfo->smart_skip != 0) && (paraminfo->weight[iwave] < weight_limit) )
		   	{
				printf("   Skipping PSF %d/%d (low weight), wavelength = %.2f microns\n",
					iwave+1, paraminfo->num_waves, paraminfo->wavelength[iwave] );
				skip_psf = 1;
		   	}
		   	else
		   	{
	    	  		printf( "   Computing PSF %d/%d for wavelength %.2f microns (weight=%f)\n", 
		    			iwave+1, paraminfo->num_waves, paraminfo->wavelength[iwave], 
					paraminfo->weight[iwave] );
				skip_psf = 0;
		   	}
	        }
                else
		{
                   	printf("   Computing PSF for wavelength = %.2f microns\n", 
				paraminfo->wavelength[iwave]);
		   	skip_psf = 0;
		}

		if ( skip_psf == 0 )
		{
		  	if ( pupil_drawn == 0 )  /* compute pupil function */
		  	{
				temp = Alloc_image( crit_dim, crit_dim );

				/* compute aperture function */

				Compute_pupil( paraminfo, crit_dim, temp );
				for ( y = 0; y < crit_dim; ++y )
				{
					for ( x = 0; x < crit_dim; ++x )
					{
						pupil[y][x].r = temp[y][x];
						temp[y][x] = 0.0;
					}
				}

				/* compute wavefront error (OPD) function; wavefront error is *
				 * in microns.						      */

				Compute_opd( paraminfo, temp, crit_dim );

				for ( y = 0; y < crit_dim; ++y )
					for ( x = 0; x < crit_dim; ++x )
						pupil[y][x].i = temp[y][x];

				Free_image( temp );

				Shift_complex_to_origin( pupil, crit_dim );

				if ( paraminfo->num_waves > 1 )
				{
					Write_image( pupil_name, pupil[0], crit_dim, crit_dim, COMPLEX );
					file_written = 1;
				}

				pupil_drawn = 1;
		  	}
		  	else
		  	{
				Read_image( pupil_name, pupil[0], crit_dim, crit_dim, COMPLEX );
		  	}

		  	Compute_mono_psf( paraminfo, iwave, monopsf, pupil );
		  	Add_mono_to_poly( paraminfo, iwave, monopsf, polypsf );
	       }
	}

	if ( file_written )
		Delete_file( pupil_name );

	/* if IRAC channel 1 or 2, and PSF is not subsampled, then *
	 * convolve PSF with charge diffusion kernel.              */

	if ( paraminfo->camera == IRAC_35 || paraminfo->camera == IRAC_45 ) 
	{
		if ( !paraminfo->is_subsampled )
		{
			if ( paraminfo->use_kernel )
			{
				printf( "   Convolving PSF with IRAC charge diffusion kernel.\n" );
				Convolve_kernel( paraminfo, polypsf );
			}
		/*
			else
			{
				printf( "   NOT convolving PSF with IRAC charge diffusion kernel.\n" );
			}
		*/
		}
		/*
		else
		{
			printf( "   NOT convolving subsampled IRAC PSF with charge diffusion kernel.\n" );
		}
		*/
	}
	     
} /* Compute_poly_psf */

