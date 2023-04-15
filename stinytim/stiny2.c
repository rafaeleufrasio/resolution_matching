#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stinytim.h"

#ifdef TT_THREADED
#include <unistd.h>
#endif

/*-----------------------------------------------------------------------------------------
*  Flip_image :
*	Flip the image horizontally or vertically.
*-----------------------------------------------------------------------------------------*/
void Flip_image( float **image, int n, int flip_x, int flip_y )
{
	int	x, y;
	float	temp;

	if ( flip_x != 0 )
	{
		for ( y = 0; y < n; ++y )
		{
			for ( x = 0; x < n/2; ++x ) 
			{
				temp = image[y][x];
				image[y][x] = image[y][n-x-1];
				image[y][n-x-1] = temp;
			}
		}
	}

	if ( flip_y != 0 )
	{
		for ( x = 0; x < n; ++x )
		{
			for ( y = 0; y < n/2; ++y )
			{
				temp = image[y][x];
				image[y][x] = image[n-y-1][x];
				image[n-y-1][x] = temp;
			}
		}
	}

}

/*-----------------------------------------------------------------------------------------
*  Set_field_position :
*	Convert from detector (x,y) pixel coordinates to telescope field (Z,Y)
*       coordinates (arcmin), accounting for detector rotation and scan offset.
*
* 	If MIPS 70 um SED or 160 um mode, then x & y are offsets along the detector
*	axes from the camera center in arcmin, with the scan offset implicitly included.
*----------------------------------------------------------------------------------------*/
void Set_field_position( paramstruct *paraminfo, float x, float y )
{
	float	arcmin_per_pix_x, arcmin_per_pix_y;
	float	xc, yc, xoff, yoff, t;
	int	camera;

	camera = paraminfo->camera;

	arcmin_per_pix_x = paraminfo->det_pixel_size_arcsec_x / 60.0;
	arcmin_per_pix_y = paraminfo->det_pixel_size_arcsec_y / 60.0;

	/* xoff, yoff = offset in arcmin from detector center in detector X-Y axis system */

	if ( camera != MIPS_70_SED && camera != MIPS_160 )
	{
		xoff = (x - paraminfo->x_pixels/2.0) * arcmin_per_pix_x;
		yoff = (y - paraminfo->y_pixels/2.0) * arcmin_per_pix_y;
	}
	else
	{
		xoff = x;
		yoff = y;
	}

	/* compute offset from detector center in CTA (Z,Y) system, *
	 * including detector rotation.				    */

	t = paraminfo->detector_rotation * M_PI / 180.0;
	paraminfo->xdelta = xoff * cos(t) - yoff * sin(t);
	paraminfo->ydelta = yoff * sin(t) + yoff * cos(t);

	xc = paraminfo->x_offset;   /* xc, yc = detector center in */ 
	yc = paraminfo->y_offset;   /* arcmin in CTA (Z,Y) system  */

	/* if MIPS 24 um or 70 um imaging modes, apply scan  *
	 * position offset to the detector center position   *
	 * (the scan position is implicit in the MIPS 70 um  *
	 * SED and 160 um modes).  The scan direction is at  *
	 * an angle for the 70 micron camera.                */

	if ( paraminfo->scan_offset != 0.0 )
	{
		/* scan_offset is in arcsec; convert to arcmin */

		t = paraminfo->scan_angle * M_PI / 180.0;
		xc += -paraminfo->scan_offset/60.0 * sin(t);
		yc +=  paraminfo->scan_offset/60.0 * cos(t);
	}

	/* compute position in CTA Z-Y field in arcmin */

	paraminfo->xfield = paraminfo->xdelta + xc;
	paraminfo->yfield = paraminfo->ydelta + yc;

} /* Set_field_position */

/*---------------------------------------------------------------------------------------*/
int main( int argc, char *argv[] )
{
	paramstruct paraminfo;
	float	**monopsf, **polypsf;
	complex **pupil;
	char	filename[MAX_STRING];
	int	i;

	if ( argc != 2 )
	{
		printf( "Call is : stiny2 paramfile\n" );
		exit(0);
	}

        NumThreads = 1;
#ifdef TT_THREADED
        NumThreads = sysconf(_SC_NPROCESSORS_ONLN);
#endif

	printf("Tiny Tim/Spitzer Version %s ", VERSION );
        if ( NumThreads > 1 )
                printf(" (Multithreaded, %d CPUs detected)", NumThreads);
	printf( "\n" );

	Read_parameters( argv[1], &paraminfo );

	monopsf = Alloc_image( paraminfo.psf_grid_size, paraminfo.psf_grid_size );
	polypsf = Alloc_image( paraminfo.psf_grid_size, paraminfo.psf_grid_size );
	pupil = Alloc_complex_image( paraminfo.nyquist_grid_size, paraminfo.nyquist_grid_size );

	for ( i = 0; i < paraminfo.num_positions; ++i )
	{
		paraminfo.current_pos = i;

		if ( paraminfo.camera != MIPS_70_SED && paraminfo.camera < IRS_LONG_LOW ) 
		{
			if ( paraminfo.num_positions > 1 )
			      printf("Computing PSF for (x,y) = (%g,%g)\n", paraminfo.x[i], paraminfo.y[i] );
			Set_field_position( &paraminfo, paraminfo.x[i], paraminfo.y[i] );
		}

		Compute_poly_psf( &paraminfo, pupil, polypsf, monopsf );

		Flip_image( polypsf, paraminfo.psf_grid_size, paraminfo.detector_x_flip, paraminfo.detector_y_flip );

		if ( paraminfo.num_positions > 1 )
			sprintf( filename, "%s%2.2d.fits", paraminfo.rootname, i+1 );
		else
			sprintf( filename, "%s.fits", paraminfo.rootname );

		printf("   Writing PSF to %s\n", filename );

		Write_FITS( filename, polypsf[0], paraminfo.psf_grid_size, paraminfo.psf_grid_size,
			&paraminfo, PSF_FILE );
	}

	Free_complex_image( pupil );
	Free_image( polypsf );
	Free_image( monopsf );

	return( 0 );
}

