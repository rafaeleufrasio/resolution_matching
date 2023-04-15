/*  
 *  File     :   psf.c
 *
 *  Contents :   
 *	Compute_nyquist_psf  :  Compute Nyquist-sampled PSF
 *
 *  Author   :   John Krist
 *
 *  Date     :   January 2000 (modified for SIRTF)
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "stinytim.h"


/*--------------------------------------------------------------------------
*  Compute_nyquist_psf :
*       Compute the Nyquist-sampled PSF from the complex pupil 
*       function "pupil" 
*
*  Inputs :
*       pupil  :  complex array of containing the pupil function
*       dim    :  dimension of "pupil" (dim by dim).
*
*  Outputs :
*       Replaces "pupil" with the Nyquist-sampled PSF.
*--------------------------------------------------------------------------*/
void Compute_nyquist_psf( complex **pupil, int dim )
{
	int     x, y;
	float   t, u, v, two_pi;

	two_pi = 2 * M_PI;

	/* Compute  Pupil.r * exp( 2*Pi*Func ) where Func = 0 + (Pupil.i)i. *
	 * The complex exponential is replaced with its sin, cos terms.     */

	for ( y = 0; y < dim; ++y )
		for ( x = 0; x < dim; ++x )
		{
			t = pupil[y][x].i * two_pi;
			pupil[y][x].i = sin(t) * pupil[y][x].r;
			pupil[y][x].r *= cos(t);
		}

	/* Compute the amplitude spread function (ASF) */

	fft2d( pupil[0], dim, 1 );

	/* Compute the PSF as the modulus squared of the ASF, place it  * 
	 * in the Pupil array.                                          */

	for ( y = 0; y < dim; ++y )
	{
		for ( x = 0; x < dim; ++x )
		{
			t = pupil[y][x].r;
			u = pupil[y][x].i;
			v = sqrt( t*t + u*u );
			pupil[y][x].r = v * v;
			pupil[y][x].i = 0.0;
		}
	}

} /* Compute_nyquist_psf */

