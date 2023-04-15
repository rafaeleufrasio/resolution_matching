#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "stinytim.h"

#ifdef TT_THREADED
#include <pthread.h>
#endif

static float **SincTable;

struct iparstruct {
	float mag;
        int   size_out;
        float **psf_out;
        int   size_in;
        complex **psf_in;
        float convergence;
	int   lower;
	int   upper;
	int   max_boxsize;
};

struct tparstruct {
        int start_row;
        struct iparstruct *ipars;
};

/*****************************************************************
  File : intpsf.c

  Contents :
       	     Sinc : sinc function = sin(x)/x
        Pixel_mtf : convolve PSF with pixel response function
         Bilinear : bilinear interpolation
  Convolve_jitter : convolve PSF with Gaussian jitter
         Sinc_int : sinc interpolation
    Integrate_psf : integrate Nyquist sampled PSF over detector pixels
  
  These routines convolve an input PSF image with the detector
  response function for an idealized square detector pixel and
  resamples the result.  The end product is a PSF integrated over
  detector pixels at detector pixel spacings.

  The input PSF is assumed to be critically sampled (or very near
  so).  Sinc interpolation is used for sampling.

  John Krist
  January 2000 (modified for SIRTF)
******************************************************************/ 

#define SQR(x) ((x)*(x))

/*---------------------------------------------------------------------------*/
float Sinc( float x )
{
	if ( x != 0.0 )
		return( sin(x) / x );
	else
		return( 1.0 );
} /* Sinc */

/*-----------------------------------------------------------------------
* Pixel_mtf :
*
*	Convolve the critically sampled PSF with the detector pixel response 
*	function.  This is equivalent to multiplying the optical transfer 
*       function (the Fourier transform of the PSF) with the detector pixel 
*       MTF in Fourier space.  Note that this assumes that the detector
*       pixel is perfectly square with uniform response over its entire area.
*
* Inputs :
*        otf :  Complex valued 2-d array containing Fourier transform of PSF.
*  psf_scale :  Scale of input PSF (units arbitrary; must be same units as
*		  "pixel_size"
*          n :  Dimension of otf array (n by n).
* pixel_size :  Pixel size in same units as psf_scale.
*
* Outputs :
*        Returns pixel-mtf-multipled otf in "otf".
*-----------------------------------------------------------------------*/
void Pixel_mtf( complex **otf, float psf_scale, int n, float pixel_size )
{
	float	p, inc, constant, *x, t;
	int	i, j;


	x = (float *)malloc( n * sizeof(float) );

	p = 0.5 * pixel_size / psf_scale;

	/* Generate a lookup table containing the pixel MTF */

	inc = 2.0 / n;
	constant = M_PI * p;

	x[0] = 0.0;
	x[n/2] = -1.;

	for ( i = 1; i < n/2; ++i )
	{
		x[i] = x[i-1] + inc;
		x[i+n/2] = x[i-1+n/2] + inc;
	}

	for ( i = 0; i < n; ++i )
		x[i] = Sinc( x[i] * constant );


	/* Convolve PSF with pixel response function     *
	 * (multiply OTF by pixel MTF in Fourier domain) */

	for ( j = 0; j < n; ++j )
	{
		for ( i = 0; i < n; ++i )
		{
			t = x[i] * x[j];
			otf[j][i].r *= t;
			otf[j][i].i *= t;
		}
	}

	free( x );

} /* Pixel_mtf */

/*------------------------------------------------------------------------------------
*  Convolve_jitter :
*	Convolve the Nyquist-sampled PSF with a Gaussian jitter function by
*	multiplying the OTF (FFT of PSF) by the jitter MTF in Fourier space.
*
*  Inputs :
*      psf_in : complex array containing OTF (FFT) of Nyquist-sampled PSF
*    scale_in : pixel scale of psf_in
*     size_in : dimensions of psf_in (size_in by size_in pixels)
*      jitter : one-sigma value of jitter in same units as scale_in (arcsec)
------------------------------------------------------------------------------------*/
void Convolve_jitter( complex **psf_in, float scale_in, int size_in, float jitter )
{
        float   constant, radius, dist, dist_sqr[4096], y_sqr, v;
        int     x, y;


        /* Create a lookup table of distances.  The table is shifted *
         * to account for the shifting of the OTF.  This essentially *
         * creates a jitter function shifted to (0,0).               */

        radius = size_in / 2;

        for ( x = 0; x < size_in/2; ++x )
        {
                dist = x / radius;
                dist_sqr[x] = dist * dist;
        }

        for ( x = size_in/2; x < size_in; ++x )
        {
                dist = (size_in - x) / radius;
                dist_sqr[x] = dist * dist;
        }

        /* Multiply the OTF (FFT of PSF) by the jitter MTF */

	constant = M_PI * jitter / (2 * scale_in);

        for ( y = 0; y < size_in; ++y )
        {
                y_sqr = dist_sqr[y];

                for ( x = 0; x < size_in; ++x )
                {
			v = sqrt(dist_sqr[x] + y_sqr) * constant;
                        v = exp( -2.0 * v * v );
                        psf_in[y][x].r *= v;
                        psf_in[y][x].i *= v;
                }
        }

} /* Convolve_jitter */

/*-----------------------------------------------------------------------------*/
void Int_row( int row, struct iparstruct *ipars )
{
	float	x_in, y_in, sum, oldsum, sum1, sum2;
	float	sinci1, sinci2, sincj1, sincj2, sincx, sincy;
	int	column, i, j, i0, j0, i1, i2, j1, j2, boxsize, ibox;

   	y_in = (row - ipars->size_out/2) * ipars->mag;
   	if ( y_in+ipars->size_in/2 < ipars->lower || y_in+ipars->size_in/2 > ipars->upper )
		return;
	
	for ( column = 0; column < ipars->size_out; ++column )
	{
		x_in = (column - ipars->size_out/2) * ipars->mag;
		if ( x_in+ipars->size_in/2 < ipars->lower || x_in+ipars->size_in/2 > ipars->upper )
			continue;

		if ( x_in < 0.0 )
			i0 = (int)(x_in - 0.5) + ipars->size_in/2;
		else
			i0 = (int)(x_in + 0.5) + ipars->size_in/2;

		if ( y_in < 0.0 )
			j0 = (int)(y_in - 0.5) + ipars->size_in/2;
		else
			j0 = (int)(y_in + 0.5) + ipars->size_in/2;

		i1 = i2 = i0;
		j1 = j2 = j0;

                sum = 0.0;
		oldsum = ipars->psf_in[j1][i1].r * SincTable[column][i1] * SincTable[row][j1];
		boxsize = 3;

                do {
                       	oldsum += sum;
                       	sum = 0.0;

		  	for ( ibox = 1; ibox <= 3; ++ibox )
		  	{
                        	--i1;
                             	++i2;
                             	--j1;
                             	++j2;

                             	sinci1 = SincTable[column][i1];
                             	sinci2 = SincTable[column][i2];
                             	sincj1 = SincTable[row][j1];
                             	sincj2 = SincTable[row][j2];

                             	/* Do next bottom and top rows */

                             	sum1 = 0.0;
                             	sum2 = 0.0;

                             	for ( i = i1; i <= i2; ++i )
                             	{
                               		sincx = SincTable[column][i];
                               		sum1 += ipars->psf_in[j1][i].r * sincx;
                               		sum2 += ipars->psf_in[j2][i].r * sincx;
                             	}

                             	sum += sum1 * sincj1 + sum2 * sincj2;

                             	/* Do next left and right columns */

			     	sum1 = 0.0;
                             	sum2 = 0.0;

                             	for ( j = j1 + 1; j < j2; ++j )
                             	{
                               		sincy = SincTable[row][j];
                               		sum1 += ipars->psf_in[j][i1].r * sincy;
                               		sum2 += ipars->psf_in[j][i2].r * sincy;
                             	}

                             	sum += sum1 * sinci1 + sum2 * sinci2;

			     	boxsize = boxsize + 2;

                          	if ( boxsize > ipars->max_boxsize )
                               		break;
		  	}
               		if ( boxsize > ipars->max_boxsize )
                      		break;
                } while ( fabs(sum/oldsum) > ipars->convergence );

                oldsum += sum;

                if ( oldsum > 0.0 )
                        ipars->psf_out[row][column] = oldsum;
                else
			ipars->psf_out[row][column] = 0.0;
	}
}

#ifdef TT_THREADED
/*---------------------------------------------------------------------------
* Integrate_thread :
*---------------------------------------------------------------------------*/
void *Integrate_thread( void *vpars )
{
        int	row;
        struct  tparstruct *tpars = vpars;

        for ( row = tpars->start_row; row < tpars->ipars->size_out; row += NumThreads )
                Int_row( row, tpars->ipars );

	return( NULL );
}
#endif

/*-----------------------------------------------------------------------------
*  Sinc_int :
*
*    Sample a 2D array ("psf_in") using adaptive sinc interpolation.
*
*  Inputs :
*	psf_in : 2D array containing PSF to be sampled (the critically sampled
*		    PSF convolved with the detector response function);
*		    This is a complex-valued array, but only the real portion
*		    is interpolated.
*     scale_in : pixel scale of "psf_in"; units are arbitrary but must be
*		    same as for "scale_out"
*      size_in : array dimension of "psf_in" (size_in by size_in pixels)
*      psf_out : 2D array to contain sampled PSF (array must have already
*		    been allocated by caller)
*    scale_out : pixel scale of "psf_out" in same units as "scale_in"
*     size_out : array dimension of "psf_out"
*  convergence : sinc interpolation summation convergence limit
*  max_boxsize : sinc interpolation summation boxsize (max)
*---------------------------------------------------------------------------*/ 
void Sinc_int( complex **psf_in, float scale_in, int size_in,
	       float **psf_out, float scale_out, int size_out,
	       float convergence, int max_boxsize, float damping )
{
#ifdef TT_THREADED
        pthread_t thread[MAX_THREADS];
	struct  tparstruct tpars[MAX_THREADS];
	int	ithread;
	void	*retval;
#endif
	int	pix_in, pix_out, y;
	float	pix, r, t;
	struct  iparstruct ipars;


	ipars.mag = scale_out / scale_in;
        ipars.size_out = size_out;
        ipars.psf_out = psf_out;
        ipars.size_in = size_in;
        ipars.psf_in = psf_in;
        ipars.convergence = convergence;
	ipars.lower = max_boxsize / 2 + 2;
	ipars.upper = size_in - max_boxsize / 2 - 2;
	ipars.max_boxsize = max_boxsize;

	/* create lookup table of sinc interpolator coefficients */

        SincTable = Alloc_image( size_in, size_out );

        for ( pix_out = 0; pix_out < size_out; ++pix_out )
        {
                for ( pix_in = 0; pix_in < size_in; ++pix_in )
                {
                        pix = (pix_out - size_out/2) * ipars.mag + size_in/2;
                        if ( pix != pix_in )
                        {
				r = pix_in - pix;
                                t = M_PI * r;
                                SincTable[pix_out][pix_in] = exp(-0.5*r*r/damping) * sin(t) / t;
                        }
                        else
                        {
                                SincTable[pix_out][pix_in] = 1.0;
                        }
                }
        }


#ifndef TT_THREADED
	for ( y = 0; y < size_out; ++y )
		Int_row( y, &ipars );
#else
        for ( ithread = 0; ithread < NumThreads; ++ithread )
        {
                tpars[ithread].start_row = ithread;
                tpars[ithread].ipars = &ipars;
                if ( pthread_create(&thread[ithread], NULL, Integrate_thread, &tpars[ithread]) )
                {
                        fprintf(stderr, "Cannot make thread %d\n", ithread);
                        exit(1);
                }
        }

        /* Join (collapse) the two threads */

        for ( ithread = 0; ithread < NumThreads; ++ithread )
        {
                if ( pthread_join(thread[ithread], &retval) )
                {
                        fprintf(stderr, "Thread join failed\n");
                        exit(1);
                }
        }
#endif

	Free_image( SincTable );

} /* Sinc_int */

/*----------------------------------------------------------------------------
*  Integrate_psf :
*
*     Integrate the Nyquist (critically) sampled PSF "psf_in" (with a sampling
*     of "scale_in") onto "scale_out" sized detector pixels.  This is done
*     by convolving the critically sampled PSF with the pixel response
*     function and then sampling the result at detector-pixel spacings.
*
*  Inputs :
*      psf_in : Complex valued 2D array containing critically sampled PSF;
*		   (only the real portion of this data is used);
*    scale_in : Spacing of "psf_in" values in same units as "scale_out"
*     size_in : Dimension of "psf_in" (size_in by size_in)
*     psf_out : 2D array allocated by the caller into which the final
*	  	   PSF is placed
*   scale_out : Spacing of "psf_out" values (detector pixel size) in same 
*		   units as "scale_in"
*    size_out : Dimension of "psf_out"
*      jitter : If jitter>0, the amount of symmetrical Gaussian 2D jitter
*	           to apply to the image; value is RMS jitter in same units
*		   as "scale_in" (arcsec)
* convergence : Sinc interpolation summation convergence limit.
* max_boxsize : Sinc interpolation summation boxsize (max) 
*
* Returns :
*    Places integrated PSF in "psf_out".
*---------------------------------------------------------------------------*/ 
void Integrate_psf( complex **psf_in, float scale_in, int size_in,
		   float **psf_out, float scale_out, int size_out,
		   float jitter, float convergence, int max_boxsize,
		   float damping )
{
	float	constant, i, r;
	int	x, y;

	/* compute Fourier transform of PSF (=OTF) */

        fft2d( psf_in[0], size_in, -1 );

	/* multiply OTF with pixel MTF */

	Pixel_mtf( psf_in, scale_in, size_in, scale_out );

	/* if jitter specified, convolve PSF with jitter  *
	 * (multiply PSF by jitter MTF in Fourier space)  */

	if ( jitter > 0.00001 )	/* jitter is given in arcsec */
		Convolve_jitter( psf_in, scale_in, size_in, jitter );

	/* Transform MTF back into PSF */

	fft2d( psf_in[0], size_in, 1 );

        constant = SQR(scale_out/scale_in) / SQR((float)size_in);

        for ( y = 0; y < size_in; ++y )
        {
                for ( x = 0; x < size_in; ++x )
                {
                        r = psf_in[y][x].r;
                        i = psf_in[y][x].i;

			/* ensure PSF is real-valued; FFT noise can often *
			 * creep into the imaginary portion               */

                        psf_in[y][x].r = sqrt(r*r + i*i) * constant;
			psf_in[y][x].i = 0.0;
                }
        }

        Shift_complex_to_center( psf_in, size_in );

	/* Sample convolved PSF using sinc interpolation */

	Sinc_int( psf_in, scale_in, size_in, psf_out, scale_out, size_out, 
		convergence, max_boxsize, damping );

} /* Integrate_psf */

