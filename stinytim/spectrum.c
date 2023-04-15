/* File     :  spectrum.c
 *
 * Contents :
 *	Integrate_spectrum
 *	Blackbody
 *	Power_law_nu
 *	Power_law_lambda
 *	Spectrum_file
 *
 * Author   :  John Krist
 * Date     :  January 2000 (modified for SIRTF)
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "stinytim.h"

#define MAX_POINTS  1000        /* Maximum number of points allowed in spectrum */
#define SQR(x) ((x)*(x))

#define DEFAULT_SPECTRUM 0
#define USER_SPECTRUM    1

#define FLUX_FLAM      1
#define FLUX_FNU       2
#define FLUX_JY        3
#define FLUX_PHOTLAM   4

/*-----------------------------------------------------------------------------------------
* Integrate_spectrum :
*	For each wavelength paraminfo->wavelength[i], compute the corresponding flux  
*	by interpolating the spectrum given in the arrays lambda and counts, then
*	integrate over the sampled wavelength span.  Note that the paraminfo->weight[] 
*	array must contain the system throughput curve before calling this routine.
*
* Inputs :
*	paraminfo : parameter info structure
*	   lambda : array containing spectrum wavelength values
*	   counts : array containing spectrum fluxes in photons/micron or directly
*		       proportional units
*         npoints : number of points in "lambda", "counts"
*-----------------------------------------------------------------------------------------*/
void Integrate_spectrum( paramstruct *paraminfo, float *lambda, float *counts, int npoints )
{
	int	i, j1, j2;
	float	dlambda, flux, slope, w, yint;

	/* lambda values must be in increasing order; spectrum fluxes must be in 
	 * photons/micron or directly proportional units */

	j1 = 0;
	j2 = 0;
	for ( i = 0; i < paraminfo->num_waves; ++i )
	{
		w = paraminfo->wavelength[i];

		/* find spectrum points bracketing wavelength "w" */

		while ( lambda[j2] < w && j2 < npoints-1 )    /* find endpoint */
			++j2;
		if ( j2 == 0 )
			j2 = 1;
		else
			j1 = j2 - 1;	/* beginning point */ 

		/* determine line between the two spectral points */

		slope = (counts[j2]-counts[j1]) / (lambda[j2]-lambda[j1]);
		yint = counts[j2] - slope * lambda[j2];

		/* crudely integrate spectrum between current and previous wavelengths */

		flux = slope * paraminfo->wavelength[i] + yint;
		if ( flux < 0.0 )
			flux = 0.0;
		if ( i == 0 )
		  	dlambda = fabs(paraminfo->wavelength[1] - paraminfo->wavelength[0]);
		else if ( i == paraminfo->num_waves-1 )
		  	dlambda = fabs(paraminfo->wavelength[i] - paraminfo->wavelength[i-1]);
		else
		  	dlambda = fabs((paraminfo->wavelength[i-1] - 
				        paraminfo->wavelength[i+1])/2.0 );
		paraminfo->weight[i] *= flux * dlambda;
	}

} /* Integrate_spectrum */

/*-----------------------------------------------------------------------------------
*  Blackbody :
*       Multiply the system throughput curve with a blackbody curve for a given
*	temperature.
*----------------------------------------------------------------------------------*/
void Blackbody( paramstruct *paraminfo )
{
	int   i;
	float dlambda, flux, h, c, k, t, w;

	t = 0.0;
	while ( t < 1.0 || t > 900000.0 )
	{	
		printf("\nEnter temperature (Kelvin) : ");
		scanf("%f", &t);
		if ( t < 1.0 || t > 900000.0 )
			printf("** Temperature must be 1 - 900000 K **\n"); 
	}

	h = 6.625e-27;
	c = 3.0e10;    		/* cm/sec */
	k = 1.38e-16;

	for ( i = 0; i < paraminfo->num_waves; ++i )
	{
		w = paraminfo->wavelength[i] * 1.0e-4;  /* convert um to cm */
		flux = 1.0 / (pow(w,4.0) * (exp((h*c)/(w*k*t))-1.0));

		/* flux is now in phot/cm^2/s/A with arbitrary normalization */

		if ( i == 0 )
		  	dlambda = fabs(paraminfo->wavelength[1] - 
				       paraminfo->wavelength[0]);
		else if ( i == paraminfo->num_waves-1 )
		  	dlambda = fabs(paraminfo->wavelength[i] - 
				       paraminfo->wavelength[i-1]);
		else
		  	dlambda = fabs((paraminfo->wavelength[i-1] - 
				        paraminfo->wavelength[i+1])/2.0 );
		paraminfo->weight[i] *= flux * dlambda;
	}

	sprintf( paraminfo->spectrum_file, "Blackbody(%gK)", t );

} /* Blackbody */

/*---------------------------------------------------------------------------
*  Power_law_nu :
*	Multiply the system throughput curve with a power law curve defined
*	as Flux(nu) = nu^alpha, where alpha is specified by the user.
*	The resulting fluxes are assumed to be directly proportional to
*	ergs/Hz (e.g. Janskys).  The normalization of the curve is not important.
*
* 	The array paraminfo->weight[] must contain the system throughput
*	curve prior to calling this routine.
*--------------------------------------------------------------------------*/
void Power_law_nu( paramstruct *paraminfo )
{
	double	alpha, c, dlambda, lambda, photlam, fnu, nu;
	int	i;


	printf( "\n\nYou have selected to use an object spectrum of the form : \n" );
	printf( "                   F(nu) = nu^alpha \n");
	printf( "where F(nu) is directly proportional to ergs/Hz.\n\n" );
	printf( "Enter the spectral index (alpha) : " );
	scanf( "%lf", &alpha );

	c = 3.0e14;	/* speed of light in microns */

	for ( i = 0; i < paraminfo->num_waves; ++i )
	{
		lambda = paraminfo->wavelength[i];
		nu = c / lambda;  /* convert microns to Hertz */
		fnu = pow(nu, alpha);

		/* "fnu" is currently in units ergs/Hz.  Convert to units directly      *
		 * proportional to photons/micron.  The normalization is not important. */

		photlam = fnu / lambda;
 
		if ( i == 0 )
		  	dlambda = fabs(paraminfo->wavelength[1] - 
				       paraminfo->wavelength[0]);
		else if ( i == paraminfo->num_waves-1 )
		  	dlambda = fabs(paraminfo->wavelength[i] - 
				       paraminfo->wavelength[i-1]);
		else
		  	dlambda = fabs((paraminfo->wavelength[i-1] - 
				        paraminfo->wavelength[i+1])/2.0 );
		paraminfo->weight[i] *= photlam * dlambda;
	}

	sprintf( paraminfo->spectrum_file, "f(nu)=nu^%g", alpha );

} /* Power_law_nu */

/*---------------------------------------------------------------------------
*  Power_law_lambda :
*	Multiply the system throughput curve with a power law curve defined
*	as Flux(lambda) = lambda^beta, where beta is specified by the user.
*	The resulting fluxes are assumed to be directly proportional to
*	ergs/micron.  The normalization of the curve is not important.
*
* 	The array paraminfo->weight[] must contain the system throughput
*	curve prior to calling this routine.
*--------------------------------------------------------------------------*/
void Power_law_lambda( paramstruct *paraminfo )
{
	float	beta, dlambda, photlam, flam;
	int	i;


	printf( "\n\nYou have selected to use an object spectrum of the form : \n" );
	printf( "               F(lambda) = lambda^beta \n" );
	printf( "where F(lambda) is directly proportional to ergs/micron.\n\n");
	printf( "Enter the spectral index (beta) : " );
	scanf( "%f", &beta );

	for ( i = 0; i < paraminfo->num_waves; ++i )
	{
		flam = pow(paraminfo->wavelength[i], beta);

		/* flam is in ergs/lambda; convert to units directly  *
		 * proportional to photons/micron.  The normalization *
		 * is not important.                                  */

		photlam = flam * paraminfo->wavelength[i];

		if ( i == 0 )
		  	dlambda = fabs(paraminfo->wavelength[1] - 
				       paraminfo->wavelength[0]);
		else if ( i == paraminfo->num_waves-1 )
		  	dlambda = fabs(paraminfo->wavelength[i] - 
				       paraminfo->wavelength[i-1]);
		else
		  	dlambda = fabs((paraminfo->wavelength[i-1] - 
				        paraminfo->wavelength[i+1])/2.0 );
		paraminfo->weight[i] *= photlam * dlambda;
	}

	sprintf( paraminfo->spectrum_file, "f(lambda)=lambda^%g", beta );

} /* Power_law_lambda */

/*-----------------------------------------------------------------------------
*  Spectrum_file :
*	Read in a user-supplied spectrum from the file specified by
*	paraminfo->spectrum_file.  The system throughput curve must be
*	in paraminfo->weight[] before calling this routine.  This routine
*	will interpolate the spectrum to find the flux at the wavelength
*	corresponding to weight[i] and integrate that value over the
*	wavelength span between weights, multiplying by the system
*	throughput.
*
*	The spectrum file must consist of wavelength & flux pairs, one pair
*	per line, up to MAX_POINTS.  The wavelengths must be in microns
*	and the fluxes in ergs/micron or directly proportional units.  
*	Comment lines in the file are denoted by a '#' as the first 
*	character in the first column.
*----------------------------------------------------------------------------*/
void Spectrum_file( paramstruct *paraminfo )
{
	FILE 	*file;
	char 	line[MAX_STRING];
	int  	i, flux_type;
	float 	lambda[MAX_POINTS], photlam[MAX_POINTS], flux;
	double	c, h, photon_energy;


	c = 3.0e14;	/* speed of light in microns/sec */
	h = 6.6266e-27; /* Planck's constant */

	file = NULL;

	while ( file == NULL )
	{
		file = fopen( paraminfo->spectrum_file, "r" );
		if ( file == NULL )
		{
			printf( "\nERROR : Could not open %s\n", paraminfo->spectrum_file );
			printf("\nEnter name of spectrum file : ");
			scanf("%s", paraminfo->spectrum_file);
		}
	}

	/* skip any comment lines (beginning with #) at beginning of the file */

	do {
		if ( fgets(line, MAX_STRING-1, file) == NULL )
		{
			printf("ERROR : Error reading %s\n", paraminfo->spectrum_file);
			fclose(file);
			exit(0);
		}
	} while ( line[0] == '#' );

	/* First non-comment line in file must be FLAM, JY, FNU, or PHOTLAM */

	if ( strstr(line,"FLAM") != NULL )   		/* ergs/micron */
		flux_type = FLUX_FLAM;
	else if ( strstr(line,"JY") != NULL )		/* Janskys */
		flux_type = FLUX_JY;
	else if ( strstr(line,"PHOTLAM") != NULL )	/* photons/micron */
		flux_type = FLUX_PHOTLAM;
	else if ( strstr(line,"FNU") != NULL )		/* ergs/Hz */
		flux_type = FLUX_FNU;
	else
	{
		printf( "ERROR : First non-comment line in spectrum file must be one\n");
		printf( "of the following : FNU, FLAM, JY, PHOTLAM\n");
		fclose(file);
		exit(0);
	}

	i = 0;
	while ( (fgets(line, MAX_STRING-1, file) != NULL) && (i < MAX_POINTS) )
	{
		/* must be wavelength(microns) & flux pairs */

		sscanf(line, "%f %f", &lambda[i], &flux);

		photon_energy = h * c / lambda[i];

		/* convert flux to photons/micron */

		if ( flux_type == FLUX_FLAM )
			photlam[i] = flux / photon_energy;
		else if ( flux_type == FLUX_JY )
			photlam[i] = 3.0e-9 * flux / SQR(lambda[i]) / photon_energy;
		else if ( flux_type == FLUX_FNU )
			photlam[i] = 3.0e14 * flux / SQR(lambda[i]) / photon_energy;
		else
			photlam[i] = flux;

		++i;
	}

	fclose( file );

	if ( i >= MAX_POINTS )
	{
	   printf( "\nERROR : Number of points in spectrum (%d) > limit (%d) ** \n",
			i, MAX_POINTS );
	   exit(0);
	}

	Integrate_spectrum( paraminfo, lambda, photlam, i );

} /* Spectrum_file */

