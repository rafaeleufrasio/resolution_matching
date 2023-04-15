#include <stdio.h>
#include <math.h>
#include "stinytim.h"

#define SQR(x) ((x)*(x))

/*--------------------------------------------------------------------------
*  Field_aberration :
*
*  Compute the specified, field dependent aberration for a given
*  position on the detector. 
*
*  Inputs :
*	xpsf, ypsf : offset in arcmin from the detector center, in the Z-Y
*		       CTA coordinate system (NOT offset from CTA center!)
*               ab : aberration index
*---------------------------------------------------------------------------*/
float Field_aberration( float xpsf, float ypsf, int ab, paramstruct *paraminfo )
{
	float	sum;
	int	i, j;

	sum = 0.0;

	for ( j = 0; j <= 3; ++j )
		for ( i = 0; i <= 3; ++i )
			sum += paraminfo->caber[ab][i][j] * 
				pow(xpsf,(float)i) * pow(ypsf,(float)j);
	return( sum );
}

/*--------------------------------------------------------------------------
*  MIPS_field_aberration :
*
*  Compute the specified, field dependent aberration for a given
*  position on a MIPS detector.  For the MIPS 70 um SED and 160 um modes,
*  Field_aberration is used instead, since the provided coordinates in 
*  those modes implicitly include the scan offset.
*
*  Inputs :
*	xpsf, ypsf : offset in arcmin from the detector center, in the Z-Y
*		       CTA coordinate system (NOT offset from CTA center!)
*       scan_index : index of reference scan position 
*               ab : aberration index
*---------------------------------------------------------------------------*/
float MIPS_field_aberration( float xpsf, float ypsf, int scan_index, 
			     int ab, paramstruct *paraminfo )
{
	float	sum;
	int	i, j;

	sum = 0.0;

	for ( j = 0; j <= 3; ++j )
		for ( i = 0; i <= 3; ++i )
			sum += paraminfo->mipscaber[scan_index][ab][i][j] * 
				pow(xpsf,(float)i) * pow(ypsf,(float)j);
	return( sum );
}

/*--------------------------------------------------------------------------
*  MIPS_aberration :
*
*  Compute the specified, field dependent aberration for a given
*  position on the detector and scan offset.
*
*  Inputs :
*	xpsf, ypsf : offset in arcmin from the detector center, in the Z-Y
*		       CTA coordinate system (NOT offset from CTA center!)
*               ab : aberration index
*---------------------------------------------------------------------------*/
float MIPS_aberration( float xpsf, float ypsf, int ab, paramstruct *paraminfo )
{
	int	s0, s1;
	float	slope, yint, center_ab, ab_offset, ab0, ab1;


	/* mipsscan is arranged from - to + offset */

	if ( paraminfo->scan_offset < paraminfo->mipsscan[1] )
	{
		s0 = 0;
		s1 = 1;
	}
	else
	{
		s0 = 1;
		s1 = 2;
	}

	/* compute aberration at detector center at specified scan position */

	slope = (paraminfo->mipsaber[s1][ab] - paraminfo->mipsaber[s0][ab]) /
		(paraminfo->mipsscan[s1] - paraminfo->mipsscan[s0]);
	yint = paraminfo->mipsaber[s0][ab] - slope * paraminfo->mipsscan[s0];
	center_ab = paraminfo->scan_offset * slope + yint;

	/* compute aberration delta at given field position and scan position */

	ab0 = MIPS_field_aberration( xpsf, ypsf, s0, ab, paraminfo );
	ab1 = MIPS_field_aberration( xpsf, ypsf, s1, ab, paraminfo );
	slope = (ab1 - ab0) / (paraminfo->mipsscan[s1] - paraminfo->mipsscan[s0]);
	yint = ab0 - slope * paraminfo->mipsscan[s0];
	ab_offset = paraminfo->scan_offset * slope + yint;

	return( center_ab + ab_offset );

} /* MIPS_aberration */

/*--------------------------------------------------------------------------
*  Compute_opd  :
*
*    Compute a 2D image which contains the optical path error (in microns RMS)
*    for each point in the pupil.  This is computed from the sum of 33% obscured
*    Zernike polynomials.  Because the field aberrations do not change 
*    significantly with scan position in MIPS, the field position used by
*    this routine is what it would be if the scan position were zero.  This
*    is required because the field-dependent aberration coefficients are
*    relative to the OTA field position.
* 
*    Since the pupil is defined only within the central n/2 by n/2 portion
*    of the image, only those points within that region are computed.
*
*    Inputs:
*        paraminfo : parameter information structure
*              opd : image array allocated with Alloc_image of size n by n.
*                n : dimension of opd array (n by n).
*
*   Returns:
*       Returns the OPD (optical path difference) function in "opd".
*-----------------------------------------------------------------------------*/
void Compute_opd( paramstruct *paraminfo, float **opd, int n )
{
	float   center, radius, dist[4096], sqr[4096], *zval, psfx, psfy;
	float   r, r2, r3, r4, r5, r6, theta, value, ydist, ysqr;
	float   sin_theta, sin_theta2, sin_theta3, sin_theta4, sin_theta5;
	float   cos_theta, cos_theta2, cos_theta3, cos_theta4, cos_theta5;
	int     camera, x, y, x1, x2, y1, y2;

	camera = paraminfo->camera;

	zval = paraminfo->zval;

	/* Assign aberrations for field center (for those aberrations   *
	 * that might vary with position); aberrations are specified as *
         * 33% obscured Zernike polynomial coefficients in microns RMS. *
	 * The angle of the aberrations is defined from the +Z axis of  *
	 * the telescope to the +Y axis.  The aberration function thus  *
	 * must be rotated to match the assumed detector orientation.   */

	/* these aberrations might vary with field */

	paraminfo->focus = zval[4];
	paraminfo->xastig = zval[5];
	paraminfo->yastig = zval[6];
	paraminfo->xcoma = zval[7];
	paraminfo->ycoma = zval[8];
	paraminfo->spherical = zval[11];

	/* Adjust certain aberrations depending on field position */

	if ( paraminfo->adjust_field_aberrations != 0 && 
	     camera != MIPS_70_SED && camera < IRS_LONG_LOW )
	{
	   /* psfx, psfy are offsets from the detector center in arcmin in  *
	    * the CTA Z-Y system (NOT offset from CTA center)               */

	   psfx = paraminfo->xdelta;
	   psfy = paraminfo->ydelta;

	   if ( camera == MIPS_24 || camera == MIPS_70_WF || camera == MIPS_70_SUPER )
	   {
		paraminfo->focus  += MIPS_aberration( psfx, psfy, 0, paraminfo );
		paraminfo->xastig += MIPS_aberration( psfx, psfy, 1, paraminfo );
		paraminfo->yastig += MIPS_aberration( psfx, psfy, 2, paraminfo );
		paraminfo->xcoma  += MIPS_aberration( psfx, psfy, 3, paraminfo );
		paraminfo->ycoma  += MIPS_aberration( psfx, psfy, 4, paraminfo );
		paraminfo->spherical += MIPS_aberration( psfx, psfy, 5, paraminfo );
	   }
	   else
	   {
		paraminfo->focus  += Field_aberration( psfx, psfy, 0, paraminfo );
		paraminfo->xastig += Field_aberration( psfx, psfy, 1, paraminfo );
		paraminfo->yastig += Field_aberration( psfx, psfy, 2, paraminfo );
		paraminfo->xcoma  += Field_aberration( psfx, psfy, 3, paraminfo );
		paraminfo->ycoma  += Field_aberration( psfx, psfy, 4, paraminfo );
		paraminfo->spherical += Field_aberration( psfx, psfy, 5, paraminfo );
	   }
	}

	center = n / 2;
	radius = n / 4.0;  /* aperture radius set for Nyquist sampling */

	/* Since the aperture is defined only within the central *
         * n/2 by n/2 region of the image, only bother with that *
         * part.                                                 */

	x1 = y1 = (int)(center - radius) - 6;
	x2 = y2 = (int)(center + radius) + 6;

	/* Compute lookup table of radii */

	for ( x = 0; x < n; ++x )
	{
		dist[x] = (x - center) / radius;
		sqr[x] = SQR(dist[x]);
	}

	for ( y = y1; y <= y2; ++y )
	{
		ydist = dist[y];
		ysqr = sqr[y];

		for ( x = x1; x <= x2; ++x )
		{
			value = 0.0;

			r = sqrt( sqr[x] + ysqr );
			r2 = r * r;
			r3 = r2 * r;
			r4 = r3 * r;
			r5 = r4 * r;
			r6 = r5 * r;

			if ( ydist != 0.0 || dist[x] != 0.0 )
				theta = atan2( ydist, dist[x] ) + 
					  M_PI * paraminfo->detector_rotation / 180.0;
			else
				theta = 0.0;

			sin_theta = sin(theta);
			cos_theta = cos(theta);
			sin_theta2 = sin(theta*2.);
			cos_theta2 = cos(theta*2.);
			sin_theta3 = sin(theta*3.);
			cos_theta3 = cos(theta*3.); 
			sin_theta4 = sin(theta*4.);
			cos_theta4 = cos(theta*4.);
			sin_theta5 = sin(theta*5.);
			cos_theta5 = cos(theta*5.);

			if ( zval[2] )    /* X tilt */
				value += zval[2] * 1.8992573 * r * cos_theta;
			if ( zval[3] )    /* Y tilt */
				value += zval[3] * 1.8992573 * r * sin_theta;
			if ( zval[4] )    /* Focus */
				value += paraminfo->focus * 3.8874443 * (r2 - 0.554450);
			if ( zval[5] )    /* 0 degree astigmatism */
				value += paraminfo->xastig * 2.3137662 * r2 * cos_theta2;
			if ( zval[6] )    /* 45 degree astigmatism */
				value += paraminfo->yastig * 2.3137662 * r2 * sin_theta2;
			if ( zval[7] )    /* X coma */
				value += paraminfo->xcoma * 8.3345629 * (r3 - 0.673796 * r) * cos_theta;
			if ( zval[8] )    /* Y coma */
				value += paraminfo->ycoma * 8.3345629 * (r3 - 0.673796 * r) * sin_theta;
			if ( zval[9] )    /* X clover */
				value += zval[9] * 2.6701691 * r3 * cos_theta3;
			if ( zval[10] )   /* Y clover */
				value += zval[10] * 2.6701691 * r3 * sin_theta3;
			if ( zval[11] )   /* 3rd order spherical */
				value += paraminfo->spherical * 16.895979 * (r4 - 1.1089 * r2 + 0.241243);
			if ( zval[12] )   /* Spherical astigmatism */
				value += zval[12] * 12.033645 * (r4 - 0.750864 * r2) * cos_theta2;
			if ( zval[13] )   /* 45 degree Spherical astigmatism */
				value += zval[13] * 12.033645 * (r4 - 0.750864 * r2) * sin_theta2;
			if ( zval[14] )   /* Ashtray */
				value += zval[14] * 2.9851527 * r4 * cos_theta4;
			if ( zval[15] )   /* Ashtray */
				value += zval[15] * 2.9851527 * r4 * sin_theta4;
			if ( zval[16] )
				value += zval[16] * 36.321412 * (r5 - 1.230566 
					* r3 + 0.323221 * r) * cos_theta;
			if ( zval[17] )
				value += zval[17] * 36.321412 * (r5 - 1.230566 
					* r3 + 0.323221 * r) * sin_theta;
			if ( zval[18] )
				value += zval[18] * 16.372202 * (r5 - 0.8001 * r3) * cos_theta3;
			if ( zval[19] )
				value += zval[19] * 16.372202 * (r5 - 0.8001 * r3) * sin_theta3;
			if ( zval[20] )
				value += zval[20] * 3.2700486 * r5 * cos_theta5;
			if ( zval[21] )
				value += zval[21] * 3.2700486 * r5 * sin_theta5;
			if ( zval[22] )   /* 5th order spherical */
				value += zval[22] * 74.82446 * (r6 - 1.663350 
					* r4 + 0.803136 * r2 - 0.104406);

			opd[y][x] = value;
		}
	}

	if ( paraminfo->use_map == 1 )
		Mirror_map( opd, paraminfo );

} /* Compute_opd */

