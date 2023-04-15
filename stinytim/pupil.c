#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "stinytim.h"

#define SUBFACTOR 8     /* subpixel resolution - DON'T CHANGE!!! */
#define	MAX_AREA	(SUBFACTOR*SUBFACTOR)
#define	MODRES(y)	((y) & 7)		/* subpixel Y modulo */
#define MAX_X	0x7FFF	/* subpixel X beyond right edge */
#define LERP(a,l,h)     ((l)+(((h)-(l))*(a)))
#define MIN(a,b)        (((a)<(b))?(a):(b))
#define MAX(a,b)        (((a)>(b))?(a):(b))
#define SQR(x)  ((x)*(x))

typedef struct VertexStruct 
{	
	int x, y;	/* subpixel coordinate */
} Vertex;

static Vertex *Vleft, *VnextLeft;	/* current left edge */
static Vertex *Vright, *VnextRight;	/* current right edge */

static struct SubPixel 
{
	/* subpixel extents for scanline */
	int xLeft, xRight;
} sp[SUBFACTOR];

static  int xLmin, xLmax;	/* subpixel x extremes for scanline */
static  int xRmax, xRmin;	/* (for optimization shortcut) */

/*---------------------------------------------------------------------------
* vLerp : Used by DrawPolygon routines
*--------------------------------------------------------------------------*/
void vLerp( double alpha, Vertex *Va, Vertex *Vb, Vertex *Vout )
{
	Vout->x = (Vb->x - Va->x) * alpha + Va->x;
	Vout->y = (Vb->y - Va->x) * alpha + Vb->y;
}

/*----------------------------------------------------------------------------
* Converge : used by DrawPolygon routines
*---------------------------------------------------------------------------*/
float Coverage( int x )
{
 	/* Compute number of subpixels covered by polygon at current pixel */
	/* x = left subpixel of pixel */

	float  area;			/* total covered area */
	int partialArea;	  	/* covered area for current subpixel y */
	int xr = x + SUBFACTOR - 1;	/*right subpixel of pixel */
	int y;

	/* shortcut for common case of fully covered pixel */

	if ( x > xLmax && x < xRmin )
		return( 0.0 );
	
	for ( area = y = 0; y < SUBFACTOR; y++ ) 
	{
		partialArea = MIN( sp[y].xRight, xr ) - MAX( sp[y].xLeft, x ) + 1;
		if ( partialArea > 0 )
			area += partialArea;
	}

	return( 1.0 - area / MAX_AREA );
}

/*------------------------------------------------------------------------------
*  RenderScanLine : used by DrawPolygon routines
*------------------------------------------------------------------------------*/
void RenderScanLine( int y, float **image, int nx, int ny )
{
	/* Vertex *Vl, *Vr = polygon vertices interpolated at scanline */
	/* y = scanline coordinate */
 
	int x;			/* leftmost subpixel of current pixel */

	if ( y < 0 || y >= ny )
		return;

	for ( x = SUBFACTOR*floor((double)(xLmin/SUBFACTOR)); x <= xRmax; x += SUBFACTOR )
		if ( x >= 0 && x < nx )  
			image[y][x/SUBFACTOR] *= Coverage(x);
}

/*-------------------------------------------------------------------------------
*  DrawPolygon :  Draw a filled, antialiased polygon.  Interior is filled with 
*                  zeroes.
*
*    Vertex polygon[] : array of numVertex Vertex structures containing
*                       polygon vertices at subpixel resolution elements
*       float **image : pointer to 2D array into which the polygon will be drawn
*          int nx, ny : dimensions of image
*
* Note : This code was adapted from the book Graphics Gems.  Code from that
*        book may be used without restriction.  The algorithm's author is
*        Jack C. Morrison.
*--------------------------------------------------------------------------------*/
void DrawPolygon( Vertex polygon[], int numVertex, float **image, int nx, int ny )
{
	Vertex *endPoly;			/* end of polygon vertex list */
	Vertex VscanLeft, VscanRight;		/* interpolated vertices at scanline */ 
	double aLeft, aRight;			/* interpolation ratios */
	struct SubPixel *sp_ptr;		/* current subpixel info */
	int xLeft, xNextLeft;			/* subpixel coordinates for */
	int  xRight, xNextRight;		/* active polygon edges */
	int i, y;

	xLeft = xRight = xNextLeft = xNextRight = 0;

	/* find vertex with minimum y (display coordinate) */

	Vleft = polygon;
	for  ( i = 1; i < numVertex; i++ )
		if ( polygon[i].y < Vleft->y )
			Vleft = &polygon[i];

	endPoly = &polygon[numVertex-1];

	/* initialize scanning edges */

	Vright = VnextRight = VnextLeft = Vleft;

	/* prepare bottom of initial scanline - no coverage by polygon */

	for ( i = 0; i < SUBFACTOR; i++ )
		sp[i].xLeft = sp[i].xRight = -1;

	xLmin = xRmin = MAX_X;
	xLmax = xRmax = -1;

	/* scan convert for each subpixel from bottom to top */

	for ( y = Vleft->y; ; y++ )
	{
		while ( y == VnextLeft->y )		/* reached next left vertex */
		{
			VnextLeft = (Vleft = VnextLeft) + 1; 	/* advance */
			if ( VnextLeft > endPoly )		/* (wraparound) */
				VnextLeft = polygon;
			if ( VnextLeft == Vright )		/* all y's same?  */
				return;				/* (null polygon) */ 
			xLeft = Vleft->x;
			xNextLeft = VnextLeft->x;
		}

		while ( y == VnextRight->y )    	/* reached next right vertex */
		{
			VnextRight = (Vright = VnextRight) - 1;
			if ( VnextRight < polygon )		/* (wraparound) */
				VnextRight = endPoly;
			xRight = Vright->x;
			xNextRight = VnextRight->x;
		}

		if ( y > VnextLeft->y || y > VnextRight->y )
		{
			/* done, mark uncovered part of last scanline */

			for (; MODRES(y); y++)
				sp[MODRES(y)].xLeft = sp[MODRES(y)].xRight = -1;
			RenderScanLine( y/SUBFACTOR, image, nx, ny );
			return;
		}

		/*
 	 	 * Interpolate sub-pixel x endpoints at this y,
 	 	 * and update extremes for pixel coherence optimization
 	 	 */
	
		sp_ptr = &sp[MODRES(y)];
		aLeft = (double)(y - Vleft->y) / (VnextLeft->y - Vleft->y);
		sp_ptr->xLeft = LERP( aLeft, xLeft, xNextLeft );
		if ( sp_ptr->xLeft < xLmin )
			xLmin = sp_ptr->xLeft;
		if ( sp_ptr->xLeft > xLmax )
			xLmax = sp_ptr->xLeft;

		aRight = (double)(y - Vright->y) / (VnextRight->y - Vright->y);
		sp_ptr->xRight = LERP( aRight, xRight, xNextRight );
		if ( sp_ptr->xRight < xRmin )
			xRmin = sp_ptr->xRight;
		if ( sp_ptr->xRight > xRmax )
			xRmax = sp_ptr->xRight;

		if ( MODRES(y) == SUBFACTOR-1 )	 	/* end of scanline */
		{
			/* interpolate edges to this scanline */

			vLerp( aLeft, Vleft, VnextLeft, &VscanLeft );
			vLerp( aRight, Vright, VnextRight, &VscanRight );
			RenderScanLine( y/SUBFACTOR, image, nx, ny );
			xLmin = xRmin = MAX_X; 		/* reset extremes */
			xLmax = xRmax = -1;
		}
  	}

} /* DrawPolygon */

/*------------------------------------------------------------------------------
*  DarkInsideCircle : Draw an antialiased circle filled with zeroes.
*
*   float **Image : 2D image into which circle is drawn
*      int Nx, Ny : dimensions of Image
*    float Xc, Yc : center of circle in Image
*    float Radius : radius of circle in Image
*------------------------------------------------------------------------------*/
void DarkInsideCircle( float **Image, int Nx, int Ny, 
			float Xc, float Yc, float Radius )
{
        float   *Row, RadiusMin2, ValueTable[200];
        float   xsub, ysub, xsub2, ysub2, SubInc, ysqr;
        int     x, y, x1, x2, y1, y2, SubPixelsSet;


        if ( Radius < 1.e-6 )
                return;

        RadiusMin2 = Radius - 2.0;
        SubInc = 1.0 / SUBFACTOR;

        /* Fill lookup table which converts number of subpixels set to  *
         * percentage of the pixel set.                                 */

        for ( x = 0; x <= SUBFACTOR*SUBFACTOR; ++x )
                ValueTable[x] = x * (SubInc*SubInc);

        y1 = (int)(Yc - Radius);
        y2 = (int)(Yc + Radius) + 1;
        if ( y1 < 0 )
                y1 = 0;
        if ( y2 >= Ny )
                y2 = Ny - 1;

        for ( y = y1; y <= y2; ++y )
        {
                Row = Image[y];
                if ( y <= Yc )
                        x2 = (int)(sqrt( SQR(Radius) - SQR(y+1-Yc) )) + 1;
                else
                        x2 = (int)(sqrt( SQR(Radius) - SQR(y-1-Yc) )) + 1;

                x1 = Xc - x2;
                x2 = x2 + Xc;

                if ( x1 < 0 )
                        x1 = 0;
                else if ( x2 >= Nx )
                        x2 = Nx - 1;

                ysqr = SQR(y - Yc);

                for ( x = x1; x <= x2; ++x )
                {
                     if ( sqrt( SQR(x-Xc) + ysqr ) > RadiusMin2 )
                     {
                          SubPixelsSet = 0;
                          for ( ysub = y; ysub < y+0.999; ysub += SubInc )
                          {
                                if ( ysub < Yc )
                                        ysub2 = SQR(ysub + SubInc - Yc);
                                else
                                        ysub2 = SQR(ysub - Yc);

                                for ( xsub = x; xsub < x+0.999; xsub+= SubInc )
                                {
                                        if ( xsub < Xc )
                                                xsub2 = SQR(xsub + SubInc - Xc);
                                        else
                                                xsub2 = SQR(xsub - Xc);

                                        if ( sqrt(xsub2 + ysub2) > Radius )
                                                ++SubPixelsSet;
                                }
                          }
                          Row[x] *= ValueTable[SubPixelsSet];
                     }
                     else
                          Row[x] = 0.0;
                }
        }

}  /* DarkInsideCircle */

/*------------------------------------------------------------------------------
*  DarkOutsideCircle : Fill an array with zeroes outside of an antialiased
*                      circle 
*
*   float **Image : 2D image into which circle is drawn
*      int Nx, Ny : dimensions of Image
*    float Xc, Yc : center of circle in Image
*    float Radius : radius of circle in Image
*------------------------------------------------------------------------------*/
void DarkOutsideCircle( float **Image, int Nx, int Ny, 
                        float Xc, float Yc, float Radius )
{
        float   xsub, ysub, xsub2, ysub2, SubInc;
        float   *Row, ValueTable[200];
        float   RadiusMin2, RadiusPlus2, ysqr, d;
        int     x, y, x1, x2, y1, y2, SubPixelsSet;


        SubInc = 1.0 / SUBFACTOR;   /* Subpixel size along one axis */

        /* Fill lookup table which converts number of subpixels set to  *
         * percentage of the pixel set.                                 */

        for ( x = 0; x <= SUBFACTOR*SUBFACTOR; ++x )
                ValueTable[x] = x * (SubInc*SubInc);

        RadiusMin2 = Radius - 2.0;
        RadiusPlus2 = Radius + 2.0;

        /* Determine top and bottom rows which contain the circle */

        y1 = (int)(Yc - Radius);
        y2 = (int)(Yc + Radius) + 1;
        if ( y1 < 0 )
                y1 = 0;
        if ( y2 >= Ny )
                y2 = Ny - 1;

        /* Fill up rows before circle with zeroes */

        for ( y = 0; y < y1; ++y )
                for ( x = 0; x < Nx; ++x )
                        Image[y][x] = 0.0;

        /* Fill up rows after circle with zeroes */

        for ( y = y2 + 1; y < Ny; ++y )
                for ( x = 0; x < Nx; ++x )
                        Image[y][x] = 0.0;

        for ( y = y1; y <= y2; ++y )
        {
                Row = Image[y];

                /* Determine the extent of the circle in this row */

                if ( y <= Yc )
                        x2 = (int)(sqrt( SQR(Radius) - SQR(y+1-Yc) )) + 1;
                else
                        x2 = (int)(sqrt( SQR(Radius) - SQR(y-1-Yc) )) + 1;

                x1 = Xc - x2;
                x2 = Xc + x2;

                if ( x1 < 0 )
                        x1 = 0;
                else if ( x2 >= Nx )
                        x2 = Nx - 1;

                /* Fill preceding and trailing columns with zeroes */

                for ( x = 0; x < x1; ++x )
                        Row[x] = 0.0;
                for ( x = x2 + 1; x < Nx; ++x )
                        Row[x] = 0.0;

                ysqr = SQR(y - Yc);

                for ( x = x1; x <= x2; ++x )
                {
                     /* If pixel is greater than Radius+2 from the center,  *
                      * the set it to 0.0; if it is less than Radius-2 from *
                      * the center, then don't change it; if it is in       *
                      * between, subpixelate the pixel to determine how     *
                      * much of the circle goes through it, and multiply    *
                      * the pixel value appropriately.                      */

                     d = sqrt( SQR(x-Xc) + ysqr );

                     if ( d > RadiusPlus2 )
                          Row[x] = 0.0;
                     else if ( d >= RadiusMin2 )
                     {
                          SubPixelsSet = 0;
                          for (ysub=y; ysub < y + 0.999;ysub += SubInc )
                          {
                                if ( ysub < Yc )
                                        ysub2 = SQR(ysub + SubInc - Yc);
                                else
                                        ysub2 = SQR(ysub - Yc);

                                for (xsub=x;xsub<x+0.999;xsub += SubInc)
                                {
                                        if ( xsub < Xc )
                                                xsub2 = SQR(xsub + SubInc - Xc);
                                        else
                                                xsub2 = SQR(xsub - Xc);

                                        if ( sqrt(xsub2 + ysub2) <= Radius )
                                                ++SubPixelsSet;
                                }
                          }
                          Row[x] *= ValueTable[SubPixelsSet];
                     }
                }
        }

} /* DarkOutsideCircle */

/*--------------------------------------------------------------------------------------
*  Compute_pupil :
*	Draw the SIRTF aperture pattern, including secondary supports and central
*	obscurations, for a given field angle.
*
*  Inputs :
*    paraminfo : pointer to parameter information structure
*            n : pupil array size (n by n pixels)
*        image : pupil array
*-------------------------------------------------------------------------------------*/
void Compute_pupil( paramstruct *paraminfo, int n, float **image )
{
	Vertex  p[4];

	/*  oldest values
	float	mbottom[] = {  0.0609555,   0.000224442, -0.00445731 };
	float   bbottom[] = {  0.00104665,  3.65148e-05, -0.000405502 };
	float	mtop[] =    { -0.0609555,  -0.000224442, -0.00445731 };
	float	btop[] =    { -0.00104665, -3.65148e-05, -0.000405502 };
	*/

	/* not-as-old old values
	float   mbottom[] = {  0.0348576,   7.12620e-05, -0.00170965 };
	float	bbottom[] = { -0.000816185, 2.30354e-05, -0.000290624 };
	float   mtop[] =    { -0.0334181,  -6.99303e-05, -0.00154696 };
	float	btop[] =    {  0.00202429, -4.89311e-06, -0.000219761 };
	*/

	float   mbottom[] = {  0.0372467,  8.32500e-05,  -0.00166610       };
	float	bbottom[] = { -0.00149612,  1.79408e-05, -0.000305957      };
	float   mtop[] =    { -0.0363804, -5.66100e-05,  -0.00151879       };
	float	btop[] =    { 0.00222571, -1.94761e-06, -0.000165466       };

	float	r, xc, yc, xv[4], yv[4];
	float	theta_field, theta_spider, xf, yf, theta_offset;
	float	slope_top, yint_top, slope_bottom, yint_bottom;
	float	topx1, topy1, topx2, topy2;
	float	bottomx1, bottomy1, bottomx2, bottomy2;
	int	i, j, x, y;

	xc = n / 2;
	yc = n / 2;
	r = n / 4.0;  /* outer aperture radius is n/4 providing Nyquist sampling */

	for ( y = 0; y < n; ++y )
		for ( x = 0; x < n; ++x )
			image[y][x] = 1.0;

	/* Draw outer aperture */

	DarkOutsideCircle( image, n, n, xc, yc, r );

	/* Draw central obscuration */

	DarkInsideCircle( image, n, n, xc, yc, r*paraminfo->secondary_radius );

	/* Draw three spider vanes */

	theta_offset = paraminfo->ota_rotation + paraminfo->detector_rotation;

	for ( i = 0; i <= 2; ++i )
	{
	        theta_field = (-120.0 + i * 120.0 + theta_offset) * M_PI / 180.0;
		if ( paraminfo->camera == IRAC_45 && i == 1 )
			theta_field = theta_field + 1.1 * M_PI / 180.0;
		else if ( (paraminfo->camera == IRAC_63 || paraminfo->camera == IRAC_35) && i == 0 )
			theta_field = theta_field + 1.1 * M_PI / 180.0;
        	theta_spider = -theta_field;

		/* paraminfo->xfield,yfield is the location of the PSF (arcmin) *
		 * in the CTA Z-Y system 				        */

        	xf = paraminfo->xfield*cos(theta_field) - paraminfo->yfield*sin(theta_field);
        	yf = paraminfo->xfield*sin(theta_field) + paraminfo->yfield*cos(theta_field);

        	slope_top = mtop[0] + mtop[1]*xf + mtop[2]*yf;
        	yint_top = btop[0] + btop[1]*xf + btop[2]*yf;
        	slope_bottom = mbottom[0] + mbottom[1]*xf + mbottom[2]*yf;
        	yint_bottom = bbottom[0] + bbottom[1]*xf + bbottom[2]*yf;

        	topx1 = -1.1;
        	topy1 = topx1 * slope_top + yint_top;
        	topx2 = -0.2;
        	topy2 = topx2 * slope_top + yint_top;

        	bottomx1 = -0.2;
        	bottomy1 = bottomx1 * slope_bottom + yint_bottom;
        	bottomx2 = -1.1;
        	bottomy2 = bottomx2 * slope_bottom + yint_bottom;

        	xv[0] = topx1*cos(theta_spider) - topy1*sin(theta_spider);
        	yv[0] = topx1*sin(theta_spider) + topy1*cos(theta_spider);
        	xv[1] = topx2*cos(theta_spider) - topy2*sin(theta_spider);
        	yv[1] = topx2*sin(theta_spider) + topy2*cos(theta_spider);

        	xv[2] = bottomx1*cos(theta_spider) - bottomy1*sin(theta_spider);
        	yv[2] = bottomx1*sin(theta_spider) + bottomy1*cos(theta_spider);
        	xv[3] = bottomx2*cos(theta_spider) - bottomy2*sin(theta_spider);
        	yv[3] = bottomx2*sin(theta_spider) + bottomy2*cos(theta_spider);

		for ( j = 0; j < 4; ++j )
		{
			p[j].x = (xv[j] * r + xc) * SUBFACTOR + SUBFACTOR/2;
			p[j].y = (yv[j] * r + yc) * SUBFACTOR + SUBFACTOR/2;
		}

		DrawPolygon( p, 4, image, n*SUBFACTOR, n*SUBFACTOR );
	}

} /* Compute_pupil */

