/* File      :  image.c
 * 
 * Contents  :
 *      Alloc_image    : Allocate a single precision image array
 *      Free_image     : Free an image allocated with Alloc_image
 *      Alloc_complex_image    : Allocate a complex valued image array
 *      Free_complex_image     : Free an image allocated with Alloc_complex_image
 *      Shift_to_center : Shift an image from 0,0 to the center
 *      Shift_complex_to_center : Shift a complex image from 0,0 to the center
 *      Shift_complex_to_origin : Shift a complex image from N/2,N/2 to the origin.
 *
 * Author    :  John Krist
 * Date      :  January 2000  (modified for SIRTF)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* #include <malloc.h> */
#include "stinytim.h"

/*--------------------------------------------------------------------------*/
float **Alloc_image( int nx, int ny )
{
	float   **image, *vector;
	int     i;


	image = (float **)malloc( ny * sizeof(float *) );
	if ( image == NULL )
	{
		printf( "ERROR (Alloc_image): Could not allocate array\n");
		exit(0);
	}

	vector = (float *)calloc( nx * ny, sizeof(float) );
	if ( vector == NULL )
	{
		printf( "ERROR (Alloc_image): Could not allocate vector\n");
		exit(0);
	}

	for ( i = 0; i < ny; ++i )
		image[i] = &vector[i*nx];

	return( image );

}  /* Alloc_image */

/*----------------------------------------------------------------------------
*  Free_image :
*       Free an image array allocated with Alloc_image.
*---------------------------------------------------------------------------*/
void Free_image( float **image )
{
	free( image[0] );
	free( image );     /* Free array of pointers */

} /* Free_image */

/*--------------------------------------------------------------------------*/
complex **Alloc_complex_image( int nx, int ny )
{
	complex **image, *vector;
	int     i;


	image = (complex **)malloc( ny * sizeof(complex *) );
	if ( image == NULL )
	{
		printf( "ERROR (Alloc_complex_image): Could not allocate array\n");
		exit(0);
	}

	vector = (complex *)calloc( nx * ny, sizeof(complex) );
	if ( vector == NULL )
	{
		printf( "ERROR (Alloc_complex_image): Could not allocate vector\n");
		exit(0);
	}

	for ( i = 0; i < ny; ++i )
		image[i] = &vector[i*nx];

	return( image );

}  /* Alloc_complex_image */

/*----------------------------------------------------------------------------
*  Free_complex_image :
*       Free a complex image array allocated with Alloc_complex_image.
*---------------------------------------------------------------------------*/
void Free_complex_image( complex **image )
{
	free( image[0] );
	free( image );     /* Free array of pointers */
}

/*-----------------------------------------------------------------------------
*  Shift_to_center :  
*       Shift an image of dimension (dim,dim) so that pixel (0,0)
*       is shifted to (dim/2,dim/2).
*
*       This routine works only on square, even dimensioned images!
*-----------------------------------------------------------------------------*/
void Shift_to_center( float **image, int dim )
{
	float   *temp;
	int     center, nhalf, ntot, y;


	temp = (float *)malloc( dim * sizeof(float) );
	if ( temp == NULL )
	{
		printf( "ERROR (Shift_to_center): Could not allocate memory\n");
		exit(0);
	}

	center = dim / 2;
	nhalf = center * sizeof(float);
	ntot = dim * sizeof(float);

	/* Shift in X direction */

	for ( y = 0; y < dim; ++y )
	{
		memcpy( &temp[center], image[y], nhalf );
		memcpy( temp, &image[y][center], nhalf );
		memcpy( image[y], temp, ntot );
	}

	/* Shift in Y direction */

	for ( y = 0; y < center; ++y )
	{
		memcpy( temp, image[y], ntot );
		memcpy( image[y], image[y+center], ntot );
		memcpy( image[y+center], temp, ntot );
	}
		
	free( temp );

} /* Shift_to_center */

/*-----------------------------------------------------------------------------
*  Shift_complex_to_center :  
*       Shift a complex valued image of dimension (dim,dim) so that pixel (0,0)
*       is shifted to (dim/2,dim/2).
*
*       This routine works only on square, even dimensioned images!
*-----------------------------------------------------------------------------*/
void Shift_complex_to_center( complex **image, int dim )
{
	complex *temp;
	int     center, nhalf, ntot, y;


	temp = (complex *)malloc( dim * sizeof(complex) );
	if ( temp == NULL )
	{
		printf( "ERROR (Shift_complex_to_center): Could not allocate memory\n");
	  	exit(0);
	}

	center = dim / 2;
	nhalf = center * sizeof(complex);
	ntot = dim * sizeof(complex);

	/* Shift in X direction */

	for ( y = 0; y < dim; ++y )
	{
		memcpy( &temp[center], image[y], nhalf );
		memcpy( temp, &image[y][center], nhalf );
		memcpy( image[y], temp, ntot );
	}

	/* Shift in Y direction */

	for ( y = 0; y < center; ++y )
	{
		memcpy( temp, image[y], ntot );
		memcpy( image[y], image[y+center], ntot );
		memcpy( image[y+center], temp, ntot );
	}
		
	free( temp );

} /* Shift_complex_to_center */

/*-----------------------------------------------------------------------------
*  Shift_complex_to_origin :  
*       Shift a complex valued image of dimension (dim,dim) so that pixel 
*	(dim/2,dim/2) is shifted to (0,0).
*
*       This routine works only on square, even dimensioned images!
*-----------------------------------------------------------------------------*/
void Shift_complex_to_origin( complex **image, int dim )
{

	Shift_complex_to_center( image, dim );

} /* Shift_complex_to_origin */

