/* File   :  misc.c
 *
 * Contents :
 *      Write_image :  Write a data array.
 *      Read_image  :  Read a data array created with Write_image.
 *
 * Author   :  John Krist
 * Date     :  January 2000 (modified for SIRTF)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
/* #include <malloc.h> */
#include "stinytim.h"


/*--------------------------------------------------------------------------
*  Strupcase :
*       Convert string to upper case.
*--------------------------------------------------------------------------*/
void Strupcase( char *string )
{
        while ( *string != '\0' )
        {
                *string = (char)toupper( *string );
                ++string;
        }
}

/*-------------------------------------------------------------------------
*  Find_string :
*       Find first occurrence of string2 in string1.  Returns pointer to
*       the location in string1 where string2 starts.  If string2 is
*       not found, a NULL pointer is returned.    Check is case insensitive
*-------------------------------------------------------------------------*/
char *Find_string( char *string1, char *string2 )
{
        char    *p0, *p1, *string1_up, *string2_up;

	/* convert input strings to upper case */

        string1_up = (char *)malloc(strlen(string1)+1);
        strcpy( string1_up, string1 );
        Strupcase( string1_up );
        string2_up = (char *)malloc(strlen(string2)+1);
        strcpy( string2_up, string2 );
        Strupcase( string2_up );

	/* look for string2 in string1 */

        p0 = string1_up;
        p1 = strstr( string1_up, string2_up );
        free( string1_up );
        free( string2_up );

        if ( p1 == NULL )
                return( NULL );			/* string2 not found */
        else
                return( &string1[p1-p0] );	/* string2 found */

} /* Find_string */

/*----------------------------------------------------------------------*/
void Delete_file( char *filename )
{
	unlink( filename );
}

/*----------------------------------------------------------------------*/
void Default_dir( char *filename )
{
        char    *temp;
        if ( (temp = getenv("STINYTIM")) == NULL )
        {
                printf( "ERROR : Environment variable STINYTIM undefined\n");
                exit(0);
        }

        strcpy( filename, temp );
	strcat( filename, "/" );
}

/*-----------------------------------------------------------------------
*  Write_image :
*       Write an image to disk.
*
*  Inputs :
*       filename : Name of file to write out to.
*       image    : Array which contains the data.
*       nx, ny   : Dimensions of image.
*      data_type : Either FLOATING or COMPLEX.
*-----------------------------------------------------------------------*/
void Write_image( char *filename, void *image, int nx, int ny, int data_type )
{
	int     size;
	FILE    *file;

	if ( data_type == FLOATING )
		size = sizeof(float);
	else
		size = sizeof(complex);

	if ( (file = fopen( filename, "wb" )) == NULL )
	{
		printf( "ERROR (Write_image) : Could not open file %s\n", filename );
		exit(0);
	}

	if ( fwrite( image, nx*ny*size, 1, file ) != 1 )
	{
		printf( "ERROR (Write_image) : Error writing to file.\n");
		fclose( file );
		exit(0);
	}

	fclose( file );

} /* Write_image */

/*-----------------------------------------------------------------------
*  Read_image :
*       Read an image written by Write_image into an array previously allocated
*       by Alloc_image or Alloc_complex_image.
*  Inputs :
*       filename : Name of file to read in.
*       image    : Array which will contain the data.
*       nx, ny   : Dimensions of image.
*      data_type : Either FLOATING or COMPLEX.
*-----------------------------------------------------------------------*/
void Read_image( char *filename, void *image, int nx, int ny, int data_type )
{
	FILE    *file;
	int     size;

	if ( (file = fopen( filename, "rb" )) == NULL )
	{
		printf( "ERROR (Read_image) : Could not open file %s\n", filename );
		exit(0);
	}

	if ( data_type == FLOATING )
		size = sizeof(float);
	else
		size = sizeof(complex);

	if ( fread( image, nx*ny*size, 1, file ) != 1 )
	{
		fclose( file );
		printf( "ERROR (Read_image) : Could not read file.\n" );
		exit(0);
	}

	fclose( file );

} /* Read_image */

