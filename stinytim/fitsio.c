#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stinytim.h"

enum types { UNDEFINED, BIG_END, LIT_END };

/*-------------------------------------------------------------------------*/
int Which_byte_order( void )
{
        /* Determines byte order of system.  If  *
         * system does not use IEEE (eg. VAX),   *
         * then the byte order is undetermined.  */

        float a = 1.4335;
        unsigned char big_end[] = {238, 124, 183, 63};
        unsigned char lit_end[] = {63, 183, 124, 238}; 
        unsigned char *b;
        int type;
        int   i;

        b = (unsigned char *)&a;
        type = UNDEFINED;

        if ( b[0] == big_end[0] )
        {
                for ( i = 1; i <= 3; ++i )
		{
                        if ( b[i] == big_end[i] )
                                type = BIG_END;
			else
			{
				type = UNDEFINED;
				break;
			}
		}
        }
        else if ( b[0] == lit_end[0] )
        {
                for ( i = 1; i <= 3; ++i )
		{
                        if ( b[i] == lit_end[i] )
                                type = LIT_END;
			else
			{
				type = UNDEFINED;
				break;
			}
		}
        }

        return( type );
}

/*-------------------------------------------------------------------------*/

void Report_error( char *message, FILE *file )
{
	printf( message );
	printf( "\n" );
	if ( file != NULL )
		fclose( file );
	exit(0);
}

/*-------------------------------------------------------------------------*/
void Read_FITS_header( FILE *file, int *nx, int *ny )
{
	char	record[2880], *p;
	int	ival, line, nfound;

	/* NOTE : This routine is intended only to read in the header  *
	 * from the FITS files used by the Tiny Tim software.  It is   *
	 * not intended as a general FITS file header reader.          */

	if ( fread( record, 2880, 1, file ) != 1 )
		Report_error( "Error reading file header", file );

	nfound = 0;
	line = 0;
	p = record;

	while ( strncmp(p,"END ",4) != 0 ) 
	{
		if ( strncmp(p,"NAXIS",5) == 0 )
		{
			sscanf( p, "%*10c%d", &ival );
			if ( p[5] == ' ' )
			{
				if ( ival != 2 )
					Report_error("NAXIS not equal to 2", file);
				++nfound;
			}
			else if ( p[5] == '1' )
			{
				*nx = ival;
				++nfound;
			}
			else if ( p[5] == '2' )
			{
				*ny = ival;
				++nfound;
			}
		}
		else if ( strncmp( p, "BITPIX", 6 ) == 0 )
		{
			sscanf( p, "%*10c%d", &ival );
			if ( ival != -32 )
				Report_error( "BITPIX is not -32", file );
			++nfound;
		}
		p += 80;
		++line;

		if ( line == 36 )
		{
			if ( fread( record, 2880, 1, file ) != 1 )
				Report_error( "Error reading file header", file );
			line = 0;
			p = record;
		}
	}

	if ( nfound != 4 )
		Report_error( "Did not find a required header keyword", file );
}

/*---------------------------------------------------------------------------*/

void Swap_bytes( float *image, long numpix )
{
	char	*p, temp;
	long	i;

	for ( i = 0; i < numpix; ++i )
	{
		p = (char *)&image[i];

		temp = p[0];
		p[0] = p[3];
		p[3] = temp;
		temp = p[1];
		p[1] = p[2];
		p[2] = temp;
	}
} 

/*---------------------------------------------------------------------------*/
float **Read_FITS( char *filename, int *nx, int *ny )
{
	FILE 	*file;
	long	numpix;
	float	**image;
	int	order;

	/* NOTE : This routine is intended to read in only those FITS files  *
	 * provided with the Tiny Tim software.  It is not meant to be a     *
	 * general FITS reader.						     */

	order = Which_byte_order();
	if ( order == UNDEFINED )
	{
		printf("ERROR : Undefined machine byte order\n");
		exit(0);
	}

	file = fopen( filename, "rb" );
	if ( file == NULL )
		Report_error( "Could not open FITS file", file );

	Read_FITS_header( file, nx, ny );

	numpix = (long)*nx * (long)*ny;

	image = Alloc_image( *nx, *ny );

	if ( fread( image[0], sizeof(float), numpix, file ) != numpix )
		Report_error( "Could not read all requested data", file );

	if ( order == BIG_END )
		Swap_bytes( image[0], numpix );

	fclose( file );

	return( image );

} /* Read_FITS */  

/*---------------------------------------------------------------------------*/
void Write_FITS_header( FILE *file, int nx, int ny, paramstruct *paraminfo, int type )
{
	char	record[2880];
	int	camera, i;


	memset( record, ' ', 2880 );

	sprintf( &record[80*0], "SIMPLE  =                    T" );
	sprintf( &record[80*1], "BITPIX  = %20d", -32 );
	sprintf( &record[80*2], "NAXIS   = %20d", 2 );
	sprintf( &record[80*3], "NAXIS1  = %20d", nx );
	sprintf( &record[80*4], "NAXIS2  = %20d", ny );
	i = 4;

	if ( type == PSF_FILE )
	{
		camera = paraminfo->camera;
		sprintf( &record[80*(i+1)], "INSTRUME= '%18s' / Simulated instrument", paraminfo->camera_name );
		sprintf( &record[80*(i+2)], "FOCUS   = %20.4f / RMS focus in microns", paraminfo->focus );
		sprintf( &record[80*(i+3)], "X_COMA  = %20.4f / RMS X-coma in microns", paraminfo->xcoma );
		sprintf( &record[80*(i+4)], "Y_COMA  = %20.4f / RMS Y-coma in microns", paraminfo->ycoma );
		sprintf( &record[80*(i+5)], "X_ASTIG = %20.4f / RMS 0 degree astigmatism in microns", 
			paraminfo->xastig );
		sprintf( &record[80*(i+6)], "Y_ASTIG = %20.4f / RMS 45 degree astigmatism in microns", 
			paraminfo->yastig );
		sprintf( &record[80*(i+7)], "SPHERICL= %20.4f / RMS 3rd order spherical in microns", 
			paraminfo->spherical );
		sprintf( &record[80*(i+8)], "PIXSCALX= %20f / Arcseconds per array X pixel", 
			paraminfo->psf_pixel_size_arcsec_x );
		sprintf( &record[80*(i+9)], "PIXSCALY= %20f / Arcseconds per array Y pixel", 
			paraminfo->psf_pixel_size_arcsec_y );

		if ( camera != MIPS_70_SED && camera != MIPS_160 )
		{
		    sprintf( &record[80*(i+10)], "X_POS   = %20d / X position of PSF center in detector pixels",
			(int)paraminfo->x[paraminfo->current_pos] );
		    sprintf( &record[80*(i+11)], "Y_POS   = %20d / Y position of PSF center in detector pixels",
			(int)paraminfo->y[paraminfo->current_pos] );
		}
		else
		{
		    sprintf( &record[80*(i+10)], "X_POS   = %20.5f / X offset of PSF from camera center (arcmin)",
			paraminfo->x[paraminfo->current_pos] );
		    sprintf( &record[80*(i+11)], "Y_POS   = %20.5f / Y offset of PSF from camera center (arcmin)",
			paraminfo->y[paraminfo->current_pos] );
		}

		if ( paraminfo->num_waves > 1 )
			sprintf( &record[80*(i+12)],"SPECTRUM= '%18s' / Spectrum type or file",
				paraminfo->spectrum_file );
		else
			sprintf( &record[80*(i+12)],"WAVELNTH= %20f / Wavelength in microns",
				paraminfo->wavelength[0] );

		sprintf( &record[80*(i+13)],"XFIELD  = %20.5f / OTA X field position (arcmin)",
			paraminfo->xfield );
		sprintf( &record[80*(i+14)],"YFIELD  = %20.5f / OTA Y field position (arcmin)",
			paraminfo->yfield );

		if ( camera == MIPS_24 || camera == MIPS_70_WF || camera == MIPS_70_SUPER )
		{
			sprintf( &record[80*(i+15)],"SCANPOS = %20.6f / Scan field offset (arcsec)",
				paraminfo->scan_offset );
			sprintf( &record[80*(i+16)],"END" );
		}
		else
		{
			sprintf( &record[80*(i+15)],"END" );
		}
	}
	else
	{
		sprintf( &record[80*(i+1)], "END" );
	}

	/* go through record and replace \0 with spaces */

	for ( i = 0; i < 2880; ++i )
		if ( record[i] == '\0' )
			record[i] = ' ';
	
	if ( fwrite( record, 1, 2880, file ) != 2880 )
		Report_error( "Could not write new header", file );  

} /* Write_FITS_header */

/*---------------------------------------------------------------------------*/
void Write_FITS( char *filename, float *image, int nx, int ny, 
		 paramstruct *paraminfo, int type )
{
	/* WARNING : This routine may write over the input image with a *
	 *           byte-swapped version                               */

	FILE	*file;
	long	numpix, remain;
	char	record[2880];
	int	order;

	order = Which_byte_order();
	if ( order == UNDEFINED )
	{
		printf( "ERROR : Undefined byte order\n" );
		exit(0);
	}

	numpix = (long)nx * (long)ny;

	file = fopen( filename, "wb" );

	Write_FITS_header( file, nx, ny, paraminfo, type );

	if ( order == BIG_END )
		Swap_bytes( image, numpix );

	if ( fwrite( image, sizeof(float), numpix, file ) != numpix )
		Report_error( "Could not write data file", file );

	/* fill any unused bytes of the final 2880 byte record */

	remain = 2880 - (numpix * sizeof(float)) % 2880;

	if ( remain != 2880 )
	{
		memset( record, 0, 2880 );
		if ( fwrite( record, 1, remain, file ) != remain )
			Report_error( "Could not finish writing data", file );
	}

	fclose( file );

} /* WriteFits */
 
