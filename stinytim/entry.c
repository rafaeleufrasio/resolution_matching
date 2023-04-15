/*  File : entry.c
 *  Contents :
 *     InitEntryList   : Initialize counter for parameter buffer.
 *     DeleteEntryList : Deallocate the parameter buffer.
 *     WriteEntryList  : Write current parameter buffer to a file.
 *     StoreEntry      : Create an entry in the parameter buffer.
 *     GetEntry        : Read a line from a file and store it in the buffer.
 *
 *  Written by John Krist, June 1993.
 *
 *  The routines in this file are used to read in values from
 *  a Tiny Tim parameter file.  Comment lines (those begining
 *  with a "#") are skipped until an uncommented line is found.
 *  All lines, including commented lines, are stored in a buffer
 *  so that the entire parameter file can be written out at once.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
/* #include <malloc.h> */
#include "stinytim.h"

char    *Entry_list[1000];
int     Entry_number;

/*---------------------------------------------------------------------------*/
void Init_entry_list( void )
{
	Entry_number = 0;
}

/*---------------------------------------------------------------------------*/
void Delete_entry_list( void )
{
	int     i;

	for ( i = 0; i < Entry_number; ++i )
		free( Entry_list[i] );
}

/*---------------------------------------------------------------------------*/
void Write_entry_list( FILE *file )
{
	int     i;

	for ( i = 0; i < Entry_number; ++i )
		fputs( Entry_list[i], file );
}

/*---------------------------------------------------------------------------*/
void Store_entry( char *entry )
{
	Entry_list[Entry_number] = (char *)malloc( MAX_STRING );
	strcpy( Entry_list[Entry_number], entry );
	++Entry_number;
}

/*--------------------------------------------------------------------------
* Routine : GetEntry
* Purpose : Read a line from a table, ignoring comment (#) lines. Stores
*           all lines read into a list.
*--------------------------------------------------------------------------*/
const char *Get_entry( FILE *file )
{
	char    *n;
	static  char entry[MAX_STRING];

	if ( fgets( entry, MAX_STRING - 1, file ) == NULL )
	{
		printf( "ERROR (Get_entry) : Unexpected EOF\n");
		exit(0);
	}

	Store_entry( entry );

	while ( entry[0] == '#' )
	{
		if ( fgets( entry, MAX_STRING - 1, file ) == NULL )
		{
			printf( "ERROR (Get_entry) : Unexpected EOF\n");
			exit(0);
		}
		Store_entry( entry );
	}

	/* Remove any comments at the end of the line by placing a \0  *
	 * after the last parameter value on the line.                 */

	n = strchr( entry, '#' );
	if ( n != NULL )
	{
		--n;

		/* Work backwards to the last visible character */

		while ( !isgraph(*n) )
			*n-- = '\0';
	}

	return( entry );

} /* Get_entry */

