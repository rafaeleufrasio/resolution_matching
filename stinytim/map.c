#include <stdio.h>
#include <string.h>
#include <math.h>
#include "stinytim.h"


/****************************************************************************
*  Primary_map :
*	Read in map of surface errors on the primary mirror; scale it to
*	match the pupil diameter
*
*  Inputs :
*	map : 2D array into which map is stored; should be previously
*		allocated
*	paraminfo : pointer to parameter information structure
****************************************************************************/ 
void Primary_map( float **map, paramstruct *paraminfo )
{
	char	filename[MAX_STRING];
	float	**primary, delta, r_in, r_out, x_in, y_in;
	float	scale, t, u, mirror_radius;
	int	n, nx, ny, x1, x2, y1, y2, x_out, y_out, xc_in, yc_in;
	int	xc_out, yc_out;

	n = paraminfo->nyquist_grid_size;

	/* read in primary mirror map;  map is oriented so that the *
	 * +X array axis is the +Z telescope axis and the +Y array  *
	 * axis is the +Y telescope axis;  the map is in microns    *
	 * of wavefront error (twice the surface error).            */

	Default_dir( filename );
	strcat( filename, "sirtf_primary.fits" );
	primary = Read_FITS( filename, &nx, &ny ); 

	delta = 1.92735;	/* spacing in mm of map in FITS file */
	mirror_radius = 425.0;  /* mirror (pupil) radius in mm */
	r_in = mirror_radius / delta;  /* pupil radius in map pixels */
	xc_in = 245;  /* pupil center in map pixels */
	yc_in = 244;

	xc_out = n / 2;   /* pupil center in aberration array */
	yc_out = n / 2;
	r_out = n / 4.0;  /* pupil radius in aberration array pixels */
	scale = r_in / r_out;

	/* scale the map using interpolation to match the pupil size */

	for ( y_out = 0; y_out < n; ++y_out )
	{
		y_in = (y_out - yc_out) * scale + yc_in;

		for ( x_out = 0; x_out < n; ++x_out )
		{ 
			x_in = (x_out - xc_out) * scale + xc_in;

			x1 = (int)x_in;
			x2 = (int)(x_in + 1);
			y1 = (int)y_in;
			y2 = (int)(y_in + 1);

			if ( x1 < 0 || x2 >= nx || y1 < 0 || y2 >= ny )
				continue;

			t = x_in - x1;
			u = y_in - y1;

			map[y_out][x_out] = (1-t)*(1-u)*primary[y1][x1] + 
					   	t*(1-u)*primary[y1][x2] +
                			      	t*u*primary[y2][x2] + 
						(1-t)*u*primary[y2][x1];
		}
	}

	Free_image( primary );

} /* Primary_map */

/****************************************************************************
*  Secondary_map :
*	Read in map of surface errors on the secondary mirror; scale it to
*	match the pupil diameter; shift it to match the off-axis projection
*	of the primary onto the secondary
*
*  Inputs :
*	map : 2D array into which map is stored; should be previously
*		allocated
*	paraminfo : pointer to parameter information structure
****************************************************************************/ 
void Secondary_map( float **map, paramstruct *paraminfo )
{
	char	filename[MAX_STRING];
	float	**secondary, delta, r_in, r_out, x_in, y_in;
	float	scale, t, u;
	float	pri_sec, sec_image, exit_dist, exit_diam, pupil_diam;
	float	y_shift, z_shift;
	int	n, nx, ny, x1, x2, y1, y2, x_out, y_out, xc_in, yc_in;
	int	xc_out, yc_out;

	n = paraminfo->nyquist_grid_size;

	/* read in secondary mirror map; map is assumed to be       *
	 * oriented so that the +X array axis is the +Z telescope   *
	 * axis and the +Y array axis is the +Y telescope axis;     *
	 * the map is in microns of wavefront error (twice the      *
	 * surface error). 					    */

	Default_dir( filename );

	if ( paraminfo->use_new_secondary == 1 )
		strcat( filename, "new_sirtf_secondary.fits" );
	else
		strcat( filename, "sirtf_secondary.fits" );

	secondary = Read_FITS( filename, &nx, &ny ); 

	delta = 0.5208;		/* spacing in mm of map in FITS file */

	pri_sec = 887.5454;	/* primary-to-secondary separation (mm) */
	sec_image = 1324.54;	/* secondary-to-image distance (mm) */
	exit_dist = 1452.0;	/* exit pupil distance (mm) */
	exit_diam = 121.0;	/* exit pupil diameter (mm) */

	/* pupil_diam = pupil diameter on secondary mirror (mm) */

	pupil_diam = exit_diam * sec_image / exit_dist;

	/* z_shift, y_shift = location of chief ray on  *
	 * secondary in mm from the secondary center.   */

	z_shift = pri_sec * tan(paraminfo->x_offset/60.0 * M_PI/180.0);
	y_shift = pri_sec * tan(paraminfo->y_offset/60.0 * M_PI/180.0);
  
	r_in = 0.5 * pupil_diam / delta;  /* pupil radius in map pixels */ 
	xc_in = nx / 2;     /* pupil center in map pixels */
	yc_in = ny / 2;

	xc_out = n / 2;
	yc_out = n / 2;
	r_out = n / 4.0;	/* pupil radius in aberration array pixels */
	scale = r_in / r_out;

	/* scale and shift map (with interpolation) */

	for ( y_out = 0; y_out < n; ++y_out )
	{
		y_in = (y_out - yc_out) * scale + yc_in + y_shift;
		for ( x_out = 0; x_out < n; ++x_out )
		{ 
			x_in = (x_out - xc_out) * scale + xc_in + z_shift;

			x1 = (int)x_in;
			x2 = (int)(x_in + 1);
			y1 = (int)y_in;
			y2 = (int)(y_in + 1);

			if ( x1 < 0 || x2 >= nx || y1 < 0 || y2 >= ny )
				continue;

			t = x_in - x1;
			u = y_in - y1;

			map[y_out][x_out] += (1-t)*(1-u)*secondary[y1][x1] + 
						   t*(1-u)*secondary[y1][x2] +
                			   	   t*u*secondary[y2][x2] + 
						   (1-t)*u*secondary[y2][x1];
		}
	}

	Free_image( secondary );

} /* Secondary_map */

/****************************************************************************
*  Mirror_map :
*	Read in maps of surface errors on the mirrors; rotate them to match
*	the detector coordinate system.
*
*  Inputs :
*	opd : 2D array containing aberration function; should already
*	        contain low-frequency aberrations
*	paraminfo : pointer to parameter information structure
****************************************************************************/ 
void Mirror_map( float **opd, paramstruct *paraminfo )
{
	float 	**map, cos_t, sin_t, x, y, x_in, y_in, t, u, theta;
	int	n, x1, x2, y1, y2, x_out, y_out, xc, yc;

	n = paraminfo->nyquist_grid_size;
	map = Alloc_image( n, n );
 
	Primary_map( map, paraminfo );
	Secondary_map( map, paraminfo );

	/* rotate map (with interpolation) to match detector coordinate system */

	theta = paraminfo->detector_rotation * M_PI / 180.0;
	cos_t = cos(theta);
	sin_t = sin(theta);
	xc = n / 2;
	yc = n / 2;

	for ( y_out = 0; y_out < n; ++y_out )
	{
		y = y_out - yc;
		for ( x_out = 0; x_out < n; ++x_out )
		{ 
			x = x_out - xc;
			x_in = x * cos_t - y * sin_t + xc;
			y_in = x * sin_t + y * cos_t + yc;

			x1 = (int)x_in;
			x2 = (int)(x_in + 1);
			y1 = (int)y_in;
			y2 = (int)(y_in + 1);

			if ( x1 < 0 || x2 >= n || y1 < 0 || y2 >= n )
				continue;

			t = x_in - x1;
			u = y_in - y1;

			opd[y_out][x_out] += (1-t)*(1-u)*map[y1][x1] + t*(1-u)*map[y1][x2] +
                			   t*u*map[y2][x2] + (1-t)*u*map[y2][x1];
		}
	}

	Free_image( map );

} /* Mirror_map */

