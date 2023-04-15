#define VERSION "2.0"

#define M_PI            3.14159265358979323846

#define MAX_POSITIONS   100
#define MAX_WAVELENGTHS 150
#define MAX_STRING      255
#define LAST_ZERNIKE    22

#define IRAC_35         1
#define IRAC_45         2
#define IRAC_63         3
#define IRAC_80         4
#define MIPS_24         5
#define MIPS_70_WF      6
#define MIPS_70_SUPER   7
#define MIPS_70_SED     8     
#define MIPS_160        9     
#define IRS_LONG_LOW    10
#define IRS_LONG_HIGH   11
#define IRS_SHORT_LOW   12
#define IRS_SHORT_HIGH  13
#define IRS_PEAKUP_BLUE 14
#define IRS_PEAKUP_RED  15

#define PSF_FILE	1
#define PUPIL_FILE	2
#define OPD_FILE	3
#define NYQUISTPSF_FILE 4

#define SEPARATOR "-----------------------------------------------------------\n"

int	NumThreads;
#define MAX_THREADS 5

typedef struct paramstruct {
	int	adjust_field_aberrations;       /* if 1, adjust for field aberrations */
	char	bandpass_file[MAX_STRING];	/* name of file containing bandpass curve */ 
        int     camera;	
	char	camera_name[MAX_STRING];
	int	current_pos;			/* current field position index */
	float	focus, xcoma, ycoma, xastig, yastig, spherical;
	int	is_subsampled;
	float	jitter;				/* RMS jitter */
	float   kernel[3][3];			/* IRAC charge diffusion kernel */
        int     num_positions;
        int     num_waves;			/* number of wavelengths in polychromatic psf */
        int     nyquist_grid_size;
        int     psf_grid_size;			/* array size of pixel-integrated PSF */
	int	smart_skip;			/* if 1, skip low weights */
        char    rootname[MAX_STRING];
	char    pupilfile[MAX_STRING];
	float	scan_offset;			/* MIPS scan position in arcsec */
        char    spectrum_file[MAX_STRING];
	int	use_map;			/* if 1, use primary mirror zonal errors */
	int	use_new_secondary;		/* if 1, use new secondary mirror map */
	int	use_kernel;
        float   wavelength[MAX_WAVELENGTHS];
	float	wavelength_mag;
	float	weight[MAX_WAVELENGTHS];
	float	weight_limit;			/* criterion for smart_skip */
	int	write_pupil, write_wave, write_nyquist_psf;
        float   x[MAX_POSITIONS];
	float	y[MAX_POSITIONS];
	float	zval[50];

	float	caber[6][4][4];			/* field dependent aberration coefficients */
	float	mipscaber[3][6][4][4];
	float	mipsaber[3][6];
	float	mipsscan[3];
	float	detector_rotation;		/* rotation of detector from CTA Z-Y axis */
	int  	detector_x_flip;		/* non-zero if flip -x to +x */
	int  	detector_y_flip;		/* non-zero if flip -y to +y */
	float	diam_telescope_mm;
	float	f_ratio;
	float	min_wavelength, max_wavelength;
	float	min_xfield, max_xfield;
	float	min_yfield, max_yfield;
	float	ota_rotation;			/* CTA spider rotation relative to CTA Z-Y axis */
        float   det_pixel_size_arcsec_x;
        float   det_pixel_size_arcsec_y;
	float	psf_pixel_size_arcsec_x;
	float	psf_pixel_size_arcsec_y;
	float	secondary_radius;		/* central obscuration radius / pupil radius */
	float	xdelta, ydelta;			/* offset from camera center in CTA Z-Y system */
	float	xfield, yfield;			/* field position after scan offset */
	float	x_offset, y_offset;		/* CTA field position of detector center (arcmin) */
	int	x_pixels, y_pixels;
	float	scan_angle;
} paramstruct;

typedef struct complex { float r, i; } complex;

#define FLOATING  1
#define COMPLEX   2

/* entry.c */
void Init_entry_list( void );
void Write_entry_list( FILE *file );
void Delete_entry_list( void );
const char *Get_entry( FILE *file );

/* fft.c */
void fft2d( complex *image, int dim, int direction );

/* fitsio.c */
float **Read_FITS( char *filename, int *nx, int *ny );
void Write_FITS( char *filename, float *image, int nx, int ny, 
	 	 paramstruct *paraminfo, int type );

/* gridsize.c */
float Max_psf_diameter( paramstruct *paraminfo );
float Get_pixel_scale( paramstruct *paraminfo );
void Set_grid_size( float psf_diam_arcsec, paramstruct *paraminfo );

/* image.c */
float **Alloc_image( int nx, int ny );
void Free_image( float **image );
complex **Alloc_complex_image( int nx, int ny );
void Free_complex_image( complex **image );
void Shift_to_center( float **image, int dim );
void Shift_complex_to_center( complex **image, int dim );
void Shift_complex_to_origin( complex **image, int dim );

/* intpsf.c */
void Integrate_psf( complex **psf_in, float scale_in, int size_in,
                   float **psf_out, float scale_out_x, float scale_out_y, 
		   int size_out, float jitter );

/* map.c */
void Mirror_map( float **opd, paramstruct *paraminfo );

/* misc.c */
char *Find_string( char *string1, char *string2 );
void Delete_file( char *filename );
void Default_dir( char *filename );
void Write_image( char *filename, void *image, int nx, int ny, int data_type );
void Read_image( char *filename, void *image, int nx, int ny, int data_type );

/* monopsf.c */
void Compute_mono_psf( paramstruct *paraminfo, int wave_index, 
	float **monopsf, complex **pupil );

/* opd.c */
void Compute_opd( paramstruct *paraminfo, float **opd, int n );

/* param.c */
void Write_parameters( char *filename, paramstruct *paraminfo );
void Read_parameters( char *filename, paramstruct *paraminfo );

/* polypsf.c */
void Compute_poly_psf( paramstruct *paraminfo,  complex **pupil, float **polypsf, float **monopsf );

/* psf.c */
void Compute_nyquist_psf( complex **pupil, int dim );

/* pupil.c */
void Compute_pupil( paramstruct *paraminfo, int n, float **image );

/* rdpupil.c */
void Read_general( FILE *file, paramstruct *paraminfo );
void Read_ota( FILE *file, paramstruct *paraminfo );
void Read_instrument( FILE *file, paramstruct *paraminfo );
void Read_aberrations( FILE *file, paramstruct *paraminfo );
void Read_pupil_file( paramstruct *paraminfo );

/* spectrum.c */
void Blackbody( paramstruct *paraminfo );
void Power_law_nu( paramstruct *paraminfo );
void Power_law_lambda( paramstruct *paraminfo );
void Spectrum_file( paramstruct *paraminfo );

