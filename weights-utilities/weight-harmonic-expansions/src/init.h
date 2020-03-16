#ifndef _H_RAY_GRID_INIT
#define _H_RAY_GRID_INIT

/* --------------------------- attenuation map utilites -------------------- */

int att_values_alloc(double**** att_values, int npixels); // allocate memory for attenuation map
int att_values_init(double*** att_values, char* att_filename, int npixels); // read attenuation map from file
int att_values_clean(double*** att_values, int npixels); // freemem
/* ------------------------------------------------------------------------- */

/* ----------------------------- beam grid utilites ----------------------- */

int ray_grid_alloc(double** phi, int nphi); // allocate memory for divergent beams (only outgoing directions from a fixed point)
int ray_grid_init(double* phi, int nphi);   // init outgoing directions (unfiform anges in [0, 2pi])
int ray_grid_clean(double* phi);            // clean allocated memory
/* --------------------------------------------------------------------------  */

/*--------------------------- threads chunks (filenames and files) utilites ------------------------ */

int chunks_filenames_alloc(char***** chunks_filenames, int nchunks, int ndegrees, int max_len);     // allocate memory for filenames (nthreads x ndegree x 2 (real/imag))
int chunks_filenames_init(char**** chunks_filenames, char*** output_filenames, int nchunks, int ndegrees); // set filenames for temporary files used by threads
int chunks_files_clean(char**** chunks_filenames, int chunks, int ndegrees);	      // clean temporary files - call 'rm -f ...'
int chunks_filenames_clean(char**** chunks_filenames, int nchunks, int ndegrees);     // clean allocated memory for filenames
/* --------------------------------------------------------------------------- */

int output_filenames_alloc(char**** output_filenames, int ndegrees, int max_len);
int output_filenames_init(char*** output_filenames, char* output_filename, int ndegrees);
int output_filenames_clean(char*** output_filenames, int ndegrees);

#endif
