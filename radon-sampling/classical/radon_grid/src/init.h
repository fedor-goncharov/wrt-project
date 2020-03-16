#ifndef __H_RADON_INIT
#define __H_RADON_INIT

/* Initial values utilites */

int function_values_alloc(int npixels, double**** function_values);
int function_values_init(char* filename, int npixels, double*** function_values);
int function_values_clean(int npixels, double*** function_values);

/* Radon grid utilites */

int radon_grid_alloc(int nphi, int ntheta, int nshift, double** phi, double** theta, double** shift);
int radon_grid_init(int nphi, int ntheta, int nshift, double* phi, double* theta, double* shift);
int radon_grid_clean(double* phi, double* theta, double* shift);


/* Threads filenames utilites */

int  chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int name_max_size);
int  chunks_filenames_init(char* filename, int nchunks, char** nchunks_filenames);
void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks);
int  chunks_files_clean(int nthreads, char** chunks_filenames);
int  chunks_filenames_clean(int nchunks, char** nchunks_filenames);

#endif
