#ifndef _H_RAY_INIT
#define _H_RAY_INIT

/* parameters utilites */

int parameters_alloc(int** parameters);
int parameters_init(char* filename, int* parameters);
int parameters_clean(int* parameters);

/* ray utilites */

int ray_grid_alloc(int nshift, int nphi, double** shift, double** phi);
int ray_grid_init(int nshift, int nphi, double* shift, double* phi);
int ray_grid_clean(double* shift, double* phi);

/* Threads filenames utilites */

int chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int chunks_filenames_init(char* filename, int nchunks, char** nchunks_filenames);
int chunks_filenames_clean(int nchunks, char** nchunks_filenames);
int chunks_files_clean(int nthreads, char** chunks_filenames);

#endif
