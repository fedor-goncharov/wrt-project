#ifndef _H_RADON_INIT
#define _H_RADON_INIT

/*depends on GSL library*/

/* Grid parameters utilites */

int parameters_alloc(int** parameters);
int parameters_init(char* filename, int* parameters);
int parameters_clean(int* parameters);

/* Radon grid utilites */

int radon_grid_alloc(int nphi, int ntheta, int nshift, double** phi, double** theta, double** shift);
int radon_grid_init(int nphi, int ntheta, int nshift, double* phi, double* theta, double* shift);
int radon_grid_clean(double* phi, double* theta, double* shift);


/* Threads filenames utilites */

int chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int chunks_filenames_init(char* filename, int nchunks, char** nchunks_filenames);
int chunks_files_clean(int nthreads, char** chunks_filenames);
int chunks_filenames_clean(int nchunks, char** nchunks_filenames);

#endif
