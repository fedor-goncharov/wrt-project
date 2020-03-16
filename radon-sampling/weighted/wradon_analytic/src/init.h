#ifndef _H_RADON_INIT
#define _H_RADON_INIT

/*depends on GSL library*/

/* Grid parameters utilites */

int init_parameters_alloc(int** parameters);
int init_parameters(char* filename, int* parameters);
int clean_parameters(int* parameters);

/* Radon grid utilites */

int init_radon_grid_alloc(int nphi, int ntheta, int nshift, double** phi, double** theta, double** shift);
  /*
   * phi   -- uniform grid on S^1 = [0,2pi]
   * theta -- Gauss's angles on [0, pi]
   * shift -- uniform points on [-1, 1]
   */
int init_radon_grid(int nphi, int ntheta, int nshift, double* phi, double* theta, double* shift);
int clean_radon_grid(double* phi, double* theta, double* shift);


/* Threads filenames utilites */

int init_chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int init_chunks_filenames(char* filename, int nchunks, char** nchunks_filenames);
int clean_chunks_files(int nthreads, char** chunks_filenames);
int clean_chunks_filenames(int nchunks, char** nchunks_filenames);

#endif