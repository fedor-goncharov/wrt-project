#ifndef __H_RADJOINT_INIT
#define __H_RADJOINT_INIT

/* Configuration parameters */

int parameters_alloc(int** parameters);
int parameters_init(char* filename, int* parameters);
int parameters_clean(int* parameters);

/* Grid and function in Radon space */

// allocate memory for grid in Radon space
int radon_all_alloc(int nphi, int ntheta, int nshift, double** phi, double** theta, double** theta_weights, double** shift, double**** radon_values);
// initialize grid in Radon space (compute angles, shifts)
int radon_grid_init(int nphi, int ntheta, int nshift, double* phi, double* theta, double* theta_weights, double* shift);
// read values of Radon transforms from file
int radon_values_init(char* filename, double*** radon_values, int nphi, int ntheta, int nshift);
// clean memory from values of Radon transforms
int radon_values_clean(double*** radon_values, int nphi, int ntheta, int nshift);
// clean memory from grid 
int radon_grid_clean(double* phi, double* theta, double* theta_weights, double* shift);


/* Chunks files */

int chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int chunks_filenames_init(char* filename, int nchunks, char** nchunks_filenames);
int chunks_files_clean(int nthreads, char** chunks_filenames);
int chunks_filenames_clean(int nchunks, char** nchunks_filenames);

#endif