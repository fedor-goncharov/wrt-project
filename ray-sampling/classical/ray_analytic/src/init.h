#ifndef __H_RAY_ANALYTIC_INIT
#define __H_RAY_ANALYTIC_INIT

/* parameters utilites */

int parameters_alloc(int** parameters); // allocate memory for input parameters (grid size, filenames...)
int parameters_init(char* filename, int* parameters); // parse input parameters, set values for variables
int parameters_clean(int* parameters); // deallocate memory for parameters 

/* ray utilites */

int ray_grid_alloc(int nshift, int nphi, double** shift, double** phi); // allocate memory for the grid of rays (z-slices, rays in each slice)
int ray_grid_init(int nshift, int nphi, double* shift, double* phi); // set values for grid coordinates
int ray_grid_clean(double* shift, double* phi); // deallocate memory for grid


/* files utilites */

int chunks_filenames_alloc(char* filename, int nchunks, char*** chunks_filenames, int size); // allocate memory for filenames for each chunk file
int chunks_filenames_init(char* filename, int nchunks, char** chunks_filenames); // initialize filenames 
void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks);
int chunks_files_clean(int nthreads, char** chunks_filenames); // remove  chunks filenames
int chunks_filenames_clean(int nchunks, char** chunks_filenames); // deallocate memory for filenames

#endif
