#ifndef __RAY_GRID_INIT_H
#define __RAY_GRID_INIT_H

/* Grid parameters utilites */

int parameters_alloc(int** parameters);
int parameters_init(char* filename, int* parameters);
int parameters_clean(int* parameters);

/* Test-function utilities */

int tfunction_values_alloc(int ngrid, double**** values);
int tfunction_values_init(char* input_filename, int ngrid, double*** values);
int tfunction_values_clean(int ngrid, double*** values);

/* Radon grid utilites */

int ray_grid_alloc(int nshift, int nphi, double** shift, double** phi);
int ray_grid_init(int nshift, int nphi, double* shift, double* phi);
int ray_grid_clean(double* shift, double* phi);


/* Threads filenames utilites */

int chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int chunks_filenames_init(char* filename, int nchunks, char** nchunks_filenames);
int chunks_files_clean(int nthreads, char** chunks_filenames);
int chunks_filenames_clean(int nchunks, char** nchunks_filenames);
void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks);

#endif
