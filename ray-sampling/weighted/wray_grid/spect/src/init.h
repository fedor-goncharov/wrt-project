#ifndef __H_RAY_GRID_INIT
#define __H_RAY_GRID_INIT

/* ----------------------------parameters utilites --------------------------*/

int parameters_alloc(int** parameters);
int parameters_init(char* filename, int* parameters);
int parameters_clean(int* parameters);
/* ------------------------------------------------------------------------- */

/* --------------------------- test-function utilities --------------------- */

int func_values_alloc(int ngrid_func, double**** values);
int func_values_init(char* filename, int ngrid_func, double*** values);
int func_values_clean(int ngrid_func, double*** values);

/* ------------------------------------------------------------------------- */


/* --------------------------- attenuation map utilites -------------------- */


int att_values_alloc(int ngrid_func, double**** values);
int att_values_init(char* filename, int ngrid_func, double*** values);
int att_values_clean(int ngrid_func, double*** values);


/* ------------------------------------------------------------------------- */


/* ----------------------------- ray grid utilites ----------------------- */

int ray_grid_alloc(int nshift, int nphi, double** shift, double** phi);
int ray_grid_init(int nshift, int nphi, double* shift, double* phi);
int ray_grid_clean(double* shift, double* phi);
/* --------------------------------------------------------------------------  */

/*--------------------------- threads chunks utilites ------------------------ */

int chunks_filenames_alloc(char* filename, int nchunks, char*** nchunks_filenames, int size);
int chunks_filenames_init(char* filename, int nchunks, char** nchunks_filenames);
int chunks_files_clean(int nthreads, char** chunks_filenames);
int chunks_filenames_clean(int nchunks, char** nchunks_filenames);

/* --------------------------------------------------------------------------- */
#endif
