#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#define MAX_COMMAND_SIZE 128

/* Grid parameters */ 

int parameters_alloc(int** parameters) {
  *parameters = (int*)malloc(4 * sizeof(int));
  
  return 0;
}

int parameters_init(char* filename, int* parameters) {
  
  FILE *f;
  char buffer[MAX_COMMAND_SIZE];
  
  f = fopen(filename, "r");
  if (f == NULL) {
    fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
    return errno;
  }
  
  int count;
  count = fscanf(f, "%d %[^\n]\n", &(parameters[0]), buffer);
  if (count == EOF) {
     if (ferror(f)) {
        perror("Init parameters : fscanf parameters 0 : EOF failure\n");
     }  else {
        fprintf(stderr, "Error: fscanf parameters 0 : reached end of file before it was expected. \n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 0 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  count = fscanf(f, "%d %[^\n]\n", &(parameters[1]), buffer);
  if (count == EOF) {
     if (ferror(f)) {
        perror("Init parameters : fscanf parameters 1 : EOF failure\n");
     }  else {
        fprintf(stderr, "Error: fscanf parameters 1 : reached end of file before it was expected. \n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 1 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  count = fscanf(f, "%d %[^\n]\n", &(parameters[2]), buffer);
  if (count == EOF) {
     if (ferror(f)) {
        perror("Init parameters : fscanf parameters 2 : EOF failure\n");
     }  else {
        fprintf(stderr, "Error: fscanf parameters 2 : reached end of file before it was expected. \n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 2 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  count = fscanf(f, "%d %[^\n]\n", &(parameters[3]), buffer);
  if (count == EOF) {
     if (ferror(f)) {
        perror("Init parameters : fscanf parameters 3 : EOF failure\n");
     }  else {
        fprintf(stderr, "Error: fscanf parameters 3 : reached end of file before it was expected. \n");
     }
     return -1;
  } else if (count != 2) {
         fprintf(stderr, "Error: fscanf parameters 3 : matching failure, fscanf expected 1 integer value.\n");
    return -1;
  }
  
  fclose(f);
  return 0;
  
}

int parameters_clean(int* parameters) {
  free(parameters);
  return 0;
}

/* Initial values of the test-function */

int radon_all_alloc(int nphi, int ntheta, int nshift, 
  double** phi, double** theta, double** theta_weights, double** shift, 
  double**** radon_values) {

  // allocate memory for grid in Radon space
  *phi = (double*)malloc(sizeof(double) * nphi);
  *theta = (double*)malloc(sizeof(double) * ntheta);
  *theta_weights = (double*)malloc(sizeof(double) * ntheta);
  *shift = (double*)malloc(sizeof(double) * nshift);

  int i,j;
  
  // allocate memory for values of Radon transforms
  *radon_values = (double***)malloc(sizeof(double**) * nphi);
  for (i = 0; i < nphi; ++i) {
      (*radon_values)[i] = (double**)malloc(sizeof(double*) * ntheta);
  }
  
  for (i = 0; i < nphi; ++i) {
      for (j = 0; j < ntheta; ++j) {
          (*radon_values)[i][j] = (double*)malloc(sizeof(double) * nshift);
      }
  }
  return 0;  
}

int radon_values_init(char* filename, double*** radon_values, int nphi, int ntheta, int nshift) {
  
  FILE *f;
  f = fopen(filename, "r");
  
  if (f == NULL) {
    fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
    return errno;
  }
  
  int i_phi, i_theta, i_shift;
  int count;
  
  for (i_shift = 0; i_shift < nshift; ++i_shift)
    for (i_phi = 0; i_phi < nphi; ++i_phi)
      for (i_theta = 0; i_theta < ntheta; ++i_theta) {
	      
        //read file line-by-line 
        double tmp1, tmp2, tmp3;
	      count = fscanf(f, "%lf, %lf, %lf, %lf\n", &tmp1, &tmp2, &tmp3, &(radon_values[i_phi][i_theta][i_shift]));
	      
        if (count == EOF) {
	        if (ferror(f)) {
	          perror("Init test-function values : fscanf read value : EOF failure\n");
	        } else {
	          fprintf(stderr, "Error : init test-function : fscanf : reached the end of file before it was expected.\n");
	        }
	        return -1;
	      } else if (count != 4) {
	        fprintf(stderr, "Error: init test-function : fscanf : matching failure, fscanf expected 4 floating-point value.\n");
	        return -1;
	      }
      }

  fclose(f);  
  return 0;
}
int radon_values_clean(double*** radon_values, int nphi, int ntheta, int nshift) {
  
  int i,j;
  for (i = 0; i < nphi; ++i) {
      for (j = 0; j < ntheta; ++j) {
          free(radon_values[i][j]);
      }
  }
  for (i = 0; i < nphi; ++i) {
      free(radon_values[i]);
  }
  free(radon_values);
  return 0;
}


int radon_grid_init(int nphi, int ntheta, int nshift, double* phi, double* theta, double* theta_weights, double* shift) {
  
  // init phi 
  const double dphi = (2 * M_PI) / nphi;
  int i_phi;
  for (i_phi = 0; i_phi < nphi; ++i_phi)
    phi[i_phi] = dphi * i_phi;
  
  // init shift
  const double dshift = 2.0 / (nshift - 1);
  int i_shift; 
  for (i_shift = 0; i_shift < nshift; ++i_shift)
    shift[i_shift] = (-1.0) + i_shift * dshift;
  
  // init theta and theta_weight (using GSL to compute Gauss-Legendre quadrature rule)
  int i_theta;
  gsl_integration_glfixed_table *gsl_quad_data = gsl_integration_glfixed_table_alloc(ntheta);
  
  for (i_theta = 0; i_theta < ntheta; ++i_theta) {
     double gauss_point, gauss_weight;
     gsl_integration_glfixed_point(-1.0, 1.0, i_theta, &gauss_point, &gauss_weight, gsl_quad_data);
     theta[i_theta] = acos(gauss_point);
     theta_weights[i_theta] = gauss_weight;
  }
  gsl_integration_glfixed_table_free(gsl_quad_data);
  
  return 0;
}
int radon_grid_clean(double* phi, double* theta, double* theta_weights, double* shift) {
   if (phi != NULL)
     free(phi);
   if (theta != NULL)
     free(theta);
   if (theta_weights != NULL)
     free(theta_weights);
   if (shift != NULL)
     free(shift);
   
   return 0;
}

/* threads filenames */

int chunks_filenames_alloc(char* filename, int nchunks, char*** chunks_filenames, int size) {
  (*chunks_filenames) = (char**)malloc(sizeof(char*) * nchunks);
  int i;
  for (i = 0; i < nchunks; ++i) {
      (*chunks_filenames)[i] = (char*)malloc(sizeof(char) * size);
  }
  
  return 0;
}
int chunks_filenames_init(char* filename, int nchunks, char** chunks_filenames) {
  int i;
  for (i = 0; i < nchunks; ++i) {
    sprintf(chunks_filenames[i], "%d_%s", i, filename);
  }
  return 0;
}

int chunks_files_clean(int nchunks, char** chunks_filenames) {
  int i_file;
  for (i_file = 0; i_file < nchunks; ++i_file) {
    char buffer[MAX_COMMAND_SIZE];
    sprintf(buffer, "rm -f %s", chunks_filenames[i_file]);
    int sys_call = system(buffer);
    if (sys_call == -1) {
      fprintf(stderr, "Warning : clean chunks files : system call (rm -f ...) : failed to clean file %s.\n", chunks_filenames[i_file]); 
    }
  }
  return 0;
}

int chunks_filenames_clean(int nchunks, char** chunks_filenames) {
  int i;
  for (i = 0; i < nchunks; ++i) {
      free(chunks_filenames[i]); 
  }
  free(chunks_filenames);
  return 0;
}
