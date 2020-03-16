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
  char buffer[ MAX_COMMAND_SIZE ];
  
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
  return 0;	//exit sucess
  
}

int parameters_clean(int* parameters) {
  free(parameters);
  
  return 0;
}

/* Radon grid */

int radon_grid_alloc(int nphi, int ntheta, int nshift, double** phi, double** theta, double** shift) {
    *phi = (double*)malloc(nphi * sizeof(double));
    *theta = (double*)malloc(ntheta * sizeof(double));
    *shift = (double*)malloc(nshift * sizeof(double));
    
    return 0;
}
int radon_grid_init(int nphi, int ntheta, int nshift, double* phi, double* theta, double* shift) {
  
  //init phi 
  const double dphi = (2 * M_PI) / nphi;
  int iphi;
  for (iphi = 0; iphi < nphi; ++iphi) {
      phi[iphi] = dphi * iphi;
  }
  
  //init shift
  const double dshift = 2.0 / (nshift - 1);
  int ishift; 
  for (ishift = 0; ishift < nshift; ++ishift) {
      shift[ishift] = (-1.0) + ishift * dshift;
  }
  
  //init theta (using GSL library to compute Gaussian points)
  int itheta;
  gsl_integration_glfixed_table *gsl_quad_data = gsl_integration_glfixed_table_alloc(ntheta);
  for (itheta = 0; itheta < ntheta; ++itheta) {
     double gauss_point, gauss_weight;
     gsl_integration_glfixed_point(-1.0, 1.0, itheta, &gauss_point, &gauss_weight, gsl_quad_data);
     theta[itheta] = acos(gauss_point);
  }
  
  gsl_integration_glfixed_table_free(gsl_quad_data);
  
  return 0;
}
int radon_grid_clean(double* phi, double* theta, double* shift) {
   if (phi != NULL)
     free(phi);
   if (theta != NULL)
     free(theta);
   if (shift != NULL)
     free(shift);
   
   return 0;
}

/* Threads filenames */

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
    char buffer[ MAX_COMMAND_SIZE ];
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
