#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define MAX_COMMAND_SIZE 64
#define MAX_COMMENT_SIZE 128

/* ----------------------------parameters utilites --------------------------*/

int parameters_alloc(int** parameters) {
  *parameters = (int*)malloc(4 * sizeof(int));
  
  return 0;
}
int parameters_init(char* filename, int* parameters) {
  
  FILE *f;
  char buffer[MAX_COMMENT_SIZE];
  
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
/* ------------------------------------------------------------------------- */

/* -------------------- initial values for the test-function --------------- */

int func_values_alloc(int ngrid, double**** values) {
  int i,j;
  
  *values = (double***)malloc(sizeof(double**) * ngrid);
  for (i = 0; i < ngrid; ++i) {
      (*values)[i] = (double**)malloc(sizeof(double*) * ngrid);
  }
  
  for (i = 0; i < ngrid; ++i) {
      for (j = 0; j < ngrid; ++j) { 
          (*values)[i][j] = (double*)malloc(sizeof(double) * ngrid);
      }
  }
  return 0;  
}

int func_values_init(char* filename, int ngrid, double*** values) {
  FILE *f;
  f = fopen(filename, "r");
  
  if (f == NULL) {
     fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
     return errno;
  }
  
  int xi, yi, zi;
  int count;
  for (xi = 0; xi < ngrid; ++xi)
    for (yi = 0; yi < ngrid; ++yi)
      for (zi = 0; zi < ngrid; ++zi) {
	
	//read file line-by-line 
	count = fscanf(f, "%lf\n", &(values[xi][yi][zi]));
	if (count == EOF) {
	   if (ferror(f)) {
	      perror("Init test-function values : fscanf read value : EOF failure\n");
	   } else {
	      fprintf(stderr, "Error : init test-function : fscanf : reached the end of file before it was expected.\n");
	   }
	   return -1;
	} else if (count != 1) {
	      fprintf(stderr, "Error: init test-function : fscanf : matching failure, fscanf expected 1 floating-point value.\n");
	      return -1;
	  }
  }
      
  fclose(f);  
  return 0;	//exist sucess
}
int func_values_clean(int ngrid, double*** values) {
  int i,j;
  for (i = 0; i < ngrid; ++i) {
      for (j = 0; j < ngrid; ++j) {
          free(values[i][j]);
      }
  }
  for (i = 0; i < ngrid; ++i) {
      free(values[i]);
  }
  free(values);
  return 0;
}
/* ------------------------------------------------------------------------- */

/* -------------------- initial values for the attenuation map ------------- */
int att_values_alloc(int ngrid, double**** att_values) {
  int i,j;
  
  *att_values = (double***)malloc(sizeof(double**) * ngrid);
  for (i = 0; i < ngrid; ++i) {
      (*att_values)[i] = (double**)malloc(sizeof(double*) * ngrid);
  }
  
  for (i = 0; i < ngrid; ++i) {
      for (j = 0; j < ngrid; ++j) { 
          (*att_values)[i][j] = (double*)malloc(sizeof(double) * ngrid);
      }
  }
  return 0;
  
}
int att_values_init(char* filename, int ngrid, double*** att_values) {
  FILE *f;
  f = fopen(filename, "r");
  
  if (f == NULL) {
     fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
     return errno;
  }
  
  int xi, yi, zi;
  int count;
  for (xi = 0; xi < ngrid; ++xi)
    for (yi = 0; yi < ngrid; ++yi)
      for (zi = 0; zi < ngrid; ++zi) {
	
	//read file line-by-line 
	count = fscanf(f, "%lf\n", &(att_values[xi][yi][zi]));
	if (count == EOF) {
	   if (ferror(f)) {
	      perror("Init test-function values : fscanf read value : EOF failure\n");
	   } else {
	      fprintf(stderr, "Error : init test-function : fscanf : reached the end of file before it was expected.\n");
	   }
	   return -1;
	} else if (count != 1) {
	      fprintf(stderr, "Error: init test-function : fscanf : matching failure, fscanf expected 1 floating-point value.\n");
	      return -1;
	}
     }
      
  fclose(f);  
  return 0;	//exist sucess
}
int att_values_clean(int ngrid, double*** att_values) {
  int i,j;
  for (i = 0; i < ngrid; ++i) {
      for (j = 0; j < ngrid; ++j) {
          free(att_values[i][j]);
      }
  }
  for (i = 0; i < ngrid; ++i) {
      free(att_values[i]);
  }
  free(att_values);
  
  return 0;
}
/* ------------------------------------------------------------------------- */

/* --------------------------- ray grid ------------------------------------ */

int ray_grid_alloc(int nshift, int nphi, double** shift, double** phi) {
  *shift = (double*)malloc(nshift * sizeof(double));  
  *phi = (double*)malloc(nphi * sizeof(double));    
  return 0;
}
int ray_grid_init(int nshift, int nphi, double* shift, double* phi) {
  
  /*
   * phi --   uniform grid on S^1 = [0,2pi]
   * theta -- Gauss's angles on [0, pi]
   * shift -- uniform points on [-1, 1]
   */
  
  //init shift
  const double dshift = 2.0 / (nshift - 1);
  int i_shift; 
  for (i_shift = 0; i_shift < nshift; ++i_shift)
    shift[i_shift] = (-1.0) + i_shift * dshift;
  
  //init phi 
  const double dphi = (2 * M_PI) / nphi;
  int i_phi;
  for (i_phi = 0; i_phi < nphi; ++i_phi)
    phi[i_phi] = dphi * i_phi;
  
  return 0;
}
int ray_grid_clean(double* shift, double* phi) {
   
   if (shift != NULL)
     free(shift);
   
   if (phi != NULL)
     free(phi);
   
   return 0;
}
/* ------------------------------------------------------------------------- */

/* --------------------------- threads filenames --------------------------- */

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
/* ------------------------------------------------------------------------- */
