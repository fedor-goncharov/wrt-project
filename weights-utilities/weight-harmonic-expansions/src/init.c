#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif 

#define MAX_COMMAND_SIZE 128

/* -------------------- initial values for the attenuation map ------------- */
int att_values_alloc(double**** att_values, int ngrid) {
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

int att_values_init(double*** att_values, char* att_filename, int ngrid) {
  FILE *f;
  f = fopen(att_filename, "r");
  
  if (f == NULL) {
     fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), att_filename);
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

int att_values_clean(double*** att_values, int ngrid) {
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

/* --------------------------- beam grid ------------------------------------ */

int ray_grid_alloc(double** phi, int nphi) {
  *phi = (double*)malloc(nphi * sizeof(double));    
  return 0;
}

int ray_grid_init(double* phi, int nphi) {

  const double dphi = (2*M_PI) / nphi;
  int i_phi;
  for (i_phi = 0; i_phi < nphi; ++i_phi)
    phi[i_phi] = dphi * i_phi;

  return 0;
}

int ray_grid_clean(double* phi) {
   if (phi != NULL)
     free(phi);
   return 0;
}
/* ------------------------------------------------------------------------- */


/* --------------------------- threads filenames --------------------------- */

int chunks_filenames_alloc( char***** chunks_filenames, int nchunks, int ndegrees, int size) {
  
  // allocate memory for chunks_filenames
  (*chunks_filenames) = (char****)malloc(sizeof(char***) * nchunks);

  int i, j;
  for (i = 0; i < nchunks; ++i) {
      (*chunks_filenames)[i] = (char***)malloc(sizeof(char**) * ndegrees);
      
      for (j = 0; j < ndegrees; ++j) {
        (*chunks_filenames)[i][j] = (char**)malloc(sizeof(char*) * 2); // allocate real/imag

        (*chunks_filenames)[i][j][0] = (char*)malloc(sizeof(char) * size);
        (*chunks_filenames)[i][j][1] = (char*)malloc(sizeof(char) * size);
      }
  }

  return 0;
}

int chunks_filenames_init(char**** chunks_filenames, char*** output_filenames, int nchunks, int ndegrees) {
  
  int i, j;
  for (i = 0; i < nchunks; ++i) {
    for(j = 0; j < ndegrees; ++j) {
      sprintf(chunks_filenames[i][j][0], "thr%d_%s", i, output_filenames[j][0]);	// real / splitting into threads
      sprintf(chunks_filenames[i][j][1], "thr%d_%s", i, output_filenames[j][1]);  // imag / splitting into threads
    }
  }

  return 0;
}

int chunks_files_clean(char**** chunks_filenames, int nchunks, int ndegrees) {
  int i, j;

  for (i = 0; i < nchunks; ++i) {
    for (j = 0; j < ndegrees; ++j) {
      
      char buffer_real[MAX_COMMAND_SIZE], buffer_imag[MAX_COMMAND_SIZE];
      sprintf(buffer_real, "rm -f %s", chunks_filenames[i][j][0]);
      sprintf(buffer_imag, "rm -f %s", chunks_filenames[i][j][1]);

      int sys_call = system(buffer_real);
      if (sys_call == -1) {
	      fprintf(stderr, "Warning : clean chunks files : system call (rm -f ...) : failed to clean file %s.\n", chunks_filenames[i][j][0]); 
      }
      sys_call = system(buffer_imag);
      if (sys_call == -1) {
        fprintf(stderr, "Warning : clean chunks files : system call (rm -f ...) : failed to clean file %s.\n", chunks_filenames[i][j][1]); 
      }
    }
  }
  return 0;

}

int chunks_filenames_clean(char**** chunks_filenames, int nchunks, int ndegrees) {
  int i, j;
  
  for (i = 0; i < nchunks; ++i) {
    for (j = 0; j < ndegrees; ++j) {
      free(chunks_filenames[i][j][1]);
      free(chunks_filenames[i][j][0]);
    }
  }
  for (i = 0; i < nchunks; ++i) {
    for (j = 0; j < ndegrees; ++j) {
      free(chunks_filenames[i][ndegrees-j-1]);
    }
  }
  for (i = 0; i < nchunks; ++i) {
    free(chunks_filenames[nchunks-i-1]);
  }

  free(chunks_filenames);
  return 0;
}
/* ------------------------------------------------------------------------- */

/* ------------------- output filenames utilities --------------------------- */

int output_filenames_alloc(char**** output_filenames, int ndegrees, int size) {

  (*output_filenames) = (char***)malloc(sizeof(char**) * ndegrees);
  
  int i;
  for (i = 0; i < ndegrees; ++i) {
      (*output_filenames)[i] = (char**)malloc(sizeof(char*) * 2); 
      (*output_filenames)[i][0] = (char*)malloc(sizeof(char) * size);
      (*output_filenames)[i][1] = (char*)malloc(sizeof(char) * size);
  }
  return 0;
}

int output_filenames_init(char*** output_filenames, char* output_filename, int ndegrees) {
  
  int i;
  for (i = 0; i < ndegrees; ++i) {
      sprintf(output_filenames[i][0], "d%d_r_%s", i, output_filename);
      sprintf(output_filenames[i][1], "d%d_i_%s", i, output_filename);
  }
  return 0;
  
}
int output_filenames_clean(char*** output_filenames, int ndegrees) {
  
  int i;
  for (i = 0; i < ndegrees; ++i) {
      free(output_filenames[i][0]);
      free(output_filenames[i][1]);
      free(output_filenames[i]);
  }
  free(output_filenames);
  return 0;
}

/* ------------------------------------------------------------------------- */
