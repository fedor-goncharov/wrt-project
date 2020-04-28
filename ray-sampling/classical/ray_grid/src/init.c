#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#define MAX_COMMAND_SIZE 128


/* Initial values of the test-function */

int function_values_alloc(int npixels, int nslices, double**** values) {
  int i,j;
  
  *values = (double***)malloc(sizeof(double**) * nslices);
  for (i = 0; i < nslices; ++i) {
      (*values)[i] = (double**)malloc(sizeof(double*) * npixels);
  }
  
  for (i = 0; i < nslices; ++i) {
      for (j = 0; j < npixels; ++j) { 
          (*values)[i][j] = (double*)malloc(sizeof(double) * npixels);
      }
  }
  return 0;  
}

int function_values_init(char* filename, int npixels, int nslices, double *** values) {
  FILE *f;
  
  f = fopen(filename, "r");
  
  if (f == NULL) {
     fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
     return errno;
  }
  
  int x_i, y_i, z_i;
  int count;
  for (z_i = 0; z_i < nslices; ++z_i)
    for (x_i = 0; x_i < npixels; ++x_i)
      for (y_i = 0; y_i < npixels; ++y_i) {
	
	//read file line-by-line 
	count = fscanf(f, "%lf\n", &(values[z_i][x_i][y_i]));
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
int function_values_clean(int npixels, int nslices, double*** values) {
  int i, j;
  for (i = 0; i < nslices; ++i) {
      for (j = 0; j < npixels; ++j) {
          free(values[i][j]);
      }
  }
  for (i = 0; i < nslices; ++i) {
      free(values[i]);
  }
  free(values);
  
  return 0;
}


/* Radon grid */

int ray_grid_alloc(int nshift, int nphi, double** shift, double** phi) {
  *shift = (double*)malloc(nshift * sizeof(double));  
  *phi = (double*)malloc(nphi * sizeof(double));    
  return 0;
}
int ray_grid_init(int nshift, int nphi, double* shift, double* phi, double radius) {
  //init shift
  const double dshift = 2.0*radius / (nshift - 1);
  int i_shift; 
  for (i_shift = 0; i_shift < nshift; ++i_shift)
    shift[i_shift] = (-1.0)*radius + i_shift * dshift;
  
  //init phi 
  const double dphi = (2 * M_PI) / nphi;
  int i_phi;
  for (i_phi = 0; i_phi < nphi; ++i_phi)
    phi[i_phi] = dphi * i_phi;
  
  return 0; // code sucess
}
int ray_grid_clean(double* shift, double* phi) {
   
   if (shift != NULL)
     free(shift);
   
   if (phi != NULL)
     free(phi);
   
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


void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks) {
  
  FILE *out; //output file
  out = fopen(output_filename, "w");
  
  int i_file;
  char *buffer = (char*)malloc(1024*1024*sizeof(char));


  for (i_file = 0; i_file < nchunks; ++i_file) {
    //write data from local files (*fp) to output file (*out)
    FILE *fp;
    fp = fopen(chunks_filenames[i_file], "r");
    
    //transmit by chunks of 1 mb
    int size;
    do {
      size = fread(buffer, 1, sizeof(buffer), fp);
      if (size <= 0) break;
      fwrite(buffer, 1, size, out);
    } while (size == sizeof(buffer));
    //reached EOF, close local file
    fclose(fp);
    //printf("  Data from %s has been sucessfully transmitted to %s\n", chunks_filenames[i_file], output_filename);
 }
 fclose(out);
}
