#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <gsl/gsl_integration.h>

#define MAX_FILENAME_LEN 128

/* read ray transforms */
// ---------------------------------------------------------------------------------------//
int ray_values_alloc(int nshift, int nphi, int ngrid, double**** ray_values) {
    int i, j;

    *ray_values = (double***)malloc(sizeof(double**) * ngrid);
    for (i = 0; i < ngrid; ++i) {
        (*ray_values)[i] = (double**)malloc(sizeof(double*) * nphi);
    }

    for (i = 0; i < ngrid; ++i) {
        for (j = 0; j < nphi; ++j) {
            (*ray_values)[i][j] = (double*)malloc(sizeof(double) * nshift);
        }
    }

    return 0;
}
int ray_values_init(char* filename, int nshift, int nphi, int ngrid, double*** ray_values) {

    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Recieved call %s while opening %s. \n", strerror(errno), filename);
        return errno;
    }
    int i_shift, i_phi, i_ngrid;
    int count; // bytes read
    int count_format = 1; // read one variable per line

    for (i_ngrid = 0; i_ngrid < ngrid; ++i_ngrid) {
        for (i_shift = 0; i_shift < nshift; ++i_shift) {
            for (i_phi = 0; i_phi < nphi; ++i_phi) {

                //read file with ray data line-by-line
                count = fscanf(f, "%lf\n", &(ray_values[i_ngrid][i_phi][i_shift]));
                
                if (count == EOF) {
                    if (ferror(f)) {
                        perror("Init test-function values : fscanf read value : EOF failure\n");
                    } else {
                        fprintf(stderr, "Error : init test-function : fscanf : reached the end of file before it was expected.\n");
                    }
		    fclose(f);
                    return -1;
                } else if (count != count_format) {
                    fprintf(stderr, "Error: init test-function : fscanf : matching failure, fscanf expected %d floating-point values.\n", count_format);
		    fclose(f);
                    return -1;
                }
            }
        }
    }
    fclose(f);
    return 0; //reading successfull
}
int ray_values_clean(int nshift, int nphi, int ngrid, double*** ray_values) {
    int i,j;
    for (i = 0; i < ngrid; ++i) {      //
        for (j = 0; j < nphi; ++j) {   //
            free(ray_values[i][j]);        // clean directions phi
        }
    }
    for (i = 0; i < ngrid; ++i) {    //
        free(ray_values[i]);         // clean shifts
    }
    free(ray_values);                    // clean z -- cross sections

    return 0;
}


/* initialize/clean Radon grid */
// ---------------------------------------------------------------------------------------//
int radon_grid_alloc(int nshift,  int nphi, int ntheta, double** shift, double** phi, double** theta) {
    *shift = (double*)malloc(nshift * sizeof(double));
    *phi =   (double*)malloc(nphi   * sizeof(double));
    *theta = (double*)malloc(ntheta * sizeof(double));

    return 0;
}
int radon_grid_init(int nshift, int nphi, int ntheta,  double* shift, double* phi, double* theta) {
    /*
     * phi --   uniform grid on S^1 = [0,2pi]
     * theta -- Gauss's angles on [0, pi]
     * shift -- uniform points on [-1, 1]
     */

    //init shift
    const double dshift = 2.0 / (nshift - 1);
    int i_shift;
    for (i_shift = 0; i_shift < nshift; ++i_shift) {
        shift[i_shift] = (-1.0) + i_shift * dshift;
    }

    //init phi
    const double dphi = (2 * M_PI) / nphi;
    int i_phi;
    for (i_phi = 0; i_phi < nphi; ++i_phi) {
        phi[i_phi] = dphi * i_phi;
    }

    //init theta (using GSL library to compute Gaussian points)
    int i_theta;
    gsl_integration_glfixed_table *gsl_quad_data = gsl_integration_glfixed_table_alloc(ntheta);
    for (i_theta = 0; i_theta < ntheta; ++i_theta) {
        double gauss_point, gauss_weight;
        gsl_integration_glfixed_point(-1.0, 1.0, i_theta, &gauss_point, &gauss_weight, gsl_quad_data);
        theta[i_theta] = acos(gauss_point);
    }

    gsl_integration_glfixed_table_free(gsl_quad_data);

    return 0;
}
int radon_grid_clean(double* shift, double* phi, double* theta) {
    if (shift != NULL)
        free(shift);
    if (phi != NULL)
        free(phi);
    if (theta != NULL)
        free(theta);

    return 0;
}

/* filenames utilities */
// ---------------------------------------------------------------------------------------//
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
        char buffer[MAX_FILENAME_LEN];
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

  //transmit by chunks by 10mb
  int buffer_size = 1024*1024*10;
  char* buffer = malloc(buffer_size);
    
  int i_file;
  for (i_file = 0; i_file < nchunks; ++i_file) {
    //write data from local files (*fp) to output file (*out)
    FILE *chunk_file;
    chunk_file = fopen(chunks_filenames[i_file], "r");

    int size;
    do {
      size = fread(buffer, 1, buffer_size, chunk_file);
      if (size <= 0) break;
      fwrite(buffer, sizeof(char), size, out);
    } while (size == sizeof(buffer));
    //reached EOF, close local file
    fclose(chunk_file);
 }

 free(buffer);
 fclose(out);
}

// --------------------------------- info utilities -----------------------------------------------//

void print_usage(FILE* stream, char* program_name) {
  fprintf(stream, "Usage: %s -p [phi] -s [shift] -t [theta] -z [slices] -i [input_file] -o (output_file) -n (number_threads)\n", program_name);
  fprintf(stream,
      "   -h --help                 Display usage information.\n"
      "   -p --nphi                 Number of projections per slice in ray data/phi angles in Radon data.\n" 
      "   -s --nshift               Number of shifts for ray and Radon data.\n"
      "   -t --ntheta               Number of theta angles for Radon data.\n"     
      "   -z --nslices              Number of z-slices in ray data.\n"
      "   -i --input filename       Read data given by ray transforms from a file.\n"
      "   -o --output filename      Write output data to the file.\n"
      "   -n --nthreads number      Number of OpenMP threads for parallelization.\n");
}


void generate_readme(char* input_filename, char* output_filename, int nphi, int ntheta, int nshift) {
  FILE *doc_file;
  char doc_filename [ MAX_FILENAME_LEN ];

  sprintf(doc_filename, "readme_%s", output_filename);

  doc_file = fopen(doc_filename, "w");
  fprintf(doc_file, "  Reduction of (weighted) ray transforms from file %s to (weighted) Radon transforms to file %s.\n", input_filename, output_filename);
  fprintf(doc_file, "  Parameters of the output data:\n");
  fprintf(doc_file, "      number of shifts    : %d\n", nshift);
  fprintf(doc_file, "      number of phi       : %d\n", nphi);
  fprintf(doc_file, "      number of theta     : %d\n", ntheta);
  fclose(doc_file);
}
