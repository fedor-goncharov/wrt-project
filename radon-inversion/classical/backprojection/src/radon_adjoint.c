#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> // measurement of time of evaluation
#include "init.h"
#include "integral.h"

#define FILENAME_MAXSIZE 128
/*
 * The program reads the grid data for Radon transforms in 3D from file (-i filename)
 * and computes the adjoint Radon transform on the grid in the unit cube [-1,1]^3.
 * Parameters of grids in Radon space and in the unit cube are taken from the 
 * configuration file (-p filename).
 * 
 * The results of computations are stored in a separate CSV file (-o filename). 
*/

/* ***************************************************************************
 * ****************************** HANDLING FILES *****************************
 * ***************************************************************************
 */
void generate_readme(char* input_filename, int nphi, int ntheta, int nshift, 
                     char* output_filename, int ngrid) {
  
  FILE *readme_file;
  char readme_filename [ FILENAME_MAXSIZE ];
  
  sprintf(readme_filename, "readme_%s", output_filename);

  readme_file = fopen(readme_filename, "w");
  fprintf(readme_file, "  Radon transforms in 3D were taken in %s\n", input_filename);
  fprintf(readme_file, "  Parameters of the input data:\n");
  fprintf(readme_file, "     number of angles of [phi]  : %d\n", nphi);
  fprintf(readme_file, "     number of angles of [theta]: %d\n", ntheta);
  fprintf(readme_file, "     number of shifts of [s]    : %d\n", nshift);
  fprintf(readme_file, "  Output of adjoint transform was stored in %s\n", output_filename);
  fprintf(readme_file, "  Parameters of the output data:\n");
  fprintf(readme_file, "     number of points on the grid: %d\n", ngrid);

  fclose(readme_file);
}

void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks) {
  
  FILE *out; //output file
  out = fopen(output_filename, "w");
  
  int i_file;
  for (i_file = 0; i_file < nchunks; ++i_file) {
    //write data from local files (*fp) to output file (*out)
    FILE *fp;
    fp = fopen(chunks_filenames[i_file], "r");
    
    //transmit by chunks of 1 mb
    char buffer[1024 * 1024];
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
}

/* ***************************************************************************
 * ****************************** HANDLING INPUT *****************************
 * ***************************************************************************
 */


const char* program_name;

/* Prints usage information to STREAM (stdout or stderr), and 
 * exits the program with EXIT_CODE. Does not return. */

void print_usage(FILE* stream, int exit_code) {
  fprintf(stream, "Usage: %s -p (filename) -i (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream, 
	  "   -h --help                 Display this usage information.\n"
	  "   -p --parameters filename  Read parameters of grids for Radon space and for the test-function.\n"
	  "   -i --input filename       Read Radon data from file.\n"
	  "   -o --output filename      Write adjoint transform to the file.\n"
	  "   -n --nthreads number      Number of OpenMP threads for parallelization.\n");
  exit(exit_code);
}

/* ***************************************************************************
 * ****************************** ENTRY POINT ********************************
 * ***************************************************************************
 */

int main(int argc, char * argv[]) {
  
  if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help'.\n", argv[0]);
     exit(EXIT_FAILURE);
  }
  
  int next_option;
  const char* short_options = "hp:i:o:n:";
  const struct option long_options[] = {
    {"help", no_argument, 0, 0}, 
    {"parameters", required_argument, 0, 0},
    {"input", required_argument, 0, 0}, 
    {"output", required_argument, 0, 0},
    {"nthreads", required_argument, 0, 0},
    {NULL, 0, NULL, 0}				/* Required at the end of array */
  };

  char* parameters_filename;
  char* input_filename; 
  char* output_filename;
  int nthreads = 1;
  
  program_name = argv[0];

  /* Reading long options */ 
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);
    
    switch (next_option) {
      case 'h':
        print_usage(stdout, 0);
      case 'p':
	      parameters_filename = optarg;
	      assert(optarg[0] != '-');
	      printf("Parameters file: %s\n", optarg);
	      break;
      case 'i':
	      input_filename = optarg;
	      assert(optarg[0] != '-');
	      printf("Input file with test-function data: %s\n", optarg);
	      break;
      case 'o':
	      output_filename = optarg;
	      assert(optarg[0] != '-');
	      printf("Output file with Radon data: %s\n", optarg);
	      break;
      case 'n':
	      nthreads = atoi(optarg);
	      assert(nthreads > 0);
	      printf("Number of threads: %d\n", nthreads);
	      break;
      case ':': /* missing argument */ 
	      fprintf(stderr, "%s: option '-%c' requires an argument\n", argv[0], optopt);
	      print_usage(stderr, 1);
      case '?':
	      print_usage(stderr, 1);
      case -1: /* Done with the options. */
	      break;
      default:
	      abort();
    }
  } while (next_option != -1);
  
  if (optind < argc) {
    printf("Warning: non-option ARGV-arguments are omitted. See '-h' or '--help'.\n"); 
  }
  
  /**************************************** START INITIALIZATIONS *****************************************/
  /* Init run parameters, grids, values of Radon transforms */

  printf("  Reading parameters from %s...", parameters_filename);
  int *parameters;
  if (parameters_alloc(&parameters) != 0) {
     perror("Aborting : init_parameters_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (parameters_init(parameters_filename, parameters) != 0) {
     perror("Aborting : init_parameters.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  int nphi   = parameters[0];
  int ntheta = parameters[1];
  int nshift = parameters[2];
  int ngrid  = parameters[3];
  
  double *phi, *theta, *shift;  // Radon grid (equatorial angle, latitude angle, shift)
  double *theta_weights;
  double ***radon_values;       // values of Radon transform
  
  /* Init grid in Radon space */
  printf("  Reading Radon transforms...");
  if (radon_all_alloc(nphi, ntheta, nshift, &phi, &theta, &theta_weights, &shift, &radon_values) != 0) {
     perror("Aborting : radon_all_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (radon_grid_init(nphi, ntheta, nshift, phi, theta, theta_weights, shift) != 0) {
     perror("Aborting : radon_grid_init.\n");
     exit(EXIT_FAILURE);
  }
  if (radon_values_init(input_filename, radon_values, nphi, ntheta, nshift) != 0) {
     perror("Aborting : radon_values_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");
  
  /************************************** END OF INITIALIZATIONS ******************************************/
  
  
  /************************* START RADON TRANSFORM WITH PARALLELIZATION VIA OPENMP ************************/
  //start timer
  struct timeval start, end;
  gettimeofday(&start, NULL);  

  
  /* Init chunks filenames (nchunks = nthreads). Chunk filename %d_(output_filename) */
  char** chunks_filenames;
  chunks_filenames_alloc(output_filename, nthreads, &chunks_filenames, FILENAME_MAXSIZE);
  chunks_filenames_init(output_filename, nthreads, chunks_filenames);
  
  omp_set_num_threads(nthreads);
  #pragma omp parallel
  {
        int thread_ID = omp_get_thread_num();

        //set equally parts of shifts for each thread
        int block_size = ngrid / nthreads + (thread_ID < (ngrid % nthreads) ? 1 : 0);
	
        int i_z_min = (thread_ID < (ngrid % nthreads) ? (ngrid / nthreads + 1) * thread_ID : (ngrid / nthreads)*thread_ID + (ngrid % nthreads)),
	          i_z_max = i_z_min + block_size;
	    
        printf("    Thread %d. Domain of work (shifts): start %d, end %d, size %d\n", thread_ID, i_z_min, i_z_max-1, block_size);
        
	
        //open ouptut file
        FILE *output_file;
        output_file = fopen(chunks_filenames[thread_ID], "w");
	      int i_z, i_x, i_y;
    
        // z-coordinate
        for (i_z = i_z_min; i_z < i_z_max; ++i_z) {      
      
            // x-coordinate
            for (i_x = 0; i_x < ngrid; ++i_x) {
	
	             // y-coordinate
	             for (i_y = 0; i_y < ngrid; ++i_y) {
		             
                 double x = (-1.0) + (1.0 / ngrid) + i_x * (2.0 / ngrid),
                        y = (-1.0) + (1.0 / ngrid) + i_y * (2.0 / ngrid),
                        z = (-1.0) + (1.0 / ngrid) + i_z * (2.0 / ngrid);

                 //adjoint Radon transform at point (x,y,z)
		             double adjoint = radon_adjoint(radon_values, x, y, z, phi, theta, theta_weights, shift, nphi, ntheta, nshift);
		             fprintf(output_file, "%lf\n", adjoint);
               }
            }
        }
	      fclose(output_file);
	      printf("    Thread %d. Job is done.\n", thread_ID);
  } //threads finish
  
  //stop timer
  gettimeofday(&end, NULL);
  double delta_t = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;
  printf("\n");
  printf("  Number of threads: %d, elapsed time: %3.1f sec\n", nthreads, delta_t);
  
  /************************* END RADON TRANSFORM WITH PARALLELIZATION VIA OPENMP **************************/
  
  
  /*********************************** START AGGREGATING DATA TO ONE FILE *********************************/
  
  printf("  Generating readme...");
  generate_readme(input_filename, nphi, ntheta, nshift, output_filename, ngrid);
  printf("  Done.\n");
  printf("  Aggregating chunks to one output file...");
  chunks_aggregate(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");
  
  /*********************************** END AGGREGATING DATA TO ONE FILE ***********************************/
  

  /******************************************* START CLEANING MEMORY **************************************/
  
  printf("  Cleaning allocated memory...");
  
  if (chunks_files_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : chunks_files_clean.\n");
  }
  if (chunks_filenames_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : chunks_filenames_clean.\n");
  }
 
  if (parameters_clean(parameters) != 0) {
     fprintf(stderr, "Warning : parameters_clean.\n");
  }

  if (radon_values_clean(radon_values, nphi, ntheta, nshift) != 0) {
     fprintf(stderr, "Warning : radon_values_clean.\n");
  }
  if (radon_grid_clean(phi, theta, theta_weights, shift) != 0) {
     fprintf(stderr, "Warning : radon_grid_clean.\n");
  }
  printf("Done.\n");
  
  /******************************************* END CLEANING MEMORY ****************************************/
  
  exit(EXIT_SUCCESS);
}
