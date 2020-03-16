#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdbool.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //measure of parallel time evaluation
#include "init.h"
#include "integral.h"
#include "test_function.h"
#include "weight.h"

#define FILENAME_MAXSIZE 64

/*
 * Program computes weighted Radon transforms along 2d-planes of a test function given by an analytical expression 
 * and for a weight also given by an analytical expression.
 * The set of planes corresponds to 'Gauss - uniform' sampling scheme, where the parameters of this scheme 
 * are given a separate file (-p filename). 
 * 
 * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
 * centered unit ball. Otherwise the generated data will be incomplete and this will produce artifacts in reconstructions.
 *
 * The result is stored in a separate CSV file (-o filename)
*/

/* ***************************************************************************
 * ****************************** HANDLING FILES *****************************
 * ***************************************************************************
 */

void generate_data_readme(char* func_filename, char* weight_filename, int nphi, int ntheta, int nshift) {
  FILE *doc_file;
  char doc_filename [ FILENAME_MAXSIZE ];
  
  sprintf(doc_filename, "rdme_%s", func_filename);

  doc_file = fopen(doc_filename, "w");
  fprintf(doc_file, "  Weighted Radon transforms in 3D for test-function in %s and weight function in %s\n", func_filename, weight_filename);
  fprintf(doc_file, "  Format of the data: [s], [phi], [theta], [RWf]\n");
  fprintf(doc_file, "  Parameters of the data:\n");
  fprintf(doc_file, "    number of angles of [phi]  : %d,\n", nphi);
  fprintf(doc_file, "    number of angles of [theta]: %d,\n", ntheta);
  fprintf(doc_file, "    number of shifts of [s]    : %d.\n", nshift);
    
  fclose(doc_file);
}

void aggregate_data_onefile(char* output_filename, char ** chunks_filenames, int nchunks) {
  
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
  fprintf(stream, "Usage: %s -p (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream, 
	  "   -h --help                 Display this usage information.\n"
	  "   -p --parameters filename  Read parameters of the grid test-function and Radon transforms from file.\n"
	  "   -o --output filename      Write output data to the file.\n"
	  "   -n --nthreads number      Number of OpenMP threads for parallelization of calculations.\n");
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
  const char* short_options = "hp:o:n:";
  const struct option long_options[] = {
    {"help", no_argument, 0, 0}, 
    {"parameters", required_argument, 0, 0},
    {"output", required_argument, 0, 0},
    {"num_threads", required_argument, 0, 0},
    {NULL, 0, NULL, 0}				/* Required at the end of array */
  };

  char* parameters_filename;
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
	fprintf(stderr, "%s: option '-%c' requires an argument\n",
		argv[0], optopt);
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
  /* Init configuration parameters of the grid and of Radon data from file */
  printf("  Reading parameters from %s...", parameters_filename);
  int *parameters;
  if (init_parameters_alloc(&parameters) != 0) {
     perror("Aborting : init_parameters_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (init_parameters(parameters_filename, parameters) != 0) {
     perror("Aborting : init_parameters.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  int nphi   = parameters[0];
  int ntheta = parameters[1];
  int nshift = parameters[2];
  int ngrid  = parameters[3];
  
  /* Init grid in Radon space */
  double *phi, *theta, *shift;
  printf("  Initializing the grid in Radon space...");
  if (init_radon_grid_alloc(nphi, ntheta, nshift, &phi, &theta, &shift) != 0) {
     perror("Aborting : init_radon_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (init_radon_grid(nphi, ntheta, nshift, phi, theta, shift) != 0) {
     perror("Aborting : init_radon_grid.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");
  
  /************************************** END OF INITIALIZATIONS ******************************************/
  
  
  /************************* START WEIGHTED RADON TRANSFORM WITH PARALLELIZATION VIA OPENMP ************************/
  //start timer
  struct timeval start, end;
  gettimeofday(&start, NULL);  

  
  /* Init chunks filenames (nchunks = nthreads). Chunk filename %d_(output_filename) */
  char** chunks_filenames;
  init_chunks_filenames_alloc(output_filename, nthreads, &chunks_filenames, FILENAME_MAXSIZE);
  init_chunks_filenames(output_filename, nthreads, chunks_filenames);
  
  omp_set_num_threads(nthreads);
  #pragma omp parallel
  {
        int thread_ID = omp_get_thread_num();

        //set equally parts of shifts for each thread
	int block_size = nshift / nthreads + (thread_ID < (nshift % nthreads) ? 1 : 0);
	
        int i_shift_min = (thread_ID < (nshift % nthreads) ? (nshift / nthreads + 1)*thread_ID : (nshift / nthreads)*thread_ID + (nshift % nthreads)),
	    i_shift_max = i_shift_min + block_size;
	    
        printf("    Thread %d. Domain of work (shifts): start %d, end %d, size %d\n", thread_ID, i_shift_min, i_shift_max-1, block_size);
        
	
        //open ouptut file
        FILE *foutput;
	foutput = fopen(chunks_filenames[thread_ID], "w");
	int i_phi, i_theta, i_shift;
    
        //iterate plane integrals over shifts along normals
        for (i_shift = i_shift_min; i_shift < i_shift_max; ++i_shift) {      
      
            //iterate plane integrals over longitude angle
            for (i_phi = 0; i_phi < nphi; ++i_phi) {
	
	        //iterate plane integrals over latitude angle
	        for (i_theta = 0; i_theta < ntheta; ++i_theta) {
		    //plane integral
		    double rt = plane_integral3d(&test_function, &weight_function, ngrid, 
						 phi[i_phi], theta[i_theta], shift[i_shift]);
		    //write it down to file
		    fprintf(foutput, "%lf, %lf, %lf, %lf\n", shift[i_shift], phi[i_phi], theta[i_theta], rt);
                }
            }
        }
	fclose(foutput);
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
  
  generate_data_readme(output_filename, nphi, ntheta, nshift);
  printf("  Aggregating chunks to one output file...");
  aggregate_data_onefile(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");
  
  /*********************************** END AGGREGATING DATA TO ONE FILE ***********************************/
  

  /******************************************* START CLEANING MEMORY **************************************/
  
  printf("  Cleaning allocated memory...");
  
  if (clean_chunks_files(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_files.\n");
  }
  if (clean_chunks_filenames(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_filenames.\n");
  }
  if (clean_parameters(parameters) != 0) {
     fprintf(stderr, "Warning : clean_parameters.\n");
  }
  if (clean_radon_grid(phi, theta, shift) != 0) {
     fprintf(stderr, "Warning : clean_radon_grid.\n");
  }
  printf("Done.\n");
  
  /******************************************* END CLEANING MEMORY ****************************************/
  
  exit(EXIT_SUCCESS);
  
}
