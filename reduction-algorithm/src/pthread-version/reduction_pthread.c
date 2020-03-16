#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>   // input options handling
#include <assert.h>
#include <pthread.h>
#include <sys/time.h> // measure of parallel time evaluation

#include "init.h"       // readme generation, memory allocation/deallocation, aggregation of results in one file
#include "threadfunc.h" // parallelization via POSIX threading (includes headers with math)


#define FILENAME_MAXSIZE 128

/*
 *
 * Version of the reduction program which uses POSIX threads for parallelization than OpenMP threds
 *
*/

/*
  * Program performes reduction of data given by (weighted) ray transforms in 3D
  * to the data given by (weighted) Radon transforms along 2d-planes in 3D. Such
  * reduction is a consequence of a representation of a plane integral through integration along
  * a set of parallel lines (even for weighted ray-Radon transforms); see [Goncharov, Novikov, Inverse Problems 2017].
  * Parameters of ray transforms on the input and of Radon transforms on the output are
  * stored in a configuration file (-p filename).
  *
  * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
  * centered unit ball. Otherwise the generated data will be incomplete and this will produce artifacts in
  * reconstructions.
  *
  * Output is a separate file in CSV format.
  */

const char* program_name;

/* Prints usage information to STREAM (stdout or stderr), and
 * exits the program with EXIT_CODE. */

int main(int argc, char ** argv) {

  //------------------------------------ARGUMENTS HANDLING------------------------------------------------------------//
  if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help'.\n", argv[0]);
     exit(EXIT_FAILURE);
  }

  program_name = argv[0];

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

  //------------------------------ Reading options-------------------------------//
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);

    switch (next_option) {
      case 'h':
	print_usage(stdout, 0);
      case 'p':
	parameters_filename = optarg;
	assert(optarg[0] != '-');
	printf("Parameters file: %s\n", parameters_filename);
	break;
      case 'i':
	input_filename = optarg;
	assert(optarg[0] != '-');
	printf("Input file with ray data: %s\n", input_filename);
	break;
      case 'o':
	output_filename = optarg;
	assert(optarg[0] != '-');
	printf("Output file with Radon data: %s\n", output_filename);
	break;
      case 'n':
	nthreads = atoi(optarg);
	assert(nthreads > 0);
	printf("Number of threads: %d\n", nthreads);
	break;
      case ':': /* missing argument */
	fprintf(stderr, "%s: option '-%c' requires an argument\n",
		argv[0], optopt);
	print_usage(stderr, EXIT_FAILURE);
      case '?':
	print_usage(stderr, EXIT_FAILURE);
      case -1: /* Done with the options. */
	break;
      default:
	abort();
    }
  } while (next_option != -1);

  if (optind < argc) {
    printf("Warning: non-option ARGV-arguments are omitted. See '-h' or '--help'.\n");
  }

  //-------------------------------------- INITIALIZATIONS ------------------------------------------------------------------------------/
  //---------------------- Read parameters of the ray/Radon grids from parameters file --------------------------------------------------/

  printf("  Reading parameters from %s...", parameters_filename);
  int *parameters;
  if (parameters_alloc(&parameters) != 0) { //allocate memory for parameters
     perror("Aborting : init_parameters_alloc.\n");
     exit(EXIT_FAILURE);
  }

  if (parameters_init(parameters_filename, parameters) != 0) {
     perror("Aborting : init_parameters.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");

  int nshift = parameters[0];
  int nphi   = parameters[1];
  int ntheta = parameters[2];
  int ngrid  = parameters[3];		// number of points on grid in one dimension


  // --------------------------------------------Read ray transform data from file --------------------------------------------------------/

  printf("  Initializing data for ray transforms...");
  double*** ray_transforms_values;
  if (ray_values_alloc(nshift, nphi, ngrid, &ray_transforms_values) != 0) {
     perror("Aborting : ray_values_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (ray_values_init(input_filename, nshift, nphi, ngrid, ray_transforms_values) != 0) {
     perror("Aborting : ray_values_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");

  //---------------------------------------------------------------------------------------------------------------------------------------/

  //------------------------------------------------------ init grid in Radon space -------------------------------------------------------/
  double *shift, *phi, *theta;
  printf("  Initializing the grid in Radon space...");
  if (radon_grid_alloc(nshift, nphi, ntheta, &shift, &phi, &theta) != 0) {		// allocate memory for arrays : shift, theta, phi
     perror("Aborting : radon_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (radon_grid_init(nshift, nphi, ntheta, shift, phi, theta) != 0) {		        // set values of : shift, theta, phi
     perror("Aborting : radon_grid_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");

  //------------------------------------ END OF INITIALIZATIONS ------------------------------------------------------------------------/

  //------------------------ START REDUCTION RAY TRANSFORMS -> RADON TRANSFORMS (PARALLELIZATION VIA POSIX THREADS) --------------------/
  //start timer
  struct timeval start, end;
  gettimeofday(&start, NULL);


  /* Init chunks filenames (nchunks = nthreads). Chunk filename %d_(output_filename) */
  char** chunks_filenames;
  chunks_filenames_alloc(output_filename, nthreads, &chunks_filenames, FILENAME_MAXSIZE);
  chunks_filenames_init(output_filename, nthreads, chunks_filenames);

  pthread_t thread[nthreads];
  struct reduction_threadf_argument thread_arg[nthreads];
  int i_thread;

  for (i_thread = 0; i_thread < nthreads; ++i_thread) {

    //set equally shifts of planes for each thread
    int block_size = nshift / nthreads + (i_thread < (nshift % nthreads) ? 1 : 0);
    int shift_min = (i_thread < (nshift % nthreads) ? (nshift / nthreads + 1)*i_thread : (nshift / nthreads)*i_thread + (nshift % nthreads)),
        shift_max = shift_min + block_size;

    thread_arg[i_thread].thread_ID = i_thread;
    thread_arg[i_thread].output_filename = chunks_filenames[i_thread];
    thread_arg[i_thread].shift_min = shift_min;
    thread_arg[i_thread].shift_max = shift_max;
    thread_arg[i_thread].nshift = nshift;
    thread_arg[i_thread].nphi = nphi;
    thread_arg[i_thread].ntheta = ntheta;
    thread_arg[i_thread].ngrid = ngrid;

    thread_arg[i_thread].ray_transforms_values = ray_transforms_values;
    thread_arg[i_thread].shift = shift;
    thread_arg[i_thread].phi = phi;
    thread_arg[i_thread].theta = theta;

    // run thread for computation of R-T over part of planes
    pthread_create(&(thread[i_thread]), NULL, &reduction_ray_radon_threadf, &(thread_arg[i_thread]));
  }

  // wait untill all threads are done
  for (i_thread = 0; i_thread < nthreads; ++i_thread) {
    pthread_join(thread[i_thread], NULL);
  }

  //stop timer
  gettimeofday(&end, NULL);
  double delta_t = ((end.tv_sec  - start.tv_sec) * 1e6u +
         end.tv_usec - start.tv_usec) / 1.e6;
  printf("\n");
  printf("  Number of threads: %d, elapsed time: %3.1f sec\n", nthreads, delta_t);

  //------------------------ END REDUCTION RAY TRANSFORMS -> RADON TRANSFORMS (PARALLELIZATION VIA POSIX)  ---------------------------/


  //------------------------------------------ START AGGREGATING DATA TO ONE FILE ----------------------------------------------------/

  generate_readme(output_filename, nshift, nphi, ntheta);
  printf("  Aggregating chunks to one output file...");
  aggregate_data_onefile(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");

  //---------------------------------------------- END AGGREGATING DATA TO ONE FILE --------------------------------------------------/


  //------------------------------------------------------- START CLEANING MEMORY ----------------------------------------------------/

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
  if (ray_values_clean(nshift, nphi, ngrid, ray_transforms_values) != 0) {
     fprintf(stderr, "Warning : ray_values_clean.\n");
  }
  if (radon_grid_clean(shift, phi, theta) != 0) {
     fprintf(stderr, "Warning : radon_grid_clean.\n");
  }
  printf("Done.\n");

  //----------------------------------------------------------- END CLEANING MEMORY --------------------------------------------------/

  exit(EXIT_SUCCESS);
}
