#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //measure of parallel time evaluation
#include "init.h"
#include "integral.h"
#include "test_function.h"
#include "weight.h"

#define FILENAME_MAXSIZE 128

/*
 * This program uses analytic expression of the test-function and of the weight in 3D
 * and computes its weighted ray transforms in 3D. The computations are performed
 * slice-by-slice in direciton of OZ. In each 'slice' the set of rays corresponds to
 * 'uniform geometry' in the plane (uniform angles, uniform shifts).
 *
 * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
 * centered unit ball. Otherwise the generated data will be incomplete and this will produce artifacts in
 * reconstructions.
 *
 * The parameters for the set of rays are taken from configuration file (-p filename).
 *
 * The result of computations is stored in a separate csv file (-o filename)
*/

/* ***************************************************************************
 * *****************  HANDLING FILES (README, AGGREGATION)********************
 * ***************************************************************************
 */

void generate_data_readme(char *func_filename, char *weight_filename, int Nz, int nshift, int nphi) {
  FILE *doc_file;
  char doc_filename [ FILENAME_MAXSIZE ];

  sprintf(doc_filename, "rdme_%s", func_filename);

  doc_file = fopen(doc_filename, "w");
  fprintf(doc_file, "  Weighted ray transforms in 3D for the test-function in %s, weight in %s\n", func_filename, weight_filename);
  fprintf(doc_file, "  Format of the data: [z], [s], [phi], [Pf]\n");
  fprintf(doc_file, "  Parameters of the data:\n");
  fprintf(doc_file, "    number of slices in OZ direction [Nz]  : %d,\n", Nz);
  fprintf(doc_file, "    number of shifts of [s]    : %d.\n", nshift);
  fprintf(doc_file, "    number of angles of [phi]  : %d,\n", nphi);

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

    //transmit data by chunks of 1 mb
    char buffer[1024 * 1024];
    int size;
    do {
      size = fread(buffer, 1, sizeof(buffer), fp);
      if (size <= 0) break;
      fwrite(buffer, 1, size, out);
    } while (size == sizeof(buffer));
    //reached EOF -> close chunk file

    fclose(fp);
 }
}

/* ***************************************************************************
 * *************************************USAGE ********************************
 * ***************************************************************************
 */


const char* program_name;

/* Prints usage information to STREAM (stdout or stderr), and
 * exits the program with EXIT_CODE. */

void print_usage(FILE* stream, int exit_code) {
  fprintf(stream, "Usage: %s -p (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream,
	  "   -h --help                 Display this usage information.\n"
	  "   -p --parameters filename  Read parameters of the sampling grid for ray transforms from file.\n"
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

  /* Reading options */
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
	printf("Output file with ray data: %s\n", optarg);
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

  int nshift = parameters[0];
  int nphi   = parameters[1];
  int ngrid  = parameters[2];		// number of z-slices; (affects density of points for and integral)

  /* Init grid in Radon space */
  double *phi, *shift;
  printf("  Initializing the grid in ray transforms space...");
  if (init_ray_grid_alloc(nshift, nphi, &shift, &phi) != 0) {		// allocate memory for arrays : phi, shift
     perror("Aborting : init_ray_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (init_ray_grid(nshift, nphi, shift, phi) != 0) {			// set values of : phi, shift
     perror("Aborting : init_ray_grid_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");

  /************************************** END OF INITIALIZATIONS ******************************************/


  /************************* START RADON TRANSFORM (PARALLELIZATION VIA OPENMP) ************************/
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

        //set equally parts of OZ slices for each thread
	int block_size = ngrid / nthreads + (thread_ID < (ngrid % nthreads) ? 1 : 0);

        int z_shift_min = (thread_ID < (ngrid % nthreads) ? (ngrid / nthreads + 1)*thread_ID : (ngrid / nthreads)*thread_ID + (ngrid % nthreads)),
	    z_shift_max = z_shift_min + block_size;

        printf("    Thread %d. Domain of work (z slices): start %d, end %d, size %d\n", thread_ID, z_shift_min, z_shift_max-1, block_size);


        //open chunk ouptut file
        FILE *foutput;
	foutput = fopen(chunks_filenames[thread_ID], "w");

	int i_slice, i_phi, i_shift;
	double zslice;

        //ray integrals in the given OZ slice
        for (i_slice = z_shift_min; i_slice < z_shift_max; ++i_slice) {
	    zslice = (-1.0) + i_slice * 2.0 / (ngrid - 1);

	    //iterate ray integrals over shifts in slice plane
	    for (i_shift = 0; i_shift < nshift; ++i_shift) {

		//iterate ray integrals over polar angle in slice plane
		for (i_phi = 0; i_phi < nphi; ++i_phi) {

		    //ray integral (z_slice -- z coordinate of the slice, (phi, shift) -- coordinates of the ray in OZ slice, ngrid -- density of the points)
		    double rt = ray_integral3d(&test_function,
					       &weight_function,
					         zslice,
						 shift[i_shift],
						 phi[i_phi],
						 ngrid);

		    //write result to file
		    fprintf(foutput, "%lf, %lf, %lf, %lf\n", zslice, shift[i_shift], phi[i_phi], rt);
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

  generate_data_readme(output_filename, ngrid, nshift, nphi);
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
  if (clean_ray_grid(shift, phi) != 0) {
     fprintf(stderr, "Warning : clean_ray_grid.\n");
  }
  printf("Done.\n");

  /******************************************* END CLEANING MEMORY ****************************************/

  exit(EXIT_SUCCESS);
}
