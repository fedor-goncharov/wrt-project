#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //measure physical time of evaluation

#include "init.h"
#include "integral.h"

#define FILENAME_MAXSIZE 128

/*
 * This program computes ray transforms of a test-function in 3D slice-by-slice. The test-function
 * is given by its values on the uniform grid in [-1,1]^3 in a separate binary file.
 * 
 * The set of rays along which the computations are performed corresponds to
 * the parallel scanner geometry : that is the support of the test-function is sliced
 * by a set of parallel planes, which are parallel to the XY plane. Planes are parametrized by
 * z-coordinate in [-1, 1].
 *
 * In each plane rays are parametrized by (s, phi), where s is in [-1,1], phi in (0,2*pi).
 *
 * The grids of s, phi are uniform in their intervals. For integration along a line the
 * Simpson's rule is used.
 *
 * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
 * centered ball of radius R cm.
 *
 * 
 *
 * The result of computations is stored in a separate binary file (-o filename) in the following order:
 *
 *   for slice = (0 : nslice-1) // z-slice
 *      for shift = (0 : nshift-1) // value for s for rays in the plane
 *         for phi = (0 : nphi-1) // value for phi for rays in the plane
 *
 *             fwrite(output, ray_transform(f, "%lf\n", slice, shift, phi)
 *
 *         endfor
 *      endfor
 *   endfor
 *
 * REMARK : In all utilities using ray transforms - reading the data will be exactly in this order.
 *
 *
*/


void readme_generate(char* filename, int nslices, int nshift, int nphi) {
  FILE *file;
  char readme_filename [ FILENAME_MAXSIZE ];

  sprintf(readme_filename, "readme_%s", filename);

  file = fopen(readme_filename, "w");
  fprintf(file, "  Ray transforms in 3D for test-function given in file %s\n", filename);
  fprintf(file, "  Format of the data: [z], [shift], [phi], [Pf]\n");
  fprintf(file, "  Parameters of the data:\n");
  fprintf(file, "      number of z-slicez OZ direction  [nslices]   : %d,\n", nslices);
  fprintf(file, "      number of angles on circle       [phi]  : %d,\n", nphi);
  fprintf(file, "      number of shifts per direciton   [s]    : %d.\n", nshift);

  fclose(file);
}

/* 
 * ------------------------- HANDLING INPUT --------------------------------------
 * 
 */


const char* program_name;

/* Prints usage information to STREAM (stdout or stderr), and
 * exits the program with EXIT_CODE. Does not return. */

void usage_print(FILE* stream, int exit_code) {
  fprintf(stream, "Usage: %s -p (filename) -i (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream,
	  "   -h --help                 Display usage information.\n"
	  "   -p --parameters filename  Read parameters of the test-function and 3D ray transforms from file.\n"
	  "   -i --input filename       Read test-functiom data from file (according to parameters in '-p').\n"
	  "   -o --output filename      Place output data to a file.\n"
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
  const char* short_options = "hp:i:o:n:";
  const struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"parameters", required_argument, 0, 0}, // deprecated argument - parameters should be taken as an argument
    {"input", required_argument, NULL, 'i'},
    {"output", required_argument, NULL, 'o'},
    {"nthreads", required_argument, NULL, 'n'},
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
	      usage_print(stdout, 0);
    case 'p': // deprecated argument
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
	      usage_print(stderr, 1);
      case '?':
	      usage_print(stderr, 1);
      case -1: /* Done with the options. */
	      break;
      default:
        printf("?? getopt returned character code 0%o??\n", next_option);
        exit(EXIT_SUCCESS);
    }
  } while (next_option != -1);

  if (optind < argc) {
    printf("Warning: non-option ARGV-arguments are omitted. See '-h' or '--help'.\n");
  }

  /**************************************** START INITIALIZATIONS *****************************************/
  /* Init configuration parameters of the grid and of Radon data from file */
  printf("  Reading parameters from %s...", parameters_filename);
  int *parameters;
  if (parameters_alloc(&parameters) != 0) {
     perror("Aborting : parameters_alloc.\n");
     exit(EXIT_FAILURE);
  }

  if (parameters_init(parameters_filename, parameters) != 0) {
    perror("Aborting : parameters_init.\n");
    exit(EXIT_FAILURE);
  }
  printf("Done.\n");

  int ngrid       = parameters[0];	// number of points per dimension (in [-1,1]) in the grid for the test-function
  int nslices     = parameters[1];  // number of z slices
  int nshift      = parameters[2];
  int nphi        = parameters[3];


  /* Init values of the test function on the grid */
  printf("  Initializing test-function data...");
  double*** values;
  if (tfunction_values_alloc(ngrid, &values) != 0) {
     perror("Aborting : tfunction_values_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (tfunction_values_init(input_filename, ngrid, values) != 0) {
     perror("Aborting : tfunction_values_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");


  /* Init grid in Radon space */
  double *phi, *shift;
  printf("  Initializing the grid in ray transforms space...");
  if (ray_grid_alloc(nshift, nphi, &shift, &phi) != 0) {
     perror("Aborting : ray_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (ray_grid_init(nshift, nphi, shift, phi) != 0 ) {
     perror("Aborting : ray_grid_init.\n");
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
	      int block_size = nslices / nthreads + (thread_ID < (nslices % nthreads) ? 1 : 0);

        int i_slice_min = (thread_ID < (nslices % nthreads) ? (nslices / nthreads + 1) * thread_ID : (nslices / nthreads)*thread_ID + (nslices % nthreads)),
	          i_slice_max = i_slice_min + block_size;

        printf("    Thread %d. Domain of work (z slices): start %d, end %d, size %d\n", thread_ID, i_slice_min, i_slice_max-1, block_size);


        //open chunk ouptut file
        FILE *foutput;
        foutput = fopen(chunks_filenames[thread_ID], "w");
	      int i_zslice, i_phi, i_shift;

	      double zslice;
        const double dzslice = 2.0 / (nslices - 1);
        //iterate ray integrals over z-slices
        for (i_zslice = i_slice_min; i_zslice < i_slice_max; ++i_zslice) {
	          zslice = (-1.0) + i_zslice * dzslice;

            for (i_shift = 0; i_shift < nshift; ++i_shift) {
                //iterate ray integrals over longitude angle
                for (i_phi = 0; i_phi < nphi; ++i_phi) {
	              //ray integral
		              double val = ray_transform(values, ngrid,
						                                zslice, shift[i_shift], phi[i_phi]);
		              fprintf(foutput, "%lf, %lf, %lf, %lf\n", zslice, shift[i_shift], phi[i_phi], val);
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

  printf("  Generating readme...");
  readme_generate(output_filename, nslices, nshift, nphi);
  printf("  Done.\n");
  printf("  Aggregating chunks to one output file...");
  chunks_aggregate(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");

  /*********************************** END AGGREGATING DATA TO ONE FILE ***********************************/


  /******************************************* START CLEANING MEMORY **************************************/

  printf("  Cleaning allocated memory...");

  // clean files
  if (chunks_files_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : chunks_files_clean.\n");
  }
  if (chunks_filenames_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : chunks_filenames_clean.\n");
  }
  // clean memory
  if (parameters_clean(parameters) != 0) {
     fprintf(stderr, "Warning : parameters_clean.\n");
  }
  if (tfunction_values_clean(ngrid, values) != 0) {		//clean data for the test-function
     fprintf(stderr, "Warning : tfunction_values_clean.\n");
  }
  if (ray_grid_clean(shift, phi) != 0) {
     fprintf(stderr, "Warning : ray_grid_clean.\n");
  }

  printf("Done.\n");

  /******************************************* END CLEANING MEMORY ****************************************/

  exit(EXIT_SUCCESS);
}
