#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //measure physical time of evaluation

#include <unistd.h>  // supress output to stdout if needed
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "init.h"
#include "integral.h"

#define FILENAME_MAXSIZE 128

/*
 * This program computes ray transforms of a test-function in 3D in slice-by-slice manner. 
 * The test-function is given by its values on the uniform grid in [-R,R]^2 in z-slices 
 * in a separate binary file.
 * 
 * The set of rays along which the computations are performed corresponds to
 * the parallel scanner geometry : that is the support of the test-function is sliced
 * by a set of parallel planes, which are parallel to the XY plane. Planes are parametrized by
 * z-coordinate in [-z_min, z_max].
 *
 * In each plane rays are parametrized by (s, phi), where s is in [-R,R], phi in (0,2*pi).
 *
 * The grids of s, phi are uniform in their intervals. For integration along a line the
 * Simpson's rule is used.
 *
 * IMPORTANT : It is always assumed that on each slice the test-function is compactly supported with support in the
 * centered ball of radius R cm. If this is not respected, data will not correspond to the full ray transform.
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
 *
*/

// don't print verbose info
int supress_stdout() {
  fflush(stdout);

  int ret = dup(1);
  int nullfd = open("/dev/null", O_WRONLY);

  dup2(nullfd, 1);
  close(nullfd);

  return ret;
}

void readme_generate(char* input_filename, char* output_filename, int npixels, int nslices, int nphi, int nshift, double radius) {
  FILE *readme_file;
  char readme_filename [ FILENAME_MAXSIZE ];

  sprintf(readme_filename, "readme_%s", output_filename);

  readme_file = fopen(readme_filename, "w");
  fprintf(readme_file, "  Slice-by-slice ray transforms in 3D for test-function from file %s to file %s.\n", input_filename, output_filename);
  fprintf(readme_file, "  Parameters of the data:\n");
  fprintf(readme_file, "    number of pixels per dimension : %d,\n", npixels);
  fprintf(readme_file, "    number of z-slices : %d,\n", nslices);
  fprintf(readme_file, "    number of projections per slice : %d,\n", nphi);
  fprintf(readme_file, "    number of shifts per projection : %d.\n", nshift);
  fprintf(readme_file, "    radius of the support in cm : %f.\n", radius);

  fclose(readme_file);
}


const char* program_name;

/* Prints usage information to STREAM (stdout or stderr), and
 * exits the program with EXIT_CODE. Does not return. */

void usage_print(FILE* stream) {
  fprintf(stream, "Usage: %s -i (filename) -o (filename) -n (number) [ARGS]\n", program_name);
  fprintf(stream,
	  "   -h --help     no_arg       Display usage information.\n"
	  "   -p --nphi     integer      Number of projections per slice.\n"
	  "   -s --nshift   integer      Number of shifts per projection.\n"
	  "   -g --npixels  integer      Number of pixels per dimension.\n"
	  "   -z --nslices  integer      Number of slices in z-direction.\n"
	  "   -r --radius   float        Radius of the support in cm in a slice.\n"
	  "   -i --input    filename     Path to file with test-functiom data.\n"
	  "   -o --output   filename     Path to output file.\n"
	  "   -n --nthreads integer      Number of OpenMP threads for parallelization.\n"
	  "   -v --verbose  no_arg       Show extra information during computations.\n");
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
  const char* short_options = "hp:s:g:z:r:i:o:n:v";
  const struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"nphi", required_argument, NULL, 'p'}, 
    {"nshift", required_argument, NULL, 's'}, 
    {"npixels", required_argument, NULL, 'g'},
    {"nslices", required_argument, NULL, 'z'},
    {"radius", required_argument, NULL, 'r'},
    {"input", required_argument, NULL, 'i'},
    {"output", required_argument, NULL, 'o'},
    {"nthreads", required_argument, NULL, 'n'},
    {"verbose", no_argument, NULL, 'v'},
    {NULL, 0, NULL, 0}				/* Required at the end of array */
  };

  char* input_filename;
  char* output_filename;
  int nphi, nshift, npixels, nslices, nthreads = 1, std_output = 0;
  double radius;

  program_name = argv[0];

  /* Reading long options */
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);

    switch (next_option) {
    case 'h':
      usage_print(stdout);
      exit(EXIT_SUCCESS);
    case 'p': 
      nphi = atoi(optarg);
      assert(nphi > 0);
      printf("  Number of projections: %d\n", nphi);
      break; 
    case 's':
      nshift = atoi(optarg);
      assert(nshift > 0);
      printf("  Number of shifts per projection: %d\n", nshift);
      break;
    case 'g':
      npixels = atoi(optarg);
      assert(npixels > 0);
      printf("  Number of pixels per dimension: %d\n", npixels);
      break;
    case 'z':
      nslices = atoi(optarg);
      assert(nslices > 0);
      printf("  Number of slices in z-direction: %d\n", nslices);
      break;
    case 'r':
      radius = atof(optarg);
      assert(radius > 0);
      printf("  Radius of the support in cm: %f\n", radius);
      break;
    case 'i':
      input_filename = optarg;
      printf("  Input file with test-function data: %s\n", optarg);
      break;
    case 'o':
      output_filename = optarg;
      assert(optarg[0] != '-');
      printf("  Output file with ray data: %s\n", optarg);
      break;
    case 'n':
      nthreads = atoi(optarg);
      assert(nthreads > 0);
      printf("  Number of OpenMP threads: %d\n", nthreads);
      break;
    case 'v':
      std_output = 1;
      break;
    case ':': /* missing argument */
      fprintf(stderr, "%s: option '-%c' requires an argument\n", argv[0], optopt);
      usage_print(stderr);
      exit(EXIT_SUCCESS);
    case '?':
      usage_print(stderr);
      exit(EXIT_SUCCESS);
    case -1: /* Done with the options. */
      break;
    default:
      printf("?? getopt returned character code 0%o?? Abort.\n", next_option);
      exit(EXIT_SUCCESS);
    }
  } while (next_option != -1);

  if (optind < argc) {
    printf("Warning: non-option ARGV-arguments are omitted. See '-h' or '--help'.\n");
  }

  /**************************************** START INITIALIZATIONS *****************************************/

  // supress stdoutput if no '-v' chosen
  int std_fd;
  if (std_output == 0) {
    std_fd = supress_stdout();
  }

  /* Init values of the test function on the grid */
  printf("  Initializing test-function data...");
  double ***function_values;
  if (function_values_alloc(npixels, nslices, &function_values) != 0) {
     perror("Aborting : function_values_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (function_values_init(input_filename, npixels, nslices, function_values) != 0) {
     perror("Aborting : function_values_init.\n");
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
  if (ray_grid_init(nshift, nphi, shift, phi, radius) != 0 ) {
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

        printf("  Thread %d. Domain of work (z slices): start %d, end %d, size %d\n", thread_ID, i_slice_min, i_slice_max-1, block_size);

	//allocate memory for local thread data
	double *local_thread_data = (double*)malloc(sizeof(double)*block_size*nphi*nshift);

        //open chunk ouptut file
        FILE *foutput;
       	foutput = fopen(chunks_filenames[thread_ID], "w");

	int i_slice, i_phi, i_shift, idx = 0; // idx -- linear index to store local data

        //iterate ray integrals over z-slices
        for (i_slice = i_slice_min; i_slice < i_slice_max; ++i_slice) {
            for (i_shift = 0; i_shift < nshift; ++i_shift) {
                //iterate ray integrals over longitude angle
                for (i_phi = 0; i_phi < nphi; ++i_phi) {
		      //ray integral
		      double ray_transf = ray_transform(function_values, npixels, nslices, radius, i_slice, shift[i_shift], phi[i_phi]);
		      //store to buffer
		      local_thread_data[idx] = ray_transf;
		      idx = idx + 1;

                }
            }
        }
	fclose(foutput);

	// open thread output file, write down thread buffer in one iteration, free memory
	FILE *thread_file_output;
	thread_file_output = fopen(chunks_filenames[thread_ID], "w");
	fwrite(local_thread_data, sizeof(double), block_size * nphi * nshift, thread_file_output);
	free(local_thread_data);
	fclose(thread_file_output);
	
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
  readme_generate(input_filename, output_filename, npixels, nslices, nphi, nshift, radius);
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

  if (function_values_clean(npixels, nslices, function_values) != 0) {		//clean data for the test-function
     fprintf(stderr, "Warning : function_values_clean.\n");
  }

  if (ray_grid_clean(shift, phi) != 0) {
     fprintf(stderr, "Warning : ray_grid_clean.\n");
  }
  printf("Done.\n");

  /******************************************* END CLEANING MEMORY ****************************************/

  exit(EXIT_SUCCESS);
}
