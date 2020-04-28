#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>   // input options handling
#include <assert.h>
#include <omp.h>      // parallelization in the set of 2d-planes
#include <sys/time.h> // measure of parallel time evaluation

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "init.h"
#include "slicing.h"


#define FILENAME_MAXSIZE 128

/*
  * Program performes reduction of data given by (weighted) ray transforms in 3D (slices by slices)
  * to the data given by (weighted) Radon transforms along two-dimensional planes in 3D.
  *
  * Such reduction is a consequence of natural fibration of plane into parallel lines (even for
  * weighted ray-Radon transforms, see Goncharov, Novikov, An analog of Chang inversion formula in 
  * multidimensions, 2016).
  * 
  * IMPORTANT : It is always assumed that the test-function is compactly supported in the
  * centered unit ball. Otherwise the generated data will be incomplete and this will produce artifacts in
  * reconstructions.
  *
  * Output is stored in CSV format in one column in the following order:
  *
  * 	for shift = shift(0 : nshift-1)
  *       for phi = phi(0 : nphi-1)
  *         for ntheta = theta(0 : ntheta-1)
  *	           Rf(shift, phi, theta)
  *	        endfor
  *	      endfor
  *	  endfor
  *
  * SOME TECHNICAL REMARKS:
  * 	shift : distance from the origin to the plane (shifts vary in [-1,1] including endpoints)
  * 	phi   : longitude angle (0, 2*pi) (phi vary in (0, 2*pi) so that the circle is partioned in equisized intervals)
  * 	theta : latitude angle (0, pi) (correspond to Gauss-Lebato grid, i.e., theta_j = arccos(tj), tj - Gauss points on [-1,1])
  *     nslices : number of z-slices for the ray transforms, also assumes to be the number of pixels in the underline image 
  */

int stdout_supress() {

    int fd_stdout = dup(1);
    int nullfd = open("/dev/null", O_WRONLY);

    dup2(nullfd, 1);
    close(nullfd);

    return fd_stdout;
}


int main(int argc, char ** argv) {

  //------------------------------------ARGUMENTS HANDLING------------------------------------------------------------//

if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help'.\n", argv[0]);
     exit(EXIT_FAILURE);
  }

  char *program_name = argv[0];

  int next_option;
  
  const char* short_options = "hs:p:z:t:i:o:n:v";
  const struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"nshift", required_argument, NULL, 's'},
    {"nphi", required_argument, NULL, 'p'},
    {"nslice", required_argument, NULL, 'z'},
    {"ntheta", required_argument, NULL, 't'},
    {"input_file", required_argument, NULL, 'i'},
    {"output_file", required_argument, NULL, 'o'},
    {"nthreads", required_argument, NULL, 'n'},
    {"verbose", no_argument, NULL, 'v'},
    {NULL, 0, NULL, 0}				/* Required at the end of array */
  };
  


  char* input_filename;
  char* output_filename;
  int   nthreads = 1, std_output = 0;

  // default values
  int nshift;
  int nphi;
  int ntheta;
  int nslice;

  //------------------------------ Reading options-------------------------------//
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);

    switch (next_option) {
      case 'h':
	     print_usage(stdout, program_name);
  	     exit(EXIT_SUCCESS);
      case 's':
	     nshift = atoi(optarg);
	     printf("Number of shifts : %d\n", nshift);
	     break;
      case 'p':
    	     nphi = atoi(optarg);
	     printf("Number of phi / projections per slice : %d\n", nphi);
	     break;
      case 'z':
	     nslice = atoi(optarg);
	     printf("Number of z-slices : %d\n", nslice);
	     break;
      case 't':
	     ntheta = atoi(optarg);
	     printf("Number of theta : %d\n", ntheta);
	     break; 
      case 'i':
	     input_filename = optarg;
	     printf("Input file with ray data: %s\n", input_filename);
	     break;
      case 'o':
	     output_filename = optarg;
	     printf("Output file with Radon data: %s\n", output_filename);
	     break;
      case 'n':
	     nthreads = atoi(optarg);
	     assert(nthreads > 0);
	     printf("Number of threads: %d\n", nthreads);
	     break;
      case 'v':
	     std_output = 1;
	     break;
      case ':': /* missing argument */
	     fprintf(stderr, "%s: option '-%c' requires an argument\n",
       argv[0], optopt);
	     exit(EXIT_SUCCESS);
      case '?': /* help */
	     print_usage(stderr, program_name);
      case -1: /* Done with the options. Cycle will quit */
	     break;

      default:
	     printf("?? getopt returned character code 0%o??\n", next_option);
	     exit(EXIT_SUCCESS);
    }
  } while (next_option != -1);

  if (optind < argc) {
    printf("Warning: non-option ARGV-arguments are omitted. Try '-h' or '--help'.\n");
  }

  //-------------------------------------- INITIALIZATIONS ------------------------------------------------------------------------------/

  int fd_stdout;
  if (std_output == 0) {
    fd_stdout = stdout_supress();
  }

  // --------------------------------------------Read ray transform data from file --------------------------------------------------------/


  printf("  Initializing data for ray transforms...");
  double*** ray_values;
  if (ray_values_alloc(nshift, nphi, nslice, &ray_values) != 0) {
     perror("Aborting : ray_values_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (ray_values_init(input_filename, nshift, nphi, nslice, ray_values) != 0) {
     perror("Aborting : ray_values_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");

  //---------------------------------------------------------------------------------------------------------------------------------------/

  //------------------------------------------------------ init slice in Radon space -------------------------------------------------------/
  double *shift, *phi, *theta;
  printf("  Initializing the slice in Radon space...");
  if (radon_grid_alloc(nshift, nphi, ntheta, &shift, &phi, &theta) != 0) {		// allocate memory for arrays shift, theta, phi
     perror("Aborting : radon_slice_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (radon_grid_init(nshift, nphi, ntheta, shift, phi, theta) != 0) {		        // set values of shift, theta, phi
     perror("Aborting : radon_grid_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");

  //------------------------------------ END OF INITIALIZATIONS ------------------------------------------------------------------------/

  //------------------------ START REDUCTION RAY TRANSFORMS -> RADON TRANSFORMS (PARALLELIZATION VIA OPENMP) ---------------------------/
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

    //set equally shifts of planes for each thread
    int block_size = nshift / nthreads + (thread_ID < (nshift % nthreads) ? 1 : 0);
    int shift_min = (thread_ID < (nshift % nthreads) ? (nshift / nthreads + 1)*thread_ID : (nshift / nthreads)*thread_ID + (nshift % nthreads)),
	shift_max = shift_min + block_size;

    printf("  Thread %d. Domain of work (shifts): start %d, end %d, size %d\n", thread_ID, shift_min, shift_max-1, block_size);


    // allocate memory for local storage (I/O operations on hard drive will be done in one iteration in the end)
    double *local_output_buffer;
    local_output_buffer = (double*)malloc(block_size * nphi * ntheta * sizeof(double));

    int i_shift, i_phi, i_theta, idx = 0;
    // ray integrals in the given OZ slice
    for (i_shift = shift_min; i_shift < shift_max; ++i_shift) {

	    // iterate ray integrals over shifts in slice plane
	    for (i_phi = 0; i_phi < nphi; ++i_phi) {

		    // iterate ray integrals over polar angle
		    for (i_theta = 0; i_theta < ntheta; ++i_theta) {

		      // reduction of ray integrals to Radon transform over the plane (shift[i_shift], phi[i_phi], theta[i_theta])
		      double rt = reduction(ray_values,
		                     shift, phi, theta,
                                     nshift, nphi, ntheta, nslice,
                                     i_shift, i_phi, i_theta);
		      local_output_buffer[idx] = rt; // write computation to buffer
		      idx = idx + 1;
                    }
            }
    }

    // open output chunk file
    FILE *file_output;
    file_output = fopen(chunks_filenames[thread_ID], "w");
    fwrite(local_output_buffer, sizeof(double), block_size * nphi * ntheta, file_output);
    fclose(file_output);
    free(local_output_buffer);
    printf("  Thread %d. Job is done.\n", thread_ID);
    
  } // OpenMP section finish

  // stop timer
  gettimeofday(&end, NULL);
  double delta_t = ((end.tv_sec  - start.tv_sec) * 1000000u +
         end.tv_usec - start.tv_usec) / 1.e6;
  printf("\n");
  printf("  Elapsed time: %3.1f sec\n", delta_t);

  //------------------------ END REDUCTION RAY TRANSFORMS -> RADON TRANSFORMS --------------------------------------------------------/


  //------------------------------------------ AGGREGATING DATA TO ONE FILE ----------------------------------------------------/

  printf("  Making a readme file...\n");
  generate_readme(input_filename, output_filename, nphi, ntheta, nshift);
  printf("  Done.\n");
  printf("  Aggregating chunks to one output file...");
  chunks_aggregate(output_filename, chunks_filenames, nthreads);
  printf("  Done.\n");

  //---------------------------------------------- END AGGREGATING DATA TO ONE FILE --------------------------------------------------/


  //------------------------------------------------------- CLEANING MEMORY ----------------------------------------------------/

  printf("  Cleaning allocated memory...");
  if (chunks_files_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : chunks_files_clean.\n");
  }
  if (chunks_filenames_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : chunks_filenames_clean.\n");
  }
  if (ray_values_clean(nshift, nphi, nslice, ray_values) != 0) {
     fprintf(stderr, "Warning : ray_values_clean.\n");
  }
  if (radon_grid_clean(shift, phi, theta) != 0) {
     fprintf(stderr, "Warning : radon_slice_clean.\n");
  }
  printf("Done.\n");

  //----------------------------------------------------------- END CLEANING MEMORY --------------------------------------------------/

  exit(EXIT_SUCCESS);
}
