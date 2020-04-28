#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //for measurement of physical evaluation time
#include "init.h"
#include "integral.h"

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>

#define FILENAME_MAXSIZE 128
/*
 * This program computes Radon transforms along 2d-planes of a test function in 3D. The test function 
 * is described by it values on the uniform grid in [-1,1]^3 which are stored in a separate binary file. 
 * The set of planes along which computations are performed has a specific distribution, which is 
 * uniform in some variables and non-uniform in others. 
 *
 * Planes have three parameters (shift, phi, theta). Angles 'phi', 'theta'
 * vary in intervals [0, 2pi], [0, pi], respectively, shifts vary in [-1,1].
 *
 * 'phi' is uniformly distributed in [0, 2pi]
 * 'theta' corresponds to Gauss-Lebato grid in [0, pi], that is theta_j = arccos(tau_j), where tau_j are 
 *  Gauss-Legendre-Lebato quadrature points in [-1,1]. 
 * 
 *  For integration along planes the zero-th order (Riemann) rule is used.  
 * 
 * 
 * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
 * centered ball of radius given in launch parameters.
 *
 * The result is stored in a separate binary file (-o filename) in the following order: 
 *
 *    for (ishift = 0 : nshift-1)
 *      for (iphi = 0 : nphi-1) 
 *        for (itheta = 0 : ntheta-1) 
 *           
 *           fwrite(output, Radon-transform(shift, phi, theta));
 *
 *        endfor
 *      endfor
 *    endfor
*/


/* 
 * ****************************** HANDLING FILES *****************************
 */
void readme_generate(char *input_filename, char *output_filename, int npixels, int nphi, int ntheta, int nshift, double radius) {
  FILE *readme_file;
  char readme_filename [ FILENAME_MAXSIZE ];
  
  sprintf(readme_filename, "readme_%s", output_filename);

  readme_file = fopen(readme_filename, "w");
  fprintf(readme_file, "  Radon transforms in 3D for test-function from  %s stored to file %s.\n", input_filename, output_filename);
  fprintf(readme_file, "  Parameters of the data:\n");
  fprintf(readme_file, "    number of pixels per dimension : %d\n", npixels);
  fprintf(readme_file, "    number of angles of phi        : %d\n", nphi);
  fprintf(readme_file, "    number of angles of theta      : %d\n", ntheta);
  fprintf(readme_file, "    number of shifts of            : %d\n", nshift);
  fprintf(readme_file, "    radius of the support in cm    : %f\n", radius);
    
  fclose(readme_file);
}

int stdout_supress() {
  int ret = dup(1);
  int nullfd = open("/dev/null", O_WRONLY);
  
  dup2(nullfd, 1);
  close(nullfd);

  return ret;
}

const char* program_name;

/* Prints usage information to STREAM (stdout or stderr) */

void print_usage(FILE* stream) {
  fprintf(stream, "Usage: %s -i (filename) -o (filename) -n (nthreads) [ARGS]\n", program_name);
  fprintf(stream, 
	  "   -h --help      no_arg     Display this usage information.\n"
	  "   -p --nphi      integer    Number of phi angles.\n"
	  "   -t --ntheta    integer    Number of theta angles.\n"
	  "   -s --nshift    integer    Number of shifts per direction \n"
	  "   -r --radius    float      Radius of the test-function support in cm.\n"
	  "   -g --npixels   integer    Number of pixels per dimension for the test-function.\n"
	  "   -i --input     string     Path to input file with test-function.\n"
	  "   -o --output    string     Path to output file.\n"
	  "   -n --nthreads  integer    Number of OpenMP threads for parallelization.\n"
	  "   -v --verbose   no_arg     Print extra information during computations.\n");
}

/* 
 * ****************************** ENTRY POINT ********************************
 */

int main(int argc, char * argv[]) {
  
  if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help'.\n", argv[0]);
     exit(EXIT_FAILURE);
  }
  
  int next_option;
  const char* short_options = "hp:t:s:g:r:i:o:n:";
  const struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"nphi", required_argument, NULL, 'p'}, //nphi
    {"ntheta", required_argument, NULL, 't'}, //ntheta
    {"nshift", required_argument, NULL, 's'}, //nshift
    {"npixels", required_argument, NULL, 'g'}, //npixels 
    {"radius", required_argument, NULL, 'r'}, //radius
    {"input", required_argument, NULL, 'i'}, 
    {"output", required_argument, NULL, 'o'},
    {"nthreads", required_argument, NULL, 'n'},
    {"verbose", no_argument, NULL, 'v'}, 
    {NULL, 0, NULL, 0}				/* Required at the end of array */
  };

  char *input_filename; 
  char *output_filename;
  int nphi, ntheta, nshift, npixels, nthreads=1, std_output = 0;
  double radius;
  
  program_name = argv[0];

  /* Reading long options */ 
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);
    
    switch (next_option) {
      case 'h':
	print_usage(stdout);
	exit(EXIT_SUCCESS);
      case 'p':
        nphi = atoi(optarg);
	assert(nphi > 0);
	printf("  Number of nphi : %d\n", nphi);
	break;
      case 't':
	ntheta = atoi(optarg);
	assert(ntheta > 0);
	printf("  Number of ntheta : %d\n", ntheta);
	break;
      case 's':
	nshift = atoi(optarg);
	assert(nshift > 0);
	printf("  Number of nshift : %d\n", nshift);
        break;
    case 'g':
        npixels = atoi(optarg);
	assert(npixels > 0);
	printf("  Number of pixels per dimension : %d\n", nshift);
        break;
      case 'r':
        radius = atof(optarg);
	assert(radius > 0);
	printf("  Radius of the support in cm : %f\n", radius);
        break;
      case 'i':
	input_filename = optarg;
	printf("Input file with test-function data: %s\n", optarg);
	break;
      case 'o':
	output_filename = optarg;
	printf("Output file with Radon data: %s\n", optarg);
	break;
      case 'n':
	nthreads = atoi(optarg);
	assert(nthreads > 0);
	printf("Number of threads: %d\n", nthreads);
	break;
      case 'v':
        std_output = 1;
      case ':': /* missing argument */ 
	fprintf(stderr, "%s: option '-%c' requires an argument\n",
		argv[0], optopt);
	print_usage(stderr);
      case '?':
	print_usage(stderr);
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

  int fd_stdout;
  if (std_output == 0) {
    fd_stdout = stdout_supress();
  }
  
  /* Init values of the test function on the grid */
  printf("  Initializing data for the test-function...");
  double ***function_values;
  if (function_values_alloc(npixels, &function_values) != 0) {
     perror("Aborting : function_values_alloc\n");
     exit(EXIT_FAILURE);
  }
  if (function_values_init(input_filename, npixels, function_values) != 0) {
     perror("Aborting : function_values_init\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  /* Init grid in Radon space */
  double *phi, *theta, *shift;
  printf("  Initializing the grid in Radon space...");
  if (radon_grid_alloc(nphi, ntheta, nshift, &phi, &theta, &shift) != 0) {
     perror("Aborting : radon_grid_alloc\n");
     exit(EXIT_FAILURE);
  }
  if (radon_grid_init(nphi, ntheta, nshift, phi, theta, shift) != 0) {
     perror("Aborting : radon_grid_init\n");
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
        // set equally parts for each thread 
        int thread_ID = omp_get_thread_num();
        int block_size = nshift / nthreads + (thread_ID < (nshift % nthreads) ? 1 : 0);	
        int i_shift_min = (thread_ID < (nshift % nthreads) ? (nshift / nthreads + 1)*thread_ID : (nshift / nthreads)*thread_ID + (nshift % nthreads)),
	    i_shift_max = i_shift_min + block_size;
	    
        printf("  Thread %d. Domain of work (shifts): start %d, end %d, size %d\n", thread_ID, i_shift_min, i_shift_max-1, block_size);

        //allocate memory for temporary data
	double *local_thread_data = (double*)malloc(sizeof(double) * block_size * nphi * ntheta);
	
        int i_phi, i_theta, i_shift, idx = 0;    
        //iterate plane integrals over shifts along normals
        for (i_shift = i_shift_min; i_shift < i_shift_max; ++i_shift) {      
      
            //iterate plane integrals over longitude angle
            for (i_phi = 0; i_phi < nphi; ++i_phi) {
	
	              //iterate plane integrals over latitude angle
	              for (i_theta = 0; i_theta < ntheta; ++i_theta) {
		               //plane integral
		               double radon_transform = plane_integral(function_values, npixels, 
								       phi[i_phi], theta[i_theta], shift[i_shift],
								       radius);
			             local_thread_data[idx] = radon_transform;
		               idx = idx + 1; 
                }
            }
        }

	//open ouptut file, write down data in one iteration, free memory
        
	FILE *thread_file_output;

        thread_file_output = fopen(chunks_filenames[thread_ID], "w");
	fwrite(local_thread_data, sizeof(double), block_size * nphi * ntheta, thread_file_output);
	free(local_thread_data);
	fclose(thread_file_output);
	
	printf("  Thread %d. Job is done.\n", thread_ID);
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
  readme_generate(input_filename, output_filename, npixels, nphi, ntheta, nshift, radius);
  printf("  Done.\n");
  printf("  Aggregating chunks to one output file...");
  chunks_aggregate(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");
  
  /*********************************** END AGGREGATING DATA TO ONE FILE ***********************************/
  

  /******************************************* START CLEANING MEMORY **************************************/
  
  printf("  Cleaning allocated memory...");
  
  if (chunks_files_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : chunk_files_clean.\n");
  }
  if (chunks_filenames_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : init_values.\n");
  }
 
  if (function_values_clean(npixels, function_values) != 0) {
     fprintf(stderr, "Warning : function_clean_values.\n");
  }
  if (radon_grid_clean(phi, theta, shift) != 0) {
     fprintf(stderr, "Warning : radon_grid_clean.\n");
  }
  printf("Done.\n");
  
  /******************************************* END CLEANING MEMORY ****************************************/
  
  exit(EXIT_SUCCESS);
}
