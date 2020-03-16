#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //physical timer 
#include "init.h"
#include "integral.h"
#include "test_function.h"

#define FILENAME_MAXSIZE 128

/*
 * This program computes ray transforms of a test-function in 3D. The test-function 
 * is described by an analytic expression in 'test_function.c, test_function.h' files. 
 * The set of rays along which the computations are performed corresponds to 
 * the parallel scanner geometry : that is the support of the test-function is sliced 
 * by a set of parallel planes, which are parallel to the XY plane. Planes are parametrized by 
 * z-coordinate in [-1, 1]. 
 *
 * In each plane the set of rays is taken over which the integrals are computed corresponds to 
 * parallel-beam geometry. That is rays are parametrized by (s, phi), x*cos(phi)+y*sin(phi)=s,
 * where s is in [-1,1] and  phi in (0,2*pi). 
 *
 * The grids of s, phi are uniform in their intervals. For integration along a line a 
 * Simpson's rule is used.
 * 
 * IMPORTANT : It is always assumed that the test-function is compactly supported with support in the
 * centered unit ball.
 * 
 * The parameters for the grids of rays (nslices, nshift, nphi)  are taken from configuration file (-p filename). 
 *
 * The result of computations is stored in a separate csv file (-o filename) in the following order: 
 *
 *   for islice = (0 : nslice-1) // z-slice
 *      for ishift = (0 : nshift-1) // value for s for rays in the plane
 *         for iphi = (0 : nphi-1) // value for phi for rays in the plane
 *
 *             fprintf(output, ray_transform(f, "%lf\n", islice, ishift, iphi)
 *         
 *         endfor
 *      endfor
 *   endfor
 *   
 *
*/



/* ***************************************************************************
 * *****************  HANDLING FILES (README, AGGREGATION)********************
 * ***************************************************************************
 */

void readme_generate(char* output_filename, int nslices, int nshift, int nphi) {

  FILE *doc_file;
  char doc_filename [ FILENAME_MAXSIZE ];
  
  sprintf(doc_filename, "readme_%s", output_filename);

  doc_file = fopen(doc_filename, "w");
  fprintf(doc_file, "  Ray transforms in 3D, output file: %s\n", output_filename);
  fprintf(doc_file, "  Format of the data: [z], [s], [phi], [Pf]\n");
  fprintf(doc_file, "  Parameters of the data:\n");
  fprintf(doc_file, "     number of z-slices [nz]   : %d,\n", nslices);
  fprintf(doc_file, "     number of shifts   [s]    : %d.\n", nshift);
  fprintf(doc_file, "     number of angles   [phi]  : %d,\n", nphi);
    
  fclose(doc_file);
}

void chunks_aggregate(char* output_filename, char ** chunks_filenames, int nchunks) {
  
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


const char* program_name;


void print_usage(FILE* stream, int exit_code) {
  fprintf(stream, "Usage: %s -p (filename) -o (filename) -n (number)\n", program_name);
  fprintf(stream, 
	  "   -h --help                 Display the usage information.\n"
	  "   -p --parameters filename  File with parameters of the ray grid from file.\n"
	  "   -o --output filename      Filename for output.\n"
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
  const char* short_options = "hp:o:n:";
  const struct option long_options[] = {
    {"help", no_argument, 0, 0}, 
    {"parameters", required_argument, 0, 0},
    {"output", required_argument, 0, 0},
    {"nthreads", required_argument, 0, 0},
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
  
  int nshift = parameters[0];
  int nphi   = parameters[1];
  int ngrid  = parameters[2];		// number of z-slices; (affects density of points for and integral)
  
  /* Init grid in Radon space */
  double *phi, *shift;
  printf("  Initializing the grid in ray transforms space...");
  if (ray_grid_alloc(nshift, nphi, &shift, &phi) != 0) {		// allocate memory for arrays : phi, shift
     perror("Aborting : ray_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (ray_grid_init(nshift, nphi, shift, phi) != 0) {			// set values of : phi, shift
     perror("Aborting : ray_grid_init.\n");
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
  chunks_filenames_alloc(output_filename, nthreads, &chunks_filenames, FILENAME_MAXSIZE);
  chunks_filenames_init(output_filename, nthreads, chunks_filenames);
  
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
		    double rt = ray_integral(&test_function, 	
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
  
  readme_generate(output_filename, ngrid, nshift, nphi);
  printf("  Aggregating chunks to one output file...");
  chunks_aggregate(output_filename, chunks_filenames, nthreads);
  printf("Done.\n");
  
  /*********************************** END AGGREGATING DATA TO ONE FILE ***********************************/
  

  /******************************************* START CLEANING MEMORY **************************************/
  
  printf("  Cleaning allocated memory...");
  
  if (chunks_files_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_files.\n");
  }
  if (chunks_filenames_clean(nthreads, chunks_filenames) != 0) {
     fprintf(stderr, "Warning : clean_chunks_filenames.\n");
  }
 
  if (parameters_clean(parameters) != 0) {
     fprintf(stderr, "Warning : clean_parameters.\n");
  }
  if (ray_grid_clean(shift, phi) != 0) {
     fprintf(stderr, "Warning : clean_ray_grid.\n");
  }
  printf("Done.\n");
  
  /******************************************* END CLEANING MEMORY ****************************************/
  
  exit(EXIT_SUCCESS);
}
