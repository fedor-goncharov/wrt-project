#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <omp.h>
#include <sys/time.h> //measure physical time of evaluation
#include <stdbool.h>  //boolean lib for flags

#include "init.h"
#include "integral.h"
#include "aggregate.h"

#define FILENAME_MAXSIZE 128


/* Circular harmonic expansions of spect weight in 3D in 2D-slices. That is in each 2D slice weight attains circlular harmonic
 * expansion and its results are stored in a separate file in binary format.
 * 
 * The program reads the data for an attenuation map from file 
 * (-a [attenuation map] filename) and computes the exponent expansion on the uniform cubic grid of the spect weight 
 * The set of angles  in each slice plane corresponds to 'uniform geometry'. 
 * SPECT weight is computed for each point of the grid (for better memory usage (and worse computation time)).
 * The parameters of the input are explained in the usage info ('-h' , '--help' for usage).
 * 
 * IMPORTANT : It is always assumed that the test-function is compactly supported inside the centered ball of radius [radius]. 
 * Otherwise the generated data will be incomplete and this will produce artifacts in reconstructions.
 *
 * INPUT : 
 *   npixels -- 
 *   nprojections -- 
 *   nradius -- 
 *   attenuation_map -- path to file with attenuation map, values should be given with separator "\n"
 *   max_degree -- maximal degree of expansions (should be non-negative)
 *   divergence_form -- computation of weight expansions in divergence form
 *   normalization -- multiply expansions with (1/2*PI) factor
 *   nthreads -- number of OpenMP threads for execution 
 *   
 * OUTPUT : 
 * The results of computations are stored in separate binary files (-o filename). 
 * 
 * NOTE : This code is not optimal at all 
 * -- many IO operations in parallel (slow)
 * -- computations are designed for array of frequenceis (code should be splitted to compute frequencies one by one) 
 * -- parallelization should be performed only for one file per evaluation
 * -- OpenMP is not scalable at all (fine withing one machine, not fine with cluster) -- better to use MPI to run on big cluster machines
 * -- code will be very reusable for Kunyansky algorithm -- should be compiled in a shared library .so
 * -- config file runs are very uncomfortable (commands are needed to run often -- changing config files is boring)
*/

/*
 * ****************************** OUTPUT README ******************************
 * 
 */

void generate_readme(char* output_filename, char* input_att_filename,
		     int max_degree, int npixels, double radius, int nphi,
		     bool div_form, bool normalization) {
  FILE *file;
  char readme_filename [FILENAME_MAXSIZE];
  
  sprintf(readme_filename, "readme_%s", output_filename);

  file = fopen(readme_filename, "w");
  fprintf(file, "  Expansions of SPECT weight (attenuation map was in %s) in circular harmonics in 2D slices.\n", input_att_filename);
  fprintf(file, "  Format of the data: [iteration over z, x, y of W_m(x,y,z)]\n");
  fprintf(file, "  Parameters of the data:\n\n");
  fprintf(file, "  \t Max degree of expansions : %d,\n", max_degree);
  fprintf(file, "  \t Number of pixels on cubic grid : %d,\n", npixels);
  fprintf(file, "  \t Radius of the support in cm : %3.1f,\n", radius);
  fprintf(file, "  \t Number of projections per slice was used : %d\n", nphi);
  fprintf(file, "  \t Divergence form of output : %s\n", div_form ? "true" : "false");
  fprintf(file, "  \t Normailization factor : %s\n", normalization ? "true" : "false");
  fclose(file);
}


/* 
 * ****************************** USAGE **************************************
 * 
 */


/* Prints usage information to STREAM (stdout or stderr) */

void print_usage(FILE* out_stream, const char* program_name) {
  fprintf(out_stream, "Usage: %s -p [file] -a [file] -m [integer] [-d] -o [file] -n [integer]\n", program_name);
  fprintf(out_stream,  
	  "   -h --help             no_argument Display usage information.\n"
	  "   -a --attenuation_map  filename    Read attenuation map from file (according to configuration parameters in '-p').\n"
	  "   -s --npixels          integer     Number of pixels per dimension in the attenuation map file.\n"
	  "   -p --nprojections     integer     Number of projections per slice.\n"
	  "   -r --radius           float       Radius of the support in cm.\n"
	  "   -m --max_degree       integer     Maximal degree of exponent expansions.\n"
	  "   -f --normalization    no_argument Normalization by multplicative factor 1/(2*M_PI).\n"
	  "   -d --divergence_form  no_argument Weight is computed in the divergence form, that is factor M_PI/2 is used in computations.\n"
	  "   -o --output           filename    Place output data to a file.\n"
	  "   -n --nthreads         integer     Number of OpenMP threads for parallelization of calculations.\n");
}



/* ***************************************************************************
 * ****************************** ENTRY POINT ********************************
 * ***************************************************************************
 */

int main(int argc, char* argv[]) {
  
  int next_option;
  const char* short_options = "ha:s:p:r:m:o:n:";
  const struct option long_options[] = {
    {"help", no_argument, NULL, 'h'}, 
    {"attenuation", required_argument, NULL, 'a'},
    {"npixels", required_argument, NULL, 's'},
    {"nprojections", required_argument, NULL, 'p'},
    {"radius", required_argument, NULL, 'r'},
    {"max_degree", required_argument, NULL, 'm'},
    {"divergence_form", no_argument, NULL, 'd'},
    {"normailization", no_argument, NULL, 'f'},
    {"output", required_argument, NULL, 'o'},
    {"nthreads", required_argument, NULL, 'n'},
    {NULL, 0, NULL, 0}
  };

  char *input_att_filename,
       *output_filename_template;
  
  int  nthreads     = 1,   // number of openmp threads (by default 1)
       max_degree   = 0,   // max degree of expansions [-max_degree, .. , max_degree]
       ndegrees     = 1,   // number of expansions = 2*max_degree + 1 (because real and imaginary part)
       npixels      = 129, // number of pixels per dimension
       nphi         = 129; // number of projections per slice
  double radius = 1.0;     // radius of the support
  bool divergence_form = false, // divergence form flag
       normalization = false; //  normalization factor flag
       
  
  const char *program_name = argv[0];

  /* reading long input options */ 
  if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help'.\n", argv[0]);
     exit(EXIT_FAILURE);
  }

  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);
    
    switch (next_option) {
      case 'h':
	      print_usage(stdout, program_name);
              exit(EXIT_SUCCESS);
      case 's':
	      npixels = atoi(optarg);
	      printf("  Number of pixels per dimension in the attenuation map : %d\n", npixels);
	      break;
      case 'p':
	      nphi = atoi(optarg); 
	      printf("  Number of projections : %d\n", nphi);
	      break;
      case 'r':
	      radius = atof(optarg);
	      printf("  Radius of the support in cm : %f\n", radius);
	      break;
      case 'a':
	      input_att_filename = optarg;
	      printf("  Input file with attenuation map data: %s\n", optarg);
	      break;
      case 'm':
	      max_degree = atoi(optarg);
              ndegrees = max_degree + 1; // frequencies for computation {0...max_degree}
	      assert(max_degree >= 0);
	      printf("  Maximal expansion degree: %d\n", max_degree);
	      break;
      case 'd':
	      divergence_form = true; 
	      if (divergence_form) {
		printf("  Divergence form is used for expansions.\n");
	      } else {
		printf("  Divergence form is off.\n");
	      }
	      break;
      case 'f':
       	      normalization = true;
	      if (normalization) {
	        printf("  Normalization factor is on.\n");
	      } else {
		printf("  Normalization factor is off.\n");
	      }
	      break;
      case 'o':
	      output_filename_template = optarg;
	      printf("  Output file with ray data: %s\n", optarg);
	      break;
      case 'n':
	      nthreads = atoi(optarg);
	      assert(nthreads > 0);
	      printf("  Number of OpenMP threads: %d\n", nthreads);
	      break;
      case ':': /* missing argument */ 
	      fprintf(stderr, "Call %s: option '-%c' requires an argument. Try '-h' or '--help' for usage info.\n", argv[0], optopt);
              exit(EXIT_SUCCESS);
      case '?':
	      print_usage(stderr, program_name);
              exit(EXIT_SUCCESS);
      case -1:     /* Done with the options. */
	      break;
      default:
	      abort();
    }
  } while (next_option != -1);
  
  if (optind < argc) {
    printf("  Warning: non-option ARGV-arguments are omitted. See '-h' or '--help'.\n"); 
  }
  
  /**************************************** START INITIALIZATIONS *****************************************/
  
  /* Init values for the attenuation map on the grid : allocate memory, read from file */
  printf("  Initializing attenuation map data...");
  double*** att_values;
  if (att_values_alloc(&att_values, npixels) != 0) {
     perror("Aborting : att_values_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (att_values_init(att_values, input_att_filename, npixels) != 0) {	// read attenuation map from file
     perror("Aborting : att_values_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n");
  
  
  /* Init grid for directions for rays : allocate memory and set angles */
  double *phi;
  printf("  Initializing the grid in ray transforms space...");
  if (ray_grid_alloc(&phi, nphi) != 0) {
     perror("Aborting : ray_grid_alloc.\n");
     exit(EXIT_FAILURE);
  }
  if (ray_grid_init(phi, nphi) != 0 ) {
     perror("Aborting : ray_grid_init.\n");
     exit(EXIT_FAILURE);
  }
  printf("Done.\n\n");
  
  /************************************** END OF INITIALIZATIONS ******************************************/
  
  
  /************************* START RADON TRANSFORM WITH PARALLELIZATION VIA OPENMP ************************/
  
  // Init output filenames 
  char*** output_filenames;  // ndegrees x 2 (real/imag)
  output_filenames_alloc(&output_filenames, ndegrees, FILENAME_MAXSIZE);
  output_filenames_init(output_filenames, output_filename_template, ndegrees);
  
  // Chunks filenames 
  char**** chunks_filenames; // nthreads x ndegrees x 2 (real/imag) i.e. each file is splitted in nthread chunks 
  chunks_filenames_alloc(&chunks_filenames, nthreads, ndegrees, FILENAME_MAXSIZE);
  chunks_filenames_init(chunks_filenames, output_filenames, nthreads, ndegrees);

  // precomputations for threds
  const double delta = 2.0*radius/npixels;     // step size in the grid
  const double half_delta = delta / 2.0;  // half step size in the grid (half of the pixel size)
  
  //start timer
  struct timeval start, end;
  gettimeofday(&start, NULL);
  
  // ----------------------------------------------------------- parallelization part -------------------------------------------------------
  omp_set_num_threads(nthreads);
  #pragma omp parallel
  {

          int thread_id = omp_get_thread_num();
          int block_size = npixels / nthreads + (thread_id < (npixels % nthreads) ? 1 : 0);
          int i_slice_min = (thread_id < (npixels % nthreads) ? (npixels / nthreads + 1) * thread_id : (npixels / nthreads)*thread_id + (npixels % nthreads)); 
          printf("  Thread %d. Domain of work (z slices): start %d, end %d, size %d\n", thread_id, i_slice_min, i_slice_min + block_size-1, block_size);
     

	  // allocate memory for local outputs
	  int i_degree = 0;
          double *local_output_real[ndegrees], *local_output_imag[ndegrees];
	  for (i_degree = 0; i_degree < ndegrees; ++i_degree) {
	    local_output_real[i_degree] = (double*)malloc(block_size * npixels * npixels * sizeof(double));
	    local_output_imag[i_degree] = (double*)malloc(block_size * npixels * npixels * sizeof(double));
	  }

          double x, y, z;                         // local points of the grid
          int i_z, i_x, i_y, i_phi;	          // local idx's
	  int idx = 0;                            // index of elements to write in file one-by-one
	  
          for (i_z = 0; i_z < block_size; ++i_z) {
               z = -radius  + half_delta + (i_slice_min + i_z) * delta;  // z coordinate=slice coordinate
	    
               for (i_x = 0; i_x < npixels; ++i_x) {
	            x = -radius + half_delta + i_x * delta;                 // x coordinate 
	        
                    for (i_y = 0; i_y < npixels; ++i_y) {
	                 y = -radius + half_delta + i_y * delta;               // y coordinate 
		      
                         double wspect_at_point[nphi];  // compute spect weight at (x,y,z) for all directions phi
	                 for (i_phi = 0; i_phi < nphi; ++i_phi) {  
			   wspect_at_point[i_phi] = exp_ray_transform(att_values, x, y, z, phi[i_phi], npixels, radius, divergence_form);
	                 }

			 for (i_degree = 0; i_degree < ndegrees; ++i_degree) {
		              double wcoeff_real = wspect_decomposition_real(wspect_at_point, i_degree, phi, nphi, normalization); // cosine 
		              double wcoeff_imag = wspect_decomposition_imag(wspect_at_point, i_degree, phi, nphi, normalization); // sine

			      local_output_real[i_degree][idx] = wcoeff_real;
			      local_output_imag[i_degree][idx] = wcoeff_imag;
			      idx = idx + 1;
			 }		 
                    }
               }
          }
	  
	  FILE *output_files_real[ndegrees], *output_files_imag[ndegrees];
	  for (i_degree = 0; i_degree < ndegrees; ++i_degree) {
               output_files_real[i_degree] = fopen(chunks_filenames[thread_id][i_degree][0], "wb"); // real part 
               output_files_imag[i_degree] = fopen(chunks_filenames[thread_id][i_degree][1], "wb"); // imaginary part

	       //write output to chunks
	       fwrite((local_output_real[i_degree]), sizeof(double), block_size * npixels * npixels, output_files_real[i_degree]);
	       fwrite((local_output_imag[i_degree]), sizeof(double), block_size * npixels * npixels, output_files_imag[i_degree]);

	       fclose(output_files_real[i_degree]);
	       fclose(output_files_imag[i_degree]);

	       // deallocate memory of local outputs
	       free(local_output_real[i_degree]);
	       free(local_output_imag[i_degree]);
	       
          }
          printf("  Thread %d. Job is done.\n", thread_id);
  }
  
  //stop timer
  gettimeofday(&end, NULL);
  double delta_t = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;
  printf("\n\n");
  printf("  Elapsed time: %3.1f sec\n", nthreads, delta_t);
  
  /************************* END WEIGHT DECOMPOSITIONS TRANSFORM WITH PARALLELIZATION VIA OPENMP **************************/
  
  printf("  Aggregating chunks to one output file...\n");
  aggregate_chunks(output_filenames, chunks_filenames, nthreads, ndegrees);
  printf("  Cleaning temporary files...\n");
  chunks_files_clean(chunks_filenames, nthreads, ndegrees);

  printf("  Generating readme...");
  generate_readme(output_filename_template, input_att_filename, max_degree, npixels, radius, nphi,
		  divergence_form, normalization);
  printf("  Done.\n");
  
  /******************************************* START CLEANING MEMORY **************************************/
  
  printf("  Cleaning allocated memory...");
  
  // clean memory 
  if (chunks_filenames_clean(chunks_filenames, nthreads, ndegrees) != 0) {
     fprintf(stderr, "Warning : chunks_filenames_clean.\n");
  }
  if (att_values_clean(att_values, npixels) != 0) {		// clean data for a test-function
     fprintf(stderr, "Warning : att_values_clean.\n");
  }
  if (ray_grid_clean(phi) != 0) {				// clean angles 'phi'
     fprintf(stderr, "Warning : ray_grid_clean.\n");
  }
  printf("Done.\n");
  
  /******************************************* END CLEANING MEMORY ****************************************/
  
  exit(EXIT_SUCCESS);
}
