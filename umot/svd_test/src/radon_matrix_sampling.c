#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>   // input options handling
#include <sys/time.h> // time of parallel evaluation

#include "integral.h"

// codes for parsing arguments
enum ARG_OPTIONS{
  ARG_VAL_HELP,
  ARG_VAL_NGRID,
  ARG_VAL_NSHIFT,
  ARG_VAL_NPHI,
  ARG_VAL_FREQ_MAX,
  ARG_VAL_FREQ_MIN,
  ARG_VAL_NFREQ,
  ARG_VAL_OUTPUT,
};

void print_usage(FILE* stream, char* program_name) {
  fprintf(stream, "Usage: %s --ngrid=[] --nshift=[] --nphi=[] --freq_min=[] --freq_max=[] --nfreq=[] --output=[]\n", program_name);
  fprintf(stream,
      "   -h --help            Display usage information.\n"
      "   --ngrid              number of the grid-points in [-1,1] (including endpoints) - must be even number.\n"
      "   --nshift             number of shifts per projection.\n"
      "   --nphi               number of projections.\n"
      "   --freq_min           minimal frequency.\n"
      "   --freq_max           maximal frequency.\n"
      "   --nfreq              number of frequncies between [min_freq, max_freq) (excluding endpoint).\n"
      "   --output             output filename.\n"
  );
}


/*
  Single thread code for computation of matrix for the weighted Radon transform with cosine weight.
  Necessary parameters for run :
  "   --ngrid              number of the grid-points in [-1,1] (including endpoints) - must be even number.\n"
  "   --nshift             number of shifts per projection.\n"
  "   --nphi               number of projections.\n"
  "   --freq_min           minimal frequency.\n"
  "   --freq_max           maximal frequency.\n"
  "   --nfreq              number of frequncies between [min_freq, max_freq) (excluding endpoint).\n"
  "   --output             output filename.\n"

  Domain of definition of the matrix - pixels inside the centered unit distk.
  Integration along lines is given by Simpson's rule. Elements of the matrix are stored in a CSV format
  separated by newline - "\n" in the following order:

  for pixel do
    for frequency do
      for phi do
        for shift do
           ...
        endfor
      endfor
    endfor
  endfor

  The reading should respect the above format.

  Future plans: add support of library of weight functions. Maybe it is possible to inject lambda functions (in C++ version)
  when weight function is analytic.

*/

int main(int argc, char ** argv) {

  if (argc == 1) {
     fprintf(stderr, "%s: arguments required. Try '-h' or '--help' for help.\n", argv[0]);
     exit(EXIT_FAILURE);
  }

  char *program_name = argv[0];

  int next_option;
  const char* short_options = "h";
  const struct option long_options[] = {
    {"help", no_argument, NULL, ARG_VAL_HELP},  // print help
    {"ngrid", required_argument, NULL, ARG_VAL_NGRID}, // number of the grid-points in [-1,1] (including endpoints)
    {"nshift", required_argument, NULL, ARG_VAL_NSHIFT}, // number of shifts per projection
    {"nphi", required_argument, NULL, ARG_VAL_NPHI}, // number of directions
    {"freq_max", required_argument, NULL, ARG_VAL_FREQ_MAX}, // minimal frequency
    {"freq_min", required_argument, NULL, ARG_VAL_FREQ_MIN}, // maximal frequency
    {"nfreq", required_argument, NULL, ARG_VAL_NFREQ}, // number of frequncies between [min_freq, max_freq] (including endpoints)
    {"output", required_argument, NULL, ARG_VAL_OUTPUT}, // output filename
    {NULL, 0, NULL, 0}				/* Required at the end of array */
  };

  char* output_filename;
  int ngrid, nshift, nphi, nfreq;
  double freq_max, freq_min;

  //------------------------------ Reading options-------------------------------//
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, NULL);

    switch (next_option) {
      case 'h':
       print_usage(stdout, program_name);
       exit(EXIT_SUCCESS);

      case ARG_VAL_HELP:
	     print_usage(stdout, program_name);
  	   exit(EXIT_SUCCESS);

      case ARG_VAL_NGRID:
	     ngrid = atoi(optarg);
	     break;

      case ARG_VAL_NSHIFT:
	     nshift = atoi(optarg);
	     break;

      case ARG_VAL_NPHI:
	     nphi = atoi(optarg);
	     break;

      case ARG_VAL_FREQ_MIN:
       freq_min = atof(optarg);
	     break;

      case  ARG_VAL_FREQ_MAX:
       freq_max = atof(optarg);
       break;

      case ARG_VAL_NFREQ:
       nfreq = atoi(optarg);
       break;

      case ARG_VAL_OUTPUT:
       output_filename = optarg;
       break;

      case ':': /* missing argument */
	     fprintf(stderr, "%s: option '-%c' requires an argument. Aborting.\n", argv[0], optopt);
	     exit(EXIT_SUCCESS);

      case '?': /* help */
	     print_usage(stderr, program_name);
       exit(EXIT_SUCCESS);
      case -1: /* end of longopts */
	     break;

      default:
	     printf("?? getopt returned character code 0%o??\n", next_option);
	     exit(EXIT_SUCCESS);
    }
  } while (next_option != -1);

  if (optind < 8) {
    printf("Not enough arguments to run. Try '-h' or '--help' for help.\n");
    exit(EXIT_SUCCESS);
  }

  const int npixels = (ngrid-1)*(ngrid-1);

  double phi[nphi], shift[nshift], freqs[nfreq];

  const double delta_grid = 2.0 / (ngrid - 1); // size of a pixel
  const double delta_phi = 2.0*M_PI / nphi;
  const double delta_shift = 2.0 / (nshift + 1);
  const double delta_freq = (freq_max - freq_min) / nfreq;

  int i_pixel, i_freq, i_phi, i_shift;

  // init phi, shift, frequency arrays
  int i_counter;
  for (i_counter = 0; i_counter < nphi; ++i_counter) {
    phi[i_counter] = i_counter * delta_phi;
  }
  for (i_counter = 0; i_counter < nshift; ++i_counter) {
    shift[i_counter] = -1.0 + (i_counter + 1) * delta_shift;
  }
  for (i_counter = 0; i_counter < nfreq; ++i_counter) {
    freqs[i_counter] = freq_min + i_counter * delta_freq;
  }


  FILE *foutput;
	foutput = fopen(output_filename, "w");

  // iterate over pixels
  for (i_pixel = 0; i_pixel < npixels; ++i_pixel) {

      // iterate over frequencies
      for (i_freq = 0; i_freq < nfreq; ++i_freq) {

          // iterate over projections
          for (i_phi = 0; i_phi < nphi; ++i_phi) {

              // iterate over shifts
              for (i_shift = 0; i_shift < nshift; ++i_shift) {

		              // ray integral for function=pixel, this is matrix element
		              double element = matrix_element(ngrid,
                                                  npixels,
                                                  i_pixel,
                                                  delta_grid,
                                                  phi[i_phi],
                                                  shift[i_shift],
                                                  freqs[i_freq]);

		              //write to file
		              fprintf(foutput, "%lf\n", element);
              }
          }
      }
  }
	fclose(foutput);

  exit(EXIT_SUCCESS);
}
