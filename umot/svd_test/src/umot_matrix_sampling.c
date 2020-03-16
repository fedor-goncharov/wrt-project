#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integral.h"


/*
 * The program computes matrix for the exponential Radon transform for
 * limited angle data and for multiple exponential factors
 * (i.e., frequencies in the weight)
 *
 * Weight = cos(2pi*freq*t)
*/


const char* output_filename = "umot_matrix_output_cosine.csv";



int main(int argc, char * argv[]) {



  int nshift = 65;
  int nphi   = 15;
  int nfreq  = 20;
  int ngrid  = 33;	// number of frequencies used (maximal frequency should be 1/4 of the Nyquist frequency)
  double max_angle = M_PI/4.0; // max angle pi/4
  int npixels = (ngrid-1) * (ngrid-1);

  double phi[nphi], shift[nshift], frequencies[nfreq];
  double delta_grid = 2.0 / (ngrid-1);
  double delta_phi = max_angle / nphi;
  double delta_shift = 2.0*sqrt(2.0) / (nshift - 1);

  double nyquist_freq = 1.0/delta_grid; // Nyquist frequency
  double delta_freq = (nyquist_freq / 4.0) / (nfreq-1); // frequencies are
  //taken in the interval [0, 0.25*Nyquist-frequency] that is full resolution
  //is not achievable

  int i_pixel, i_phi, i_shift, i_freq;
  int i_counter;

  for (i_counter = 0; i_counter < nphi; ++i_counter) {
    phi[i_counter] = i_counter * delta_phi;
  }
  for (i_counter = 0; i_counter < nshift; ++i_counter) {
    shift[i_counter] = -sqrt(2.0) + i_counter * delta_shift;
  }
  for (i_counter = 0; i_counter < nfreq; ++i_counter) {
    frequencies[i_counter] = i_counter * delta_freq;
  }


  FILE *foutput;
	foutput = fopen(output_filename, "w");


  // iterate over pixels in the image
  for (i_pixel = 0; i_pixel < npixels; ++i_pixel) {

      //iterate over projections
      for (i_phi = 0; i_phi < nphi; ++i_phi) {

          //iterate over shifts in one projection
          for (i_shift = 0; i_shift < nshift; ++i_shift) {

              //iterate over frequencies in the weight
              for (i_freq = 0; i_freq < nfreq; ++i_freq) {

		              //compute ray integral of the function which gives an element of the matrix
		              double element = umot_matrix_element(ngrid, npixels, i_pixel, delta_grid,
                                                       phi[i_phi],
                                                       shift[i_shift],
                                                       frequencies[i_freq]);

		              //write element to file
		              fprintf(foutput, "%lf\n", element);
              }
          }
      }
  }
	fclose(foutput);

  exit(EXIT_SUCCESS);
}
