#include <stdio.h>
#include "integral.h"


struct spect_threadf_argument {
	int thread_ID;
	char* output_filename;
	int z_slice_min, z_slice_max;	// interval of z-slices [z_slice_min, z_slice_max-1]
	int nphi, nshift, nslice, ngrid;
	double ***func_values, ***att_values;
	double *phi, *shift;
};

void* spect_ray_threadf(void* thread_argument) {

	
	struct spect_threadf_argument* pargument = (struct spect_threadf_argument*)(thread_argument);

	int thread_ID         = pargument -> thread_ID;
	char *output_filename = pargument -> output_filename;

	int	z_slice_min = pargument -> z_slice_min,
	 	z_slice_max = pargument -> z_slice_max,
	 	ngrid       = pargument -> ngrid,
	 	nphi        = pargument -> nphi, 
	 	nshift      = pargument -> nshift, 
	 	nslice      = pargument -> nslice;

	double ***func_values = pargument -> func_values,
		   ***att_values  = pargument -> att_values,
		   *phi           = pargument -> phi,
		   *shift         = pargument -> shift;

    fprintf(stdout, "  POSIX thread %d. Z slices: start %d, end %d, size %d\n", thread_ID, z_slice_min, z_slice_max-1, 
    	z_slice_max - z_slice_min);
        
	
        //open ouptut file
    FILE *foutput;
    foutput = fopen(output_filename, "w");
	int i_z_slice, i_phi, i_shift;
	
	double z_slice;
    const double dz_slice = 2.0 / (nslice - 1);
        
    //iterate ray integrals over z-slices
    for (i_z_slice = z_slice_min; i_z_slice < z_slice_max; ++i_z_slice) {      
	    z_slice = (-1.0) + i_z_slice * dz_slice;
	    for (i_shift = 0; i_shift < nshift; ++i_shift) {
                for (i_phi = 0; i_phi < nphi; ++i_phi) {
	
	        		//spect ray integral
		    		double spect_rt = spect_ray_transform(func_values, 
		    										 att_values,
		    										 ngrid,
		    										 z_slice,
		    										 shift[i_shift], 
		    										 phi[i_phi]);
		    		fprintf(foutput, "%lf, %lf, %lf, %lf\n", z_slice, shift[i_shift], phi[i_phi], spect_rt);	// spect ray transform
                }
        }
    }
        
	fclose(foutput);
	printf("  POSIX thread %d. Job is done.\n", thread_ID);
	return NULL;
}