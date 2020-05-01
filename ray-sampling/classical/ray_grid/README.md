## General description

Program computes ray integrals of a test-function given on a uniform  
rectangular grid in 3D in a slice-by-slice, parallel-beam geometry.

Grid is parameterized by 2 parameters - **npixels**, **nslices**, that is the  
grid is represented as a three dimensional array (npixels, npixels, nslices)  
(one may think of coordinates XYZ, where grid has **nslices** in z-direction).  
Test-function is given by its values on a grid in a separate binary file,  
such that reading a linear array of floats corresponds to the following
order:

		`for z in nslices
		  for x in npixels
		   for y in npixels
		 	 test-function(x,y,z)
		   endfor
		  endfor
		 endfor`

In each z-slice function must be supported only in a disk inscribed in the  
square (npixels, npixels).

In each z-slice **nphi** projections are taken uniformly spread over [0, 2pi].  
Each projection consists of **nshift** parallel rays varying uniformly in  
[-**radius**, **radius**], where **radius** is the length of the side  
of the square slice.

Integration along each ray is performed via Simpson's rule and function  
is assumed to be continuous linear between nodes - nodes of the grid.  
By default **2D linear interpolation** is used when computing values  
of the function between the nodes.

Output is stored in a separate binary file in the following order:  

	 `for z in nslices
	   for shift in nshifts
	    for phi in nphi
		  ray-transform(f)(z, shift, phi)
		endfor
	   endfor
	  endfor`

Code can be run in **parallel**, where work for z-slices is splitted  
uniformly over OpenMP threads.

## Input arguments
         -h --help     no_arg       Display usage information.
	     -p --nphi     integer      Number of projections per slice.
	     -s --nshift   integer      Number of shifts per projection.
	     -g --npixels  integer      Number of pixels per dimension.
	     -z --nslices  integer      Number of slices in z-direction.
	     -r --radius   float        Radius of the support in cm in a slice.
	     -i --input    filename     Path to file with test-function data.
	     -o --output   filename     Path to output file.
	     -n --nthreads integer      Number of OpenMP threads for parallelization.
	     -v --verbose  no_arg       Show extra information during computations.

## Output
	binary file named with argument of '-o', '--output'

## Compilation
	make
	make clean

## Requirements
	gcc, openmp

## Usage / Examples
>:$ ray_grid -p 1024 -s 512 -g 512 -z 128 -r 1.0 -i test-func.bin -o ray-tr.bin -n 8
