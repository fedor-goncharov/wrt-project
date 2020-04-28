## General information

Program computes Radon transforms in 3D (along planes) of a **test-function**
given on a uniform rectangular grid.

Rectangular grid is parameterized by 1 parameter - **npixels**, that is the
test-function on the grid is represented as a three dimensional array
**(npixels, npixels, npixels)**.

Test-function is given by its values on a grid in a separate binary file,
such that reading a linear array of floats corresponds to the following
order:

	for z in npixels
		for x in npixels
			for y in npixels
		 		test-function(x,y,z)
			endfor
		endfor
	endfor

Inside the grid cube test-function must me supported in the maximal inscribed
ball that touches borders of the cube.  

Grid of 2d-planes is characterised by three parameters: **shift**,
**phi**, **theta**, which stand for the usual coordinates of the Radon
transform.

* **shift** stands for distance between the plane and the origin,
	shifts vary uniformly in [-**radius**, **radius**] with **nshifts** in
	total.
* **phi** stands for the equatorial angle of the normal to the plane,
 angles vary uniformly in [0, 2pi) with **nphi** angles in total.
* **theta** stands for the azimuth angle of the normal to the plane,
vary in (0, pi) with **ntheta** in total (generated according
Gauss-Lebato quadrature rule : theta_k = arccos(t_k), t_k from [-1, 1]
by the quadrature of degree **ntheta**)

Integration over the plane is performed via Simpson's rule and function is
assumed to be linear continuous between the nodes (i.e., linear interpolation)
is used between nodes.

Output is stored in a separate binary file in the following order:

	for shift in shift_array
		for phi in phi_array
			for theta in theta_array
				radon-transform(f)(shift, phi, theta)
			endfor
		endfor
	endfor


 Code can be run in **parallel**, where work for shifts is splitted
 uniformly over OpenMP threads.


## Input arguments

	 -h --help      no_arg     Display this usage information.
	 -p --nphi      integer    Number of phi angles.
	 -t --ntheta    integer    Number of theta angles.
	 -s --nshift    integer    Number of shifts per direction.
	 -r --radius    float      Radius of the test-function support in cm.
	 -g --npixels   integer    Number of pixels per dimension for the test-function.
	 -i --input     string     Path to input file with test-function.
	 -o --output    string     Path to output file.
	 -n --nthreads  integer    Number of OpenMP threads for parallelization.
	 -v --verbose   no_arg     Print extra information during computations.

## Compilation
	make
	make clean

## Requirements
	gcc, openmp

## Usage/Examples

>:$ radon_grid -p 128 -s 256 -t 128 -g 256 -r 1.0 -i test-func.bin -o ray-tr.bin -n 8 -v
