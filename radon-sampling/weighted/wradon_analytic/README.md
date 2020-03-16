## General information

This program computes classical Radon transforms in 3D (along 2D planes) of analytic functions whose prototype should be realized in "test_function.c".  

File "test_function.c" contains a template of such realization. Note that only the function with name "test_function" will be used for computations. After realizaiton of your test function you have to compile the code so it can be used.  

## Requirements 

The programs here are designed to work under Unix operating systems.  
To compile the project on your computer you need to have installed:  

GCC compiler, OpenMP libraries, GNU GSL libraries (+2.5)

## Compilation / Installation
  1) Go to 'src' directory:  
        ```
          cd src
        ```
  2) Open Makefile and set the name of the output file:
        open Makefile in any text editor and set
        ```
          ONAME=(output name of your binary)
        ```
  
  3) Run Makefile
      ```
        make install
      ```
  4) Clean directory from object files (optional):
  
      ```
        make clean 
      ```
  If you want to generate data for other test-function then you have to change the file
  'test_function.c' and repeat steps (1-4), possibly setting a new name in 'ONAME'.
  
## Usage / Examples

(binary) -h --help &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: Display usage information  
(binary) -p --parameters filename &nbsp;&nbsp;&nbsp;&nbsp;: Read parameters of the grid in Radon space from config file  
(binary) -o --output filename &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: Write output data to a file  
(binary) -n --nthreads number &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: Use (number) of OpenMP threads for parallelization  

### Config file requirements

The purpose of the config file for parameters -p (--parameters) is to provide to the program information about the grid  
in Radon space. Recall that data in Radon space is parametrized by directions on the unit sphere (two angles: latitude, longitude) and  by shifts along each direction (shifts vary uniformly in [-1,1]). Also, though the test-function is
given by an analytical expression its integration over planes necessarily requires discretization which is also must be 
given in the config file. More precisely, a test-function is always assumed to be supported in a unit ball in 3D which 
lies inside the unit cube [-1,1]^3. So the parameter to be specified is a number of points per dimension in the uniform rectangular grid inside this unit cube.

1) The first line contains a number of longitude angles which are positioned uniformly along [0,2pi].  
2) The second line contains a number of latitude angles which are positioned according to Gauss-Legendre quadrature 
   rule on [0, pi].  
3) The third line contains a number of shifts which are positioned uniformly along [-1,1].  
4) The fourth lines containes a number of points per dimension in a unit cube.  

Some comments are allowed after each line, however the length of each line should not exceed 128 symbols.

### Example of a config file

> 256			: number of longitude angles  
> 128			: number of latitude angles  
> 129			: number of steps per fixed direction  
> 129			: number of points on the grid per dimension

### Output

The output stored in a specified output file (-o --output) in a CSV format in the following order:  
**[shift], [longitude angle], [latitude angle], [value of Radon transform]**  

Example:  
> -1.000000, 0.000000, 3.122878, -0.009362  
> -1.000000, 0.000000, 3.098635, -0.009187  
> -1.000000, 0.000000, 3.074249, -0.008878

### Examples of test-functions

One can find in folder 'test_functions' some examples of realizations of some simple test-functions.  
In order to try them, rename any of those files to 'test_function.c' and copy them to the main directory of 'radon_analytic'.     Then proceed with steps in 'Compilation / Installation' in order to obtain a compiled binary for a given
test-function. 


