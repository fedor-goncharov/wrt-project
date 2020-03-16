# Weighted ray / Radon transforms in 3D

<p float="center">
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/k_comparison_output.gif" width="360" />
  <img src="https://github.com/fedor-goncharov/Weighted-ray-Radon-transforms-in-3D/blob/master/pictures/shepp_logan_reduction.gif" width="360" />
</p>

This is project is a part of my phd thesis conducted under the supervision of professor [Roman Novikov](http://www.cmap.polytechnique.fr/~novikov/). The thesis is now [here](http://www.theses.fr/2019SACLX029).

We aim to develop new inversion methods for weighted (generalized) Radon transforms in Euclidean space. 
The latter are of particular importance in various applications in the domain of inverse 
problems (e.g., in tomographies, geophysics). In particular, we work on methods which 
could be numerically more stable against the noise in tomographical data in SPECT. 

The very precise theoretical explanation of given algorithms can be found in [[2,3]](http://www.cmap.polytechnique.fr/~fedor.goncharov/publications.html).
Also my personal page is [here](http://www.cmap.polytechnique.fr/~fedor.goncharov/).

## Some information about ray / Radon transforms

Here it will be the link to a pdf file with minimal information about the subject. 

## Sructure of the project

Here I briefly explain for what each program is intended. Details about their input, output, parameters, usage, etc., 
you can find further in respective folders.

  * #### radon_inversion (Matlab / Octave)
        Octave/Matlab scripts for inversion of classical Radon transforms in 2D/3D using projection theorem
        or backprojection algorithms.
  
  * #### radon_analytic (C) 
        computations of Radon transforms in 3D of a test-function with compact support in 
        a three-dimensional unit ball; the main feature is that the expression of the test-function 
        must be given by an analytical expression in a separate C file
  
  * #### radon_grid (C)
        computations of Radon transforms in 3D of a test-function with compact support in 
        a three-dimensional unit ball; here the test-function is given by its values 
        on a discrete grid in a three-dimensional cube 
  
  * #### ray_analytic (C) 
        computations of ray transforms in 3D in a layer-by-layer sampling scheme of a test-function 
        with compact support in a three-dimensional unit ball (see also the README.md inside project 
        for information about 'layer-by-layer sampling scheme');
        the test-function must be given by an analytical expression in a separate C file
  
  * #### ray_grid (C)
        computations of ray transforms in 3D in a layer-by-layer sampling scheme of a test-function 
        with compact support in a three-dimensional unit ball (see the README.md inside project 
        for information about 'layer-by-layer sampling scheme');
        test-function is given by its values on a discrete grid in a three-dimensional cube
  
  * #### radon_reduction (C)
        reduces the data given by ray transforms in 3D in layer-by-layer sampling scheme to the 
        data given by Radon transforms in 3D
  * #### wradon_inversion/iterative (Matlab / Octave)
        iterative inversions of weighted Radon transforms in 3D. Takes on input initial starting point 
        and data about the weight (coefficients in spherical harmonics expansion) and solves iteratively 
        the related integral equation (see link [3] given above).

## Future plans

  * The programs above are designed for computations of the classical ray/Radon transforms in 3D and their inversion, however, 
  the goal of the project is to implement some iterative inversion algorithms for weighted Radon transforms (see [[2,3]](http://www.cmap.polytechnique.fr/~fedor.goncharov/publications.html)). 
  
  * The current version of the programs can be used to test Chang-type approximate inversions (described 
    precisely in [[2]](http://www.cmap.polytechnique.fr/~fedor.goncharov/publications.html)) and Kunyansky-type inversions
      
  * Add more options for computations of classical Radon transforms (for example sampling grids: 
      right now in Radon space the grid is 'Gauss-uniform', however, the set of directions on the sphere can have many 
      different realizations (possibly with different properties of stability of algorithms))
      
  * Add options for tests with noise as for classical Radon transforms and for weighted Radon transforms (like SPECT, PET)
      
      
      

