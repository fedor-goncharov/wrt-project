### Radon inversion in 2D/3D

Octave/Matlab scripts for inversion of classical Radon transforms in 2D/3D via projection theorem. 
Generally speaking, functions here take an input a CSV file containing Radon transforms and inverts them via 
projection theorem (1D Fourier transform + 3D inverse Fourier transform).

  * **RtFt_2d.m** - Fourier transform of Radon transforms in 2D along shift variable
  * **RtFt_3d.m** - Fourier transforms of Radon transforms in 3D along shift variable
  * **lgwt.m** - script for computation of Gaussian nodes in interval [-1, 1] (used in RtFt_3d.m)
  * **nfft_reconstruct_3d.m** - application of inverse Fourier transform to the data returned by RtFt_3d.m


#### Dependencies

   You must have Matlab/Octave interface for [NFFT](https://www-user.tu-chemnitz.de/~potts/nfft/) installed. 
   The link to this interface (mex-libraries) can be found [here](https://www-user.tu-chemnitz.de/~potts/nfft/download.php).

   To install NFFT download Matlab/Octave binaries from the link given above (precompiled mex-files), unzip it 
   and configure the path your Matlab/Octave to these files.

#### Usage

 Before to use RtFt_3d.m load lgwt.m to your Matlab/Octave.

 * **RtFt_3d.m**
 
       [nodes, values, jacobian_weights] = RtFt_2d(filename, nphi, ntheta, nshift, rsupp, padding_coeff = 4)
       
       Script reads Radon transforms in 3D from file and performs 1D Fourier transforms 
       along shift variable. 
       
       Returns values
         nodes            : points in 3D frequency space where Fourier transforms is evaluated (size Nx3)
         values           : values of Fourier transform in nodes (size Nx1complex)
         jacobian_weights : volume associated to each node in frequency space

       Usage of the script
         filename          : file where the data is stored in CSV format
                             Data is expected in the following format : "[shift], [phi], [theta], [value]\n",
                             Variables 'shift, phi, theta' vary in the following order : 
                 
                                 for (shift) 
                                    for (phi) 
                                       for (theta)
                                        ....
                                       end
                                    end
                                 end
                                      
         nphi              : number of projections in azimuth angle [0, 2*pi)
         ntheta            : number of projections in polar angle (0, pi)
        
                             Angles 'phi' are uniform on the circle and angles 'theta' 
                             correspond to Gaussian quadrature points, i.e. theta_j = arccos(t_j), 
                             (t_j, j = 1, ntheta) - Gauss-Lebato points on [-1, 1]. 

         nshift            : number of hyperplanes per one direction
                             Shifts are uniform along [-1,1]
         rsupp             : radius of the support of the test function
         padding_coeff (default=4) : parameter to padd Radon transforms with zeros along shift
 
 * **RtFt_2d.m**
 
       [nodes, values, jacobian_weights] = RtFt_2d(filename, nphi, nshift, rsupp, padding_coeff = 4)
       
       Script reads Radon transforms in 2D from file and performs 1D Fourier transforms 
       along shift variable. 
       
       Returns values
         nodes            : points in 2D frequency domain where Fourier transforms is evaluated (size Nx3)
         values           : values of Fourier transform in nodes (size Nx1complex)
         jacobian_weights : volume associated to each node in frequency space

       Usage of the script
         filename          : file where the data is stored in CSV format
                             Data is expected in the following format : "[shift], [phi], [value]\n",
                             Variables 'shift, phi' vary in the following order : 
                 
                                 for (shift) 
                                    for (phi) 
                                        ...
                                    end
                                 end
                                      
         nphi              : number of projections in azimuth angle [0, 2*pi)
                             Angles 'phi' are uniform on the circle.

         nshift            : number of hyperplanes per one direction
                             Shifts are uniform along [-1,1].
         rsupp             : radius of the support of the test function
         padding_coeff(default=4) : parameter to padd Radon transforms with zeros along shift
 
 * **nfft_reconstruct_3d.m**
 
       test_function = nfft_reconstruct_3d(ngrid, nodes, values, jacobian_weights)
       
       Script performs reconstruction of a function from its Fourier transforms at 'nodes' 
       with values 'values'. Mathematically it works as a discretized version of (inverse) Fourier integral 
       in 3D. Result is given as a 3D matrix of size : ngrid x ngrid x ngrid, which, in turn, is a 
       grid on [-1,1)x[-1,1)x[-1,1). 

       Returns values 
         test_function : real-valued matrix of size (ngrid x ngrid x ngrid)

       Usage of the script
        ngrid : number of points in [-1.0, 1.0); ngrid must be even;
                Note that it is number of points on non-closed interval. 
                NFFT assumes that your signal is periodic, so the values on the missing edge points
                of the grid is can reconstructed from periodicity.  
        
        nodes            : matrix of size (Nnodesx3); these are the points 
                           in space where Fourier transform of signal is known. 
        values           : vector of size (Nnodesx1(complex)); these are the values of 
                           Fourier transform of the signal at nodes;
        jacobian_weights : volumes of cells related to nodes in the Riemann summ of discretized Fourier integral
