# Thin-Vessel-Filtering
Files for my third year university project based upon thin vessel filtering. The objective is to obtain a crisp image of the blood vessels following an MRI scan. This will be done by using the second derivative of the Gaussian, convoluted with the image of the blood vessels.

Note: The file 'new_file_second_differential.m' is no longer in use.

Completed tasks:
- Second differentials of Gaussians with respect to dx^2, dy^2 and dxdy are plotted.
- Gaussians are convoluted with a fake image consisting of 5 columns of ones down the middle, a 5 x 5 square of ones in the corner, and zeros elsewhere.
- Eigenvalues are computed using the Hessian Matrix.
- Filter equation is complete and ready to be used.

Current tasks:
- Scale factor is not working, this must be fixed.
- Eigenvalue filtering equation does not appear to be filtering. Suspect it is a problem with the scale factor.
