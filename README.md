# Thin-Vessel-Filtering
Files for my third year university project based upon thin vessel filtering. The objective is to obtain a crisp image of the blood vessels following an MRI scan. This will be done by using the second derivative of the Gaussian, convoluted with the image of the blood vessels.

The 2D filter is now complete and works well on test data. It struggles to work on real data, however I suspect this is due to the resolution of the images used (as they weren't official MRI scans, just screenshots of scans). I came to this conclusions as another filter also did not work on the real image.

The 3D filter is mostly complete. I believe it's as fast as I can get it working. The loops have been removed, and it's fairly barebones, i.e there are no default arguments or anything, and it only returns the vesselness, nothing else. It does not yet allow you to filter those results by an amount, however this would probably be done in the script calling this function anyway.

Note: The file 'new_file_second_differential.m' is no longer in use.

Completed tasks:
- Second differentials of Gaussians with respect to dx^2, dy^2 and dxdy are plotted.
- Gaussians are convoluted with a fake image consisting of 5 columns of ones down the middle, a 5 x 5 square of ones in the corner, and zeros elsewhere.
- Eigenvalues are computed using the Hessian Matrix.
- Filter equation is complete and ready to be used.
- Scale factor is not working, this must be fixed. FIXED
- Eigenvalue filtering equation does not appear to be filtering. Suspect it is a problem with the scale factor. FIXED
- Tested on test image, successful
- Tested on real image, not ideal results
-

Current tasks:
- Compute Dice score of 3D filter
- Test on cluster data
- See if the algorithm can be optimised (Currently I think it's as good as it's going to get).


