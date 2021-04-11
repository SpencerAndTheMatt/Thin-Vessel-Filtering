# Thin-Vessel-Filtering
Files for my third year university project based upon thin vessel filtering. The objective is to obtain a crisp image of the blood vessels following an MRI scan. This will be done by using the second derivative of the Gaussian, convoluted with the image of the blood vessels.

The 2D filter is now complete and works well on test data. It struggles to work on real data, however I suspect this is due to the resolution of the images used (as they weren't official MRI scans, just screenshots of scans). I came to this conclusions as another filter also did not work on the real image.

*Big update*
The filter is complete, and can be seen in Filter3D.m

This filter works extremely effectively upon MRI data, and also the test data found in 'test3Dfilter7' and 'test3Dfilter8'

Some changes made from the last iteration
- Sorting algorithm now used a logical index sort, i.e [~, I] = sort('abs', [], 4) which then uses this to sort the eigenvalues
- A new constraint was added, namely if ev1/ev2 < 0.25, vesselness = 0. This distinguishes between a plate/vessel structure.
- Memory usage is extremely high (at least ~ 80gb for an MRI image). Due to the nature of MRI data, I don't think this is to be avoided.
- The dice scores and SNR is useful, however it was found that the dice score probability was meaningless when applied to MRI data.


Summation:
Filter3D.m is a matlab file which takes a 3 dimensional array, and computes likelihood of vesselness at every point within the image at a single scale factor. It uses Hessian based eigenvalue decomposition for this task. To make it multiscale, see example in the test3Dfilter files, but essentially make an array to store data, call a loop at different scale factors, then compute the maximum of all values within the holder array via the fourth dimension.

*Note*
The filter will cut out vessels at intersections (As seen on test files). This is just filter behaviour. It is expected, and due to the nature of how the filter works, it cannot be helped. Essentially, 3 crossed vessels do not look like a vessel at the intersection.
