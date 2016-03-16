# test_ground_truth

The main program is test_ground_truth.c
It takes a huge blurred image (of dimension z*WOUT, z*HOUT for some z integer) and an homography to compare several warping methods.


The only purpose of zoom.c is to zoom-in an image to obtain, once for all, a huge blurred image that will be test_ground_truth's input. Thus zoom.c can be deleted without altering the main program's functionment, if the user already has huge blurred images.
