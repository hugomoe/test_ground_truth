/*
 * zoom in an image by a given integer factor, using zero-padding
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.c"
#include "fft_zoom.h"
#include "parameters.h"
#include "aux_fun.h"



int main(int argc,char *argv[]){

	if (argc != 4) {
		printf("usage :\n\t[image_in.png] [image_out.png] z\n"); 
		return 1;
	}


	
//read the input
	char *filename_in = argv[1];
	char *filename_out = argv[2];
	int z = strtol(argv[3],NULL,10);

	float *img;
	int w,h,pd;

	img = iio_read_image_float_vec(filename_in, &w, &h, &pd);

//zoom-in
	float *img_zoomed = malloc(WOUT*z*HOUT*z*pd*sizeof(float));
	//symmetrize the input
	float *img_sym = malloc(2*w*2*h*pd*sizeof(float));
	int i_sym, j_sym;
	for(int l=0;l<pd;l++){
		for(int i=0;i<2*w;i++){
			i_sym = i-w/2;
			if(i_sym<0){i_sym = -1-i_sym;}
			else if(i_sym>w-1){i_sym = 2*w-1-i_sym;}
			for(int j=0;j<2*h;j++){
				j_sym = j-h/2;
				if(j_sym<0){j_sym = -1-j_sym;}
				else if(j_sym>h-1){j_sym = 2*h-1-j_sym;}
				img_sym[pd*(i+2*w*j)+l] = img[pd*(i_sym+w*j_sym)+l];
			}
		}
	}
	free(img);
	//zoom-in by zero-padding
	float *img_sym_zoomed = malloc(2*WOUT*z*2*HOUT*z*pd*sizeof(float));
	zoom_in(img_sym,2*w,2*h,pd,2*WOUT*z,2*HOUT*z,img_sym_zoomed);
	free(img_sym);
	//erase the symmetrization
	for(int l=0;l<pd;l++){
		for(int i=0;i<WOUT*z;i++){
			for(int j=0;j<HOUT*z;j++){
				img_zoomed[pd*(i + WOUT*z * j) + l] = img_sym_zoomed[pd*(i+WOUT*z/2 + 2*WOUT*z * (j+HOUT*z/2)) + l];
			}
		}
	}
	free(img_sym_zoomed);

//output
	iio_save_image_float_vec(filename_out,img_zoomed,WOUT*z,HOUT*z,pd);
	free(img_zoomed);
	return 0;
}
