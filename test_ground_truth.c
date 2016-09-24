/**
  * warp a huge image to a small image by three methods and compare the results
  * first method : zoom-out and warp using decomposition
  * second method : zoom-out and warp using ripmap
  * third method : bilinear warping then zoom-out (ground truth)
  * 
  * the input must be of size z*WOUT * z*HOUT for some integer z
  *
  * the input 3x3 matrix corresponds to the homography that decomp and ripmap will apply on the zoomed-out image
  * (i.e. from a WOUT*HOUT image to a WOUT*HOUT image)
  * it does NOT correspond to the homography that will be applied to the large blurred image
  * 
  * 
  * 
  * compilation
  * 	gcc-5 -fopenmp -O3 test_ground_truth.c -lfftw3 -ltiff -ljpeg -lpng -o test_ground_truth
  * 	gcc -std=c99 -fopenmp -O3 test_ground_truth.c -lfftw3 -ltiff -ljpeg -lpng -o test_ground_truth
  * usage
  * 	./test_ground_truth [image.png] a b p c d q r s t
  * 	(input image and homography matrix)
  *		or
  * 	./test_ground_truth [image.png] a b p c d q r s t i1 j1 i2 j2
  * 	(input image, homography matrix and two pixels between which the error will be computed)
  */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.c"
#include "decomp.h"
#include "ground_truth.h"
#include "ripmap.h"
#include "parameters.h"
#include "aux_fun.h"



int main(int argc,char *argv[]){

	if (argc != 11 && argc != 15) {
		printf("usage :\n\t[image.png] a b p c d q r s t\nor\n\t[image.png] a b p c d q r s t i1 j1 i2 j2"); 
		return 1;
	}


	
//read the input
	char *filename_in = argv[1];

	double H[3][3];
	H[0][0]=strtod(argv[2],NULL);
	H[0][1]=strtod(argv[3],NULL);
	H[0][2]=strtod(argv[4],NULL);
	H[1][0]=strtod(argv[5],NULL);
	H[1][1]=strtod(argv[6],NULL);
	H[1][2]=strtod(argv[7],NULL);
	H[2][0]=strtod(argv[8],NULL);
	H[2][1]=strtod(argv[9],NULL);
	H[2][2]=strtod(argv[10],NULL);

	int i1,j1,i2,j2; //coordinates of two opposite vertices of a sub-image of the output
	switch(argc){
		case 15:
		i1 = strtol(argv[11],NULL,10);
		j1 = strtol(argv[12],NULL,10);
		i2 = strtol(argv[13],NULL,10);
		j2 = strtol(argv[14],NULL,10);
		break;
		case 11:
		default:
		i1 = 0;
		j1 = 0;
		i2 = WOUT-1;
		j2 = HOUT-1;
		break;
	}

	float *img_big;
	int w,h,pd;

	img_big = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	if (h/HOUT != w/WOUT) {printf("image size is incompatible with WOUT and HOUT.\n");exit(1);}
	float *img = malloc(WOUT*HOUT*pd*sizeof(float));
	if(img==NULL){
		printf("@compare_with_truth.c @main : img is too large");
		exit(1);
	}
	zoom_out(img_big,w,h,pd,WOUT,HOUT,img);

//Ripmap
	float *img_rip = malloc(3*WOUT*HOUT*sizeof(float));
	if(img_rip==NULL){
        printf("@compare_with_truth.c @main : img_rip is too large");
        exit(1);
    }
	
	if(pd==3){
        apply_homo_ripmap(img,img_rip,WOUT,HOUT,WOUT,HOUT,H);
	}else{//suppose pd=1
        float *img3 = malloc(3*WOUT*HOUT*sizeof(float));
        for(int i=0;i<WOUT*HOUT;i++){
            for(int l = 0;l<3;l++){
                img3[3*i+l]=img[i];
            }
        }
        apply_homo_ripmap(img3,img_rip,WOUT,HOUT,WOUT,HOUT,H);
	}

//decomposition
	float *img_dec = malloc(WOUT*HOUT*3*sizeof(float));
	if(img_dec==NULL){
        printf("@compare_with_truth.c @main : img_dec is too large");
        exit(1);
    }

	if(pd==3){
        apply_homography_decomp(img,img_dec,WOUT,HOUT,WOUT,HOUT,H);
	}else{//suppose pd=1
        float *img3 = malloc(3*WOUT*HOUT*sizeof(float));
        for(int i=0;i<WOUT*HOUT;i++){
            for(int l = 0;l<3;l++){
                img3[3*i+l]=img[i];
            }
        }
        apply_homography_decomp(img3,img_dec,WOUT,HOUT,WOUT,HOUT,H);
	}
free(img);

//ground_truth
	float *img_warped_zoomed = malloc(pd*w*h*sizeof(float));
	if(img_warped_zoomed==NULL){
		printf("@compare_with_truth.c @main : img_warped_zoomed is too large");
		exit(1);
	}
	//naive warping
	double fzoom=(double)h/(double)HOUT;
	double HH[3][3];
	HH[0][0]=H[0][0];
	HH[0][1]=H[0][1];
	HH[0][2]=H[0][2]*fzoom;
	HH[1][0]=H[1][0];
	HH[1][1]=H[1][1];
	HH[1][2]=H[1][2]*fzoom;
	HH[2][0]=H[2][0]/fzoom;
	HH[2][1]=H[2][1]/fzoom;
	HH[2][2]=H[2][2];
	for(int l=0;l<pd;l++)
	for (int j = 0; j < h; j++)
	for (int i = 0; i < w; i++)
	{
		double p[2] ={i,j};
		
		apply_homography_1pt(p, HH, p);
		int idx = pd*(w * j + i)+l;
		img_warped_zoomed[idx] = bilinear_interpolation_at(img_big, w, h, pd, p[0], p[1], l);
	}

	float *img_grd = malloc(pd*WOUT*HOUT*sizeof(float));
	if(img_grd==NULL){
		printf("@compare_with_truth.c @main : img_grd is too large");
		exit(1);
	}
	zoom_out(img_warped_zoomed,w,h,pd,WOUT,HOUT,img_grd);
	free(img_warped_zoomed);



//output
	iio_save_image_float_vec("img_dec.png",img_dec,WOUT,HOUT,3);
	iio_save_image_float_vec("img_grd.png",img_grd,WOUT,HOUT,pd);
	iio_save_image_float_vec("img_rip.png",img_rip,WOUT,HOUT,3);



//error on a sub-image
	int i_tl, j_tl, //coordinates of the top left pixel of the sub-image
		i_br, j_br; //coordinates of the bottom right pixel of the sub-image
	if(i1<i2){
		i_tl = i1;
		i_br = i2;
	}else{
		i_tl = i2;
		i_br = i1;
	}
	if(j1<j2){
		j_tl = j1;
		j_br = j2;
	}else{
		j_tl = j2;
		j_br = j1;
	}

	//the sub-image must be a subset of the output image
	if(i_tl<0){i_tl=0;}
	if(j_tl<0){j_tl=0;}
	if(i_br>WOUT-1){i_br=WOUT-1;}
	if(j_br>WOUT-1){j_br=WOUT-1;}

	float l1_dg = 0., l2_dg = 0., //compare decomposition and ground truth
		l1_dr = 0., l2_dr = 0., //compare decomposition and ripmap
		l1_rg = 0., l2_rg = 0.; //compare ripmap and ground truth
	float diff;
	int idx, idxx;

	for(int l=0 ; l<pd ; l++)
		for(int i=i_tl ; i<=i_br ; i++)
			for(int j=j_tl ; j<=j_br ; j++){
				idx = pd*(i+WOUT*j) + l;
				idxx = 3*(i+WOUT*j) + l;

				diff = fabs(img_dec[idxx]-img_grd[idx]);
				l1_dg += diff;
				l2_dg += pow(diff,2);

				diff = fabs(img_dec[idxx]-img_rip[idxx]);
				l1_dr += diff;
				l2_dr += pow(diff,2);

				diff = fabs(img_rip[idxx]-img_grd[idx]);
				l1_rg += diff;
				l2_rg += pow(diff,2);
	}

	//normalize
	int N = (i_br-i_tl+1)*(j_br-j_tl+1)*pd; //size of the sub-image
    l1_dg = l1_dg/(float) N;
    l2_dg = sqrt(l2_dg/(float) N);
    l1_dr = l1_dr/(float) N;
    l2_dr = sqrt(l2_dr/(float) N);
    l1_rg = l1_rg/(float) N;
    l2_rg = sqrt(l2_rg/(float) N);

    //print
    printf("l1-error :\n\t\t\tGround Truth\tRipMap\n\tDecomposition\t%f\t%f\n\tRipMap\t\t%f\n\n",l1_dg,l1_dr,l1_rg);
    printf("l2-error :\n\t\t\tGround Truth\tRipMap\n\tDecomposition\t%f\t%f\n\tRipMap\t\t%f\n\n",l2_dg,l2_dr,l2_rg);

    free(img_dec);
    free(img_grd);
    free(img_rip);
	return 0;
}
