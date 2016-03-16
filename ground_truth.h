//contain functions for bilinear interpolation and gaussian zoom-out
//warning : this file is totally different from other ground_truth.h that can be found in other repositories



#include "parameters.h"
#include "aux_fun.h"



//extrapolate by a constant value
float getsample_cons(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h)
		return 0;
	if (l >= pd)
		l = pd - 1;
	return x[(i+j*w)*pd + l];
}

// auxiliary function for bilinear interpolation
static float evaluate_bilinear_cell(float a, float b, float c, float d,
							float x, float y)
{
	return a * (1-x) * (1-y)
	     + b * ( x ) * (1-y)
	     + c * (1-x) * ( y )
	     + d * ( x ) * ( y );
}

// bilinear interpolation
float bilinear_interpolation_at(float *x, int w, int h, int pd,
		float p, float q, int l)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = getsample_cons(x, w, h, pd, ip  , iq  , l);
	float b = getsample_cons(x, w, h, pd, ip+1, iq  , l);
	float c = getsample_cons(x, w, h, pd, ip  , iq+1, l);
	float d = getsample_cons(x, w, h, pd, ip+1, iq+1, l);
	return evaluate_bilinear_cell(a, b, c, d, p-ip, q-iq);
}



int zoom_out(float *img,int w,int h,int pd,int w_f,int h_f,float *img_f){	
	int zoom = w/w_f;

	//zoom-out by convolution with a gaussian kernel
	double sigma = 0.8 * zoom;
	int taps = ceil(6*sigma);
	double *gauss = malloc((2*taps+1)*sizeof(double));
	double tot = 0;
	for(int k=-taps;k<=taps;k++){
			tot += (gauss[k+taps] = exp(-(pow(k,2))/(2*pow(sigma,2))));
	}
	for(int k=-taps;k<=taps;k++){gauss[k+taps] /= tot;}


	float *img_aux = malloc(w_f*h_f*zoom*pd*sizeof(float)); //size w_f x h_f*zoom
	if(img_aux==NULL){
        printf("img_aux is too large");
        exit(1);
    }
	int idx, i0, j0;
	for(int l=0;l<pd;l++){
		//horizontal convolution
		for (int i = 0; i < w_f; i++){
			for (int j = 0; j < h_f*zoom; j++){
				float v = 0;
				for(int k=-taps;k<=taps;k++){
					i0 = good_modulus(i*zoom+k,w_f*zoom);
					v += gauss[k+taps]*img[(i0+j*zoom*w_f)*pd+l];
				}
				idx = l+pd*(w_f*j+i);
				img_aux[idx] = v;
			}
		}
		//vertical convolution
		for (int i = 0; i < w_f; i++){
			for (int j = 0; j < h_f; j++){
				float v = 0;
				for(int k=-taps;k<=taps;k++){
					j0 = good_modulus(j*zoom+k,h_f*zoom);
					v += gauss[k+taps]*img_aux[(i+j0*w_f)*pd+l];
				}
				idx = l+pd*(w_f*j+i);
				img_f[idx] = v;
			}
		}
	}

	free(img_aux);
	free(gauss);
	return 0;
}
