/***********************************************************
 * sample_discrete.c
 *   Refer to sample_discrete.m for documentation.
 * 
 * Author: David Ross, 2004.
 * $Id$
 ***********************************************************/

#include "mex.h"
#include <math.h>
#include <stdlib.h>

void sample_discrete(double *prob, int length, double *samples, 
		     int num_samples);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
		 const mxArray *prhs[]) {
   int Nx, Ny, num_samples;
   double *probability, *samples;

   /* check for the proper number of arguments */
   if (nrhs != 2) {
      mexErrMsgTxt("wrong number of input arguments");
   }
   if (nlhs != 1) {
      mexErrMsgTxt("wrong number of output arguments");
   }

   /* get size of probability */
   Ny = mxGetM(prhs[0]);
   Nx = mxGetN(prhs[0]);

   /* do some simple checks on probability */
   if (!mxIsDouble(prhs[0]) || Nx < 1 || Ny < 1) {
      mexErrMsgTxt("probability has invalid size");
   }
   if (Nx != 1 && Ny != 1) {
      mexErrMsgTxt("probability must be a vector");
   }

   /* get probability */
   probability = mxGetPr(prhs[0]);


   /* make sure the # of samples to get is okay */
   if (!mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || 
       mxGetN(prhs[1]) != 1) {
       mexErrMsgTxt("num_samples must be a scalar");
   }

   /* get the number of samples */
   num_samples = (int)mxGetScalar(prhs[1]);
   
   if (num_samples < 1) {
      mexErrMsgTxt("num_samples must be at least 1");
   }
   

   /* create the result */
   plhs[0] = mxCreateDoubleMatrix(num_samples, 1, mxREAL);
   samples = mxGetPr(plhs[0]);
   
   /* do it */
   sample_discrete(probability, Nx*Ny, samples, num_samples);
}


void sample_discrete(double *prob, int length, double *samples, 
		     int num_samples) {
    int n=0,b=0;
    double urand = 0.0;
    double *cumsum = NULL;
    
    cumsum = malloc(sizeof(double)*length);
    cumsum[0] = prob[0];
    for(b=1; b < length; b++) {
        cumsum[b] = cumsum[b-1] + prob[b];
    }
    
    for(n=0; n < num_samples; n++) {
        urand = ((double)rand()) / RAND_MAX;
        for(b=0; b < length-1 && urand >= cumsum[b]; b++);
        samples[n] = (double)(b+1);
    }
    
    free(cumsum);
}
