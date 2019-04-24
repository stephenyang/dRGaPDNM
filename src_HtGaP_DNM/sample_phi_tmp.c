
#include "mex.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>

/*
 * #include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
 * phi_ct = sample_phi_tmp(phi_ct, phi_pt, phi_ft, lambda_ct, Y_n_k_ft, m_n_k_ct, c_n, P_n_k_ft); */

/*
double rgam(const double a, const double b)
{
  
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);
    mxArray *rsize,*rNorm; double *rsizeP;
    double *randnorm;
    while (1)
      {
        do
          {
            rsize = mxCreateDoubleMatrix(1,2,mxREAL);
            rsizeP= mxGetPr(rsize); rsizeP[0] = 1;rsizeP[1] = 1;
            mexCallMATLAB(1, &rNorm, 1, &rsize, "randn");  
            randnorm = mxGetPr(rNorm);
            x = randnorm[0];
            v = 1.0 + c * x;
          }
        while (v <= 0);

        v = v * v * v;
        u = rand();

        if (u < 1 - 0.0331 * x * x * x * x) 
          break;

        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }
    
    return b * d * v;
}
*/

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /*phi_ct, phi_pt, phi_ft, lambda_ct, Y_n_k_ft, m_n_k_ct, c_n, P_n_k_ft*/
    mwSize N, K;
    double *phi_ct, *phi_pt, *phi_ft, *lambda_ct, *Y_n_k_ft, *m_n_k_ct, *c_n, *P_n_k_ft;
    double *phi_lambda_ct, *sum_phi_lambda_ct;
    phi_ct = mxGetData(prhs[0]);
    phi_pt = mxGetData(prhs[1]);
    phi_ft = mxGetData(prhs[2]);
    lambda_ct = mxGetData(prhs[3]);
    Y_n_k_ft = mxGetData(prhs[4]);
    m_n_k_ct = mxGetData(prhs[5]); 
    c_n = mxGetData(prhs[6]); 
    P_n_k_ft = mxGetData(prhs[7]);
    
    N = mxGetM(prhs[0]);
    K = mxGetN(prhs[0]);
    
    phi_lambda_ct = mxCreateDoubleMatrix(N, K, mxREAL);
    sum_phi_lambda_ct = mxCreateDoubleMatrix(1, K, mxREAL);
    
    phi_ct = mxGetPr(plhs[0]);

    mwIndex i, k, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            phi_lambda_ct[i + j * N] = 0;
            for (k = 0; k < K; k++)
            {
                phi_lambda_ct[i + j * N] = phi_lambda_ct[i + j * N] + phi_ct[i + k * N]*lambda_ct[k + j * K];
            }
        }
    }
    
    for (j = 0; j < K; j++)
    {
        sum_phi_lambda_ct[j] = 0;
        for (i = 0; i < N; i++)
        {
            sum_phi_lambda_ct[j] = sum_phi_lambda_ct[j] + phi_lambda_ct[i + j * N];
        }
    }
    
    /* phi*/
    /*
    double delta;
    for (i = 0; i < N; i++)
    {
        for (k = 0; k < K; k++)
        {sum_phi_lambda_ct[k] = sum_phi_lambda_ct[k] - phi_lambda_ct[i + k * N];}
        
        for (k = 0; k < K; k++)
        {
            if (1 - P_n_k_ft[i + k*N] > -1000)
                delta = 1 - P_n_k_ft[i + k*N];
            else
                delta = -1000;
            double gam_a = phi_pt[i + k*N] + m_n_k_ct[i + k*N] + Y_n_k_ft[i+k*N];
            double gam_b = 1./(c_n[i] + sum_phi_lambda_ct[k] - log(delta));
            phi_ct[i + k*N] = 0;
          }
        
        for (k = 0; k < K; k++)
        {
            sum_phi_lambda_ct[k] = sum_phi_lambda_ct[k] + phi_ct[i + k*N] * lambda_ct[k+k*K];
        }
     
        
    }
    */
    
}

