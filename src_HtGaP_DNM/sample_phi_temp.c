#include "mex.h"
#include "string.h"
#include <math.h>
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    /*phi_ct, phi_pt, phi_ft, lambda_ct, Y_n_k_ft, m_n_k_ct, c_n, P_n_k_ft*/
    mwSize N, K;
    double *phi_ct, *phi_pt, *phi_ft, *lambda_ct, *Y_n_k_ft, *m_n_k_ct, *c_n, *P_n_k_ft;
    double *phi_lambda_ct, *sum_phi_lambda_ct,*phi_ct_update;
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
    
    plhs[0] = mxCreateDoubleMatrix(N, K, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, K, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N, K, mxREAL);
    
    phi_lambda_ct = mxGetPr(plhs[0]);
    sum_phi_lambda_ct = mxGetPr(plhs[1]);
    phi_ct_update = mxGetPr(plhs[2]);
    
    mwIndex i, k, j;
    /*
    mxArray *gam_params, *gam_rnd;double *gam_paramsP;
    double *randnorm;
    */
    mxArray *rsize,*rNorm; double *rsizeP;
    double *randnorm;
    
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
            
            
            /*
             * gam_paramsP = mxGetPr(gam_params);
             * gam_paramsP[0] = 5;
             * gam_paramsP[1] = 2;
             * mexCallMATLAB(1, &gam_rnd, 1, &gam_params, "gamrnd");
             * randnorm = mxGetPr(gam_rnd);
             */
            rsize = mxCreateDoubleMatrix(1,2,mxREAL);
            rsizeP= mxGetPr(rsize); 
            rsizeP[0] = phi_pt[i + k*N] + m_n_k_ct[i + k*N] + Y_n_k_ft[i+k*N];
            rsizeP[1] = 1;
            mexCallMATLAB(1, &rNorm, 1, &rsize, "randg");  
            randnorm = mxGetPr(rNorm);
            phi_ct_update[i + k*N]  = randnorm[0];
            phi_ct_update[i + k*N]  = phi_ct_update[i + k*N]/(c_n[i] + sum_phi_lambda_ct[k] - log(delta));
          }
        
        for (k = 0; k < K; k++)
        {
            sum_phi_lambda_ct[k] = sum_phi_lambda_ct[k] + phi_ct_update[i + k*N] * lambda_ct[k+k*K];
        }
     
        
    }
    
        
}