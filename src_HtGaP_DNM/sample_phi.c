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
    plhs[0] = mxCreateDoubleMatrix(N, K, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, K, mxREAL);
    phi_lambda_ct = mxGetPr(plhs[0]);
    sum_phi_lambda_ct = mxGetPr(plhs[1]);
    
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
        
}

