#include <cmath>
#include <cstring>
#include <algorithm>
#include "pthread.h"
#include "math.h"
#include "mex.h"
#include "libs.h"
#include <time.h>
#include <cmath>
#include <stdlib.h>   
#define MAX_SUPPORTED_THREADS 100


//xx, preDerivative, derivative, niter
void doIter(double *HH, double *hh, int n, int k, double tolorance, int maxIter, int verbose,
        double *xx, double *iter, int max_threads)
{
    double **H = convert(HH, k, k);
    double **h = convert(hh, n, k);
    double **x = convert(xx, n, k);
    double *sqrtH = createV(k);
    double *iSqrtH = createV(k);
    double e=0, pre=0, der=0;
    
    *iter = 0;
    For(i, k) {
        sqrtH[i] = (H[i][i] > 0 ? sqrt(H[i][i]) : 1.0); 
        iSqrtH[i] = 1.0/sqrtH[i];
    }
    
    For(i, k) For(j, k) H[i][j] *= iSqrtH[i]*iSqrtH[j];
    
    double *its = createV(max_threads);
    int start = 0, end;
    pthread_t threads[MAX_SUPPORTED_THREADS];
    thread_data* args[MAX_SUPPORTED_THREADS];
    
    For(t, max_threads){
        end = start + (n - start) / (max_threads - t);
        args[t] = new thread_data(H, h, k, tolorance, maxIter, verbose, x, &its[t], start, end, sqrtH, iSqrtH);
        if (pthread_create(&threads[t], NULL, thread_iter, (void*)args[t]))
            mexErrMsgTxt("problem with return code from pthread_create()");
        start = end;
    }
    
    For(t, max_threads){
        pthread_join(threads[t], NULL);
        *iter += its[t];
    }
    
    free(its);
    free(H);
    free(h);
    free(x);
    free(sqrtH);
    free(iSqrtH);
    For(t, max_threads) delete args[t];
}

///[xx, preDerivative, derivative, iter] = doIter (HH, hh, xx, maxIter, tolorance, verbose)
void mexFunction(int noOuts, mxArray *outs[], int noIns, const mxArray *inps[])
{
	double *xx, *HH, *hh, *values;
	double tolorance, *error, *preDerivative, *derivative, *iter;
	int maxIter, k, n, verbose, a, b;
    int max_threads;

    //mexPrintf("%d, %d\n", noOuts, noIns);
	// Input Arguments
	getArray(inps[0], &HH, &k, &k);
    getArray(inps[1], &hh, &a, &n);
    if (a != k){
        usage();
		mexPrintf("k=%d != k=%d\n", k, k);
    }
    getArray(inps[2], &xx, &a, &b);
    if (a != k || b != n){
        usage();
		mexPrintf("(%d, %d) != (%d, %d)\n", a, b, k, n);
    }

	maxIter = (int)getDouble(inps[3]);
    tolorance = (double)getDouble(inps[4]);
    verbose = (int)getDouble(inps[5]);
    max_threads = (int)getDouble(inps[6]);
    
    if (verbose){
        mexPrintf("number threads = %d\n", max_threads);
        printf("%d, %d\n", k, k);
        print(&HH[k], 10);
        printf("%d, %d\n", n, k);
        print(&hh[k], 10);
        printf("%d, %d\n", n, k);
        print(&xx[k], 10);
        printf("%d,  %d,  %d,  %.6E, %d\n", n, k, maxIter, tolorance, verbose);
    }

	/// Output arguments
	outs[0] = mxCreateDoubleMatrix(k, n, mxREAL);
	double* newX =  mxGetPr(outs[0]);
    
    outs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    iter  = mxGetPr(outs[1]);
    *iter = 0.0;
    
    //For(i, n*k) newX[i] = xx[i];
    doIter(HH, hh, n, k, tolorance, maxIter, verbose, xx, iter, max_threads);
    *iter /= (n + 0.0);
    For(i, n*k) newX[i] = xx[i];
}
