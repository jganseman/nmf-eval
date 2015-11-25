#include <algorithm>
#include <time.h>
#include <cmath>
#include <stdlib.h>  

#include "pthread.h"
#include "math.h"
#include "mex.h"
#include "libs.h" 
#define MAX_SUPPORTED_THREADS 100

double maxGrad = 0;
double maxIterGlobal = 0;

void exact_antilopsided_nqp(double** H, double *h, int k, double tolorance, int maxIter, int verbose,
        double *x, double* iter){
    double *df = createV(k);
    double *dbarf = createV(k);
    double *newX = createV(k);
    
    double den, a;
    double error, preDerivative, derivative;
    For(i, k) x[i] = 0;
    dotMV(H, x, k, k, df);
    addVV(df, h, k, df);
    removeKTTElements(df, x, k, dbarf);
    preDerivative = square(dbarf, k);
    *iter = 0;
    double previous = preDerivative;
    double optimal = preDerivative;
    double *sol = createV(k);
    
    while (*iter < 2*k){
        *iter = *iter + 1;
        //antilopsided
        den = dotVMV(dbarf, H, k);
        a = square(dbarf, k) / den;
        if (a != a) break;
        For(i, k) if (dbarf[i] != 0) {
            newX[i] = max(0.0, x[i] - a * dbarf[i]) - x[i];
            addkV(newX[i], H[i], k, df);
            x[i] = max(0.0, x[i] - a * dbarf[i]);
        }
        //fast coordinate descent
        //FCD(H, h, k, df, x);
        FCD(H, h, k, df, x);
        removeKTTElements(df, x, k, dbarf);
        derivative = square(dbarf, k);
        if (derivative < optimal){
            optimal = derivative;
            assign(sol, x, k);
        }
        //mexPrintf("iter %d = %e\n", (int)(*iter), derivative);
        if ((previous < derivative && *iter > 2 && derivative < 1) || (derivative < 1e-22))
            break;
        previous = derivative;
        //if (*iter == 3) break;
    }
    if (verbose && derivative > maxGrad){
        maxGrad = derivative;
        mexPrintf("iter %d = %e\n", (int)(*iter), derivative);
    }
    if (verbose && *iter > maxIterGlobal){
        maxIterGlobal = *iter;
        mexPrintf("maxIter = %d\n", (int)maxIterGlobal);
    }
    assign(x, sol, k);
    free(sol);
    free(df);
    free(dbarf);
    free(newX);
    //if (verbose) mexPrintf("out antilopsided\n");
}

void* thread_exact_antilopsided_nqp(void* thread_args){
    thread_data *args = (thread_data *) thread_args;
    *(args->iter) = 0;
    for (int i=args->start; i < args->end; i++) {
        double it = 0;
        exact_antilopsided_nqp(args->HH, args->hh[i], args->k, args->tolorance, 
                args->maxIter, args->verbose, args->xx[i], &it);
        *(args->iter) += it;
    }
    return NULL;
}

//xx, preDerivative, derivative, niter
void nqp4nnls(double *HH, double *hh, int n, int k, double tolorance, int maxIter, int verbose,
        double *xx, double *iter, int max_threads)
{
    //struct stableParams params;
    double **H = convert(HH, k, k);
    double **h = convert(hh, n, k);
    double **x = convert(xx, n, k);
    double *sqrtH = createV(k);
    //double e=0, pre=0, der=0;
    
    *iter = 0;
    For(i, k) sqrtH[i] = (H[i][i] > 0 ? sqrt(H[i][i]) : 1); 
    For(i, n) For(j, k){
        x[i][j] *= sqrtH[j];
        h[i][j] /= sqrtH[j];
    }
    For(i, k) For(j, k) H[i][j] /= sqrtH[i]*sqrtH[j];
    
    double *its = createV(MAX_SUPPORTED_THREADS);
    int start = 0, end;
    pthread_t threads[MAX_SUPPORTED_THREADS];
    thread_data args[MAX_SUPPORTED_THREADS];
    if (max_threads > 1 || true) {
        For(t, max_threads){
            end = start + (n - start) / (max_threads - t);
            args[t] = thread_data(H, h, k, tolorance, maxIter, verbose, x, &its[t], start, end);
            if (pthread_create(&threads[t], NULL, thread_exact_antilopsided_nqp, (void*)&args[t]))
                mexErrMsgTxt("problem with return code from pthread_create()");
            start = end;
        }
        For(t, max_threads) pthread_join(threads[t], NULL);
        //For(t, max_threads) pthread_join(threads[t], NULL);
        For(t, max_threads) *iter += its[t];
    } else {
        //mexPrintf("Come 1 thread\n");
        args[0] = thread_data(H, h, k, tolorance, maxIter, verbose, x, &its[0], 0, n);
        thread_exact_antilopsided_nqp((void*)&args[0]);
        *iter += its[0];
    }
    
    For(i, n) For(j, k)
            x[i][j] /= sqrtH[j];
    
    free(its);
    free(H);
    free(h);
    free(x);
    free(sqrtH);
}

///[xx, preDerivative, derivative, iter] = doIter (HH, hh, xx, maxIter, tolorance, verbose)
void mexFunction(int noOuts, mxArray *outs[], int noIns, const mxArray *inps[])
{
    double *HH, *hh, *values;
	double tolorance, *error, *preDerivative, *derivative, *iter;
	int maxIter, k, n, verbose, a, b;
    int max_threads;
    
	getArray(inps[0], &HH, &k, &k);
    getArray(inps[1], &hh, &a, &n);
    
    if (a != k){
        usage();
		mexPrintf("k=%d != k=%d\n", k, k);
    }
    
    maxIter = (int)getDouble(inps[2]);
    verbose = (int)getDouble(inps[3]);
    max_threads = (int)getDouble(inps[4]);
    
    if (verbose){
        mexPrintf("number threads = %d\n", max_threads);
        printf("%d, %d\n", k, k);
        print(&HH[0], k);
        printf("%d, %d\n", n, k);
        print(&hh[0], k);
        printf("%d,  %d,  %d,  %.6E, %d\n", n, k, maxIter, tolorance, verbose);
    }

	/// Output arguments
	outs[0] = mxCreateDoubleMatrix(k, n, mxREAL);
	double* newX =  mxGetPr(outs[0]);
    For(i, n*k) newX[i] = 0.0;
    
    outs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    iter  = mxGetPr(outs[1]);
    *iter = 0;
    tolorance = 1e-30;
    
    nqp4nnls(HH, hh, n, k, tolorance, maxIter, verbose, newX, iter, max_threads);
    
    if (verbose) mexPrintf("Entered Finish\n");
}
