#ifndef __LIBS__
//#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include "pthread.h"
#include "math.h"
#include "mex.h"
#include <time.h>
#include <cmath>
#include <stdlib.h>   
#define MAX_SUPPORTED_THREADS 100

#define For(i,n) for(int i=0; i<n; i++)
#define max(a, b) (a > b ?  a : b)
#define abs(x) (x >= 0 ? x : -x)
#define eps 1e-30

using namespace std;

//inner product of two vectors
double dotVVFrom(double *a, double *b, int from, int n){
    double res = 0; 
    double *ai, *bi, *end=a+n;
    for (ai=(a+from), bi=(b+from); ai!=end; ai++,bi++) res += (*ai)*(*bi);
    return res;
}
void assign(double *x, double *y, int k){
    For(i, k) x[i] = y[i];
}
//inner product of two vectors
double dotVV(double *a, double *b, int n){
    return dotVVFrom(a, b, 0, n);
}

//sum of two vectors
void addVV(double *a, double *b, int n, double *res){
    For(i,n) res[i] = a[i] + b[i];
}

void addkV(double alpha, double *a, int n, double *res){
    For(i, n) res[i] += alpha*a[i];
}

//subtract of two vectors
void subVV(double *a, double *b, int n, double *res){
    For(i,n) res[i] = a[i] - b[i];
}

//a product of matrix and vector
void dotMV(double **a, double *b, int m, int n, double *res){
    For(i,m) res[i] = dotVV(a[i], b, n);
}

//a product of matrix and vector
void dotVM(double **a, double *b, int m, int n, double *res){
    For(i,m) res[i] = dotVV(a[i], b, n);
}



//a product of v^TQv
double dotVMV(double *x, double **M, int k){
    double res = 0, sub;
    for (int i = 0; i < k; i++) {
        res += x[i] * x[i] * M[i][i];
        sub = 0;
        for (int j = i + 1; j < k; j++)
            sub += M[i][j] * x[j];
        res += 2 * x[i] * sub;
    }
    return res;
}

//Remove negative elements of a vector
void removeVVNegative(double* a, int n, double *b){
    For(i, n) b[i] = (a[i] >= 0 ? a[i] : 0);
}

//Remove KKT-satisfied elements of a vector
void removeKTTElements(double* df, double* x, int n, double *dbarf){
    For(i, n) dbarf[i] = (df[i] < 0 || x[i] > 0) ? df[i] : 0;
}


//fro norm ||a||^2_2 
double square(double *a, int n){
    double res = 0;
    For(i, n) res += a[i]*a[i];
    return res;
}

//get index of max element
int maxV(double *a, int n){
    int k = 0;
    For(i, n) k = (a[i] > a[k] ? i : k);
    return k;
}

//convert double* into double**
double **convert(double *a, int m, int n){
    double **res = (double **)malloc(sizeof(double*)*m);
    For(i,m) res[i] = &a[i*n];
    return res;
}

//create a vector
double *createV(int n){
    double *res = (double *)malloc(sizeof(double)*n);
    For(i,n) res[i] = 0;
    return res;
}

double getMaxReduce(double **H, double *h, int k, double* df, double* x){
    double dx, maxValue=-1e30, dfx;
    int pos = -1, v;
    For(i, k)
        if ((df[i] < 0 || x[i] > 0) && H[i][i] > 0){
            dx = max(0.0, x[i] - df[i]) - x[i];
            dfx = (-0.5 * dx - df[i]) * dx;
            if (maxValue < dfx){
                pos = i;
                maxValue = dfx;
            }
        }
    return maxValue;
}

class thread_data {
public:
    double **HH;
    double **hh;
    int k;
    double tolorance;
    int maxIter;
    int verbose;
    double **xx;
    double *iter;
    int start;
    int end;
    double *sqrtH;
    double *iSqrtH;
    
    double *df;
    double *dbarf;
    double *newX;
    double *saveDf;
    double *saveX;
    double *dx;
    
    thread_data(){
        df = createV(k);
        dbarf = createV(k);
        newX = createV(k);

        saveDf = createV(k);
        saveX = createV(k);
        dx = createV(k);
    }
    
    thread_data(double **HH, double **hh, int k, 
            double tolorance, int maxIter, int verbose,
            double **xx, double *iter, int start, int end, double *sqrtH, double *iSqrtH){
        this->HH = HH;
        this->hh = hh;
        this->k = k;
        this->tolorance = tolorance;
        this->maxIter = maxIter;
        this->verbose = verbose;
        this->xx = xx;
        this->iter = iter;
        this->start = start;
        this->end = end;
        this->sqrtH = sqrtH;
        this->iSqrtH = iSqrtH;
        
        df = createV(k);
        dbarf = createV(k);
        newX = createV(k);

        saveDf = createV(k);
        saveX = createV(k);
        dx = createV(k);
    }
    
    ~thread_data(){
        free(df);
        free(dbarf);
        free(newX);

        free(saveDf);
        free(saveX);
        free(dx);
    }
};

int FCD(double **H, double *h, int k, double* df, double* x){
    double dx, maxValue = 0, dfx, first = 1e100;
    int pos = -1;
    For(iter, k){
        pos = -1;
        maxValue = 0;
        For(i, k) 
          if ((df[i] < 0 || x[i] > 0) && H[i][i] > eps) {
              dx = max(0.0, x[i] - df[i]) - x[i];
              dfx = (- 0.5 * H[i][i] * dx - df[i]) * dx;
              if (maxValue < dfx){
                  pos = i;
                  maxValue = dfx;
              }
          }
        if (iter == 1) first = maxValue;
        //if (maxValue < first * 1e-3) return iter;
        //mexPrintf("v = %d, max = %.5e\n", pos, maxValue);
        if (pos == -1) return iter;
        //pos = iter;
        dx = max(0.0, x[pos] - df[pos]) - x[pos];
        x[pos] += dx;
        double *h = H[pos];
        For(i, k) df[i] += h[i] * dx;
    }
    return k;
}

void solveNQP(double** H, double *h, int k, double tolorance, int maxIter, int verbose,
        double *x, double* iter, double *stopMax, double* sqrtH, double* iSqrtH, thread_data* args){
    double *df = args->df;
    double *dbarf = args->dbarf;
    double *newX = args->newX;
    
    double *saveDf = args->saveDf;
    double *saveX = args->saveX;
    double *dx = args->dx;
    For(i, k) df[i] = dbarf[i] = newX[i] = saveDf[i] = saveX[i] = dx[i] = 0;
    
    For(i, k) x[i] *= sqrtH[i], h[i] *= iSqrtH[i];
    double den, a;
    double error, preDerivative, derivative;
    dotMV(H, x, k, k, df);
    addVV(df, h, k, df);
    removeKTTElements(df, x, k, dbarf);
    preDerivative = square(dbarf, k);
    *iter = 0;
    //double dfInit = getMaxReduce(H, h, k, df, x);
    while (*iter < maxIter){
        *iter = *iter + 1;
        assign(saveDf, df, k);
        assign(saveX, x, k);
        
        //antilopsided
        den = dotVMV(dbarf, H, k);
        a = dotVV(dbarf, dbarf, k) / den;
        if (a != a || abs(a) < eps || abs(a) > 1e30) break;
        For(i, k) newX[i] = max(0.0, x[i] - a * dbarf[i]) - x[i];
        For(i, k) if (newX[i] != 0) addkV(newX[i], H[i], k, df);
        addVV(x, newX, k, x);
        //derivative = square(dbarf, k);
        
        //Greedy coordinate descent
        int niter = FCD(H, h, k, df, x);
        removeKTTElements(df, x, k, dbarf);
        derivative = square(dbarf, k);
        if (derivative <= tolorance * preDerivative || derivative <= *stopMax)
            break;
        
        //Accelerated step
        For(i, k) dx[i] = saveX[i] - x[i];
        den = dotVMV(dx, H, k);
        a = dotVV(saveDf, dx, k) / den;
        if (a != a || abs(a) < eps || abs(a) > 1e30) break;
        For(i, k) dx[i] = max(0.0, saveX[i] - a * dx[i]) - saveX[i];
        assign(df, saveDf, k);
        assign(x, saveX, k);
        For(i, k) if (dx[i] != 0) addkV(dx[i], H[i], k, df);
        addVV(x, dx, k, x);
        //removeKTTElements(df, x, k, dbarf);
        //derivative = square(dbarf, k);
        
        //Greedy coordinate descent
        FCD(H, h, k, df, x);
        removeKTTElements(df, x, k, dbarf);
        derivative = square(dbarf, k);
        if (derivative <= tolorance * preDerivative || derivative <= *stopMax)
            break;
    }
    if (*stopMax < derivative) *stopMax = derivative;
    For(i, k) x[i] *= iSqrtH[i];// h[i] /= sqrtH[i];
    if (verbose) mexPrintf("out antilopsided\n");
}


void antilopsided(double** H, double *h, int k, double tolorance, int maxIter, int verbose,
        double *x, double* iter, double *stopMax, double* sqrtH){
    double *df = createV(k);
    double *dbarf = createV(k);
    double *newX = createV(k);
    For(i, k) x[i] *= sqrtH[i], h[i] /= sqrtH[i];
    double den, a;
    double error, preDerivative, derivative;
    dotMV(H, x, k, k, df);
    addVV(df, h, k, df);
    removeKTTElements(df, x, k, dbarf);
    preDerivative = square(dbarf, k);
    *iter = 0;
    //double dfInit = getMaxReduce(H, h, k, df, x);
    while (*iter < maxIter){
        *iter = *iter + 1;
        
        //antilopsided
        den = dotVMV(dbarf, H, k);
        a = dotVV(dbarf, dbarf, k) / den;
        if (a != a || abs(a) < eps || abs(a) > 1e30) break;
        For(i, k) newX[i] = max(0.0, x[i] - a * dbarf[i]) - x[i];
        For(i, k) if (newX[i] != 0) addkV(newX[i], H[i], k, df);
        addVV(x, newX, k, x);
        removeKTTElements(df, x, k, dbarf);
        derivative = square(dbarf, k);
        
        //fast coordinate descent
        FCD(H, h, k, df, x);
        removeKTTElements(df, x, k, dbarf);
        derivative = square(dbarf, k);
        
        //double dfCurrent = getMaxReduce(H, h, k, df, x);
        //if (dfCurrent < tolorance * dfInit)
        //if (derivative <= tolorance * preDerivative || derivative <= *stopMax)
        break;
    }
    if (*stopMax < derivative) *stopMax = derivative;
    For(i, k) x[i] /= sqrtH[i];// h[i] /= sqrtH[i];
    free(df);
    free(dbarf);
    free(newX);
    if (verbose) mexPrintf("out antilopsided\n");
}


void* thread_iter(void* thread_args){
    thread_data *args = (thread_data *) thread_args;
    double stopMax = 0;
    for (int i=args->start; i < args->end; i++) {
        double it = 0;
        //antilopsided(args->HH, args->hh[i], args->k, args->tolorance, 
        //        args->maxIter, args->verbose, args->xx[i], &it, &stopMax, args->sqrtH);
        solveNQP
        (args->HH, args->hh[i], args->k, args->tolorance, args->maxIter, args->verbose, args->xx[i], &it, &stopMax, args->sqrtH, args->iSqrtH, args);
        *(args->iter) += it;
    }
    return NULL;
}


void print(double* a, int k){
    For(i,k) mexPrintf("%.4f ", a[i]);
    mexPrintf("\n");
}

void usage()
{
	printf("Error calling doiter.\n");
	printf("Usage: Wnew = doiter(GW, HH^T, W, tol, maxinner)\n");
}

void getArray(const mxArray* a, double **res, int *m, int *n){
    *res = mxGetPr(a);
    *m = mxGetM(a);
    *n = mxGetN(a);
}

double getDouble(const mxArray* a){
    double *values = mxGetPr(a);
    return values[0];
}

#endif
