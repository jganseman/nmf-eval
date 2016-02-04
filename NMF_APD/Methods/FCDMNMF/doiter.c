#include "math.h"
#include "mex.h" 
#include <time.h>
#include <stdlib.h>

void goWiter(double *GW,double *HH, double *W, double* iter, double tol, int n, int k, int maxinner)
{
	// initial maximum function value decreasing over all coordinates. 
	double init=0; 

	// Diagonal of HH
	double *HH_d = (double *)malloc(sizeof(double)*k);
	for ( int i=0 ; i<k ; i++ )
		HH_d[i] = HH[i*k+i];

	// Create SWt : store step size for each variables 
	double *SWt = (double *)malloc(sizeof(double)*k);

	// Get init value 
	for ( int i=0, nowidx=0 ; i<n ; i++ )
	{
		for ( int j=0 ; j<k ; j++, nowidx++ )
		{
			double s = GW[nowidx]/HH_d[j];
			s = W[nowidx]-s;
			if ( s< 0)
				s=0;
			s = s-W[nowidx];
			double diffobj = (-1)*s*GW[nowidx]-0.5*HH_d[j]*s*s;
			if ( diffobj > init )
				init = diffobj;
		}
	}

	// stopping condition
    *iter = 0;
	// coordinate descent 
	for ( int p=0 ; p<n ; p++)
	{
		double *GWp = &(GW[p*k]);
		double *Wp = &(W[p*k]);
        int counter = 0;
		for ( int winner = 0 ; winner < maxinner ; winner++)
		{
            //(*iter)++;
            counter++;
			// find the best coordinate 
			int q = -1;
			double bestvalue = 0;

			for ( int i=0; i<k ; i++ )
			{
				double ss = GWp[i]/HH_d[i];
				ss = Wp[i]-ss;
				if (ss < 0)
					ss=0;
				ss = ss-Wp[i];
				SWt[i] = ss;
				double diffobj = (-1)*(ss*GWp[i]+0.5*HH_d[i]*ss*ss);
				if ( diffobj > bestvalue ) 
				{
					bestvalue = diffobj;
					q = i;
				}
			}
			if ( q==-1 )
				break;

			Wp[q] += SWt[q];
			int base = q*k;
			for ( int i=0 ; i<k ; i++ )
				GWp[i] += SWt[q]*HH[base+i];
			if ( bestvalue < init*tol)
				break;
		}
        *iter += (counter + k - 1) / k;
        //printf("counter = %d, %f\n", counter, tol);
	}
    *iter /= n;
	free(HH_d);
	free(SWt);
}


void usage()
{
	printf("Error calling doiter.\n");
	printf("Usage: Wnew = doiter(GW, HH^T, W, tol, maxinner)\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *GW, *HH, *W, *values;
	double tol;
	int maxinner, k, n;

	// Input Arguments

	GW = mxGetPr(prhs[0]);
	k = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);

	HH = mxGetPr(prhs[1]);
	if ( (mxGetM(prhs[1]) != k) || (mxGetN(prhs[1])!=k) ) {
		usage();
		printf("Error: %d %d HH^T should be a %d by %d matrix\n", mxGetM(prhs[2]), mxGetN(prhs[2]), k, k);
	}

	W = mxGetPr(prhs[2]);
	if ( (mxGetM(prhs[2])!=k) || (mxGetN(prhs[2])!=n) ) {
		usage();
		printf("Error: W should be a %d by %d matrix", k, n);
	}

	values = mxGetPr(prhs[3]);
	tol = values[0];

	values = mxGetPr(prhs[4]);
	maxinner = values[0];

	/// Output arguments
	plhs[0] = mxCreateDoubleMatrix(k, n, mxREAL);
	double *Wout = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* iter  = mxGetPr(plhs[1]);
    *iter = 0;
		
	goWiter(GW, HH, W, iter, tol, n, k, maxinner);

	for ( int i=0 ; i<n*k ; i++ )
		Wout[i] = W[i];
}
