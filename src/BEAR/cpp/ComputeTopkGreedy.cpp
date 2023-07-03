

#include "mex.h"
#include "matrix.h"
#include "string.h"

#include <algorithm>

#include <fstream>

#define FUNC_NAME					"ComputeTopkGreedy"
#define ERRMSG_WRONG_INPUT_NUM		FUNC_NAME" requires exactly 2 inputs."
#define ERRMSG_WRONG_OUTPUT_NUM		FUNC_NAME" produces at most 1 output."
#define ERRMSG_WRONG_INPUT_TYPE		FUNC_NAME" requires a sparse matrix with the same row and column size as the input."
#define ERRMSG_TOO_LARGE_K			"In " FUNC_NAME ", the graph size should be at least the same as k."

//#define TEST

class ComputeTopkGreedy{
public:
	static void run(double* pr, mwIndex* ir, mwIndex* jc, mwSize n, mwSize nz, double* topk, mwSize k){
		ComputeTopkGreedy::pr = pr;
		ComputeTopkGreedy::ir = ir;
		ComputeTopkGreedy::jc = jc;
		ComputeTopkGreedy::n = n;
		ComputeTopkGreedy::nz = nz;
		ComputeTopkGreedy::topk = topk;
		ComputeTopkGreedy::k = k;

		rho = new mwSize[n];
		sorted = new mwSize[n];
		cumul = new mwSize[n];
		
		run();

		delete [] rho;
		delete [] sorted;
		delete [] cumul;
	}
private:
	static void run(){
		init_rho_cumul();
		sort_rho();
		for(mwSize i=0; i<k; i++){
			double val = double(rho[sorted[i]] - cumul[sorted[i]]);
			if( i != n && val < rho[sorted[i+1]] ){
				update_rho_cumul();
				sort_rho();
			}

#ifdef TEST
			degsize[i] = double(rho[sorted[i]] - cumul[sorted[i]]);
#endif

			topk[i] = double(sorted[i]);
			rho[sorted[i]] = 3*n;

			cumulate(sorted[i]);
		}

		for(mwSize i=0; i<k; i++){
			topk[i]++;
		}
	}
	static void init_rho_cumul(){
		memset(rho, 0, n*sizeof(*rho));
		for(mwSize c=0; c<n; c++){
			mwSize lr, rr;
			get_row_range(c,lr,rr);
			rho[c] = rr - lr;

			cumul[c] = 0;
		}
	}
	static void update_rho_cumul(){
		for(int i=0; i<n; i++){
			rho[i] -= cumul[i];
			cumul[i] = 0;
		}
	}
	static void sort_rho(){
		for(int i=0; i<n; i++){
			sorted[i] = i;
		}
		qsort(sorted, n, sizeof(*sorted), ComputeTopkGreedy::compare_descend);
	}
	static int compare_descend(const void* a, const void* b){
		return rho[*(mwSize*)b] - rho[*(mwSize*)a];
	}
	static void cumulate(mwSize id){
		mwSize lr, rr;
		get_row_range(id, lr, rr);
		for(mwSize i=lr; i<rr; i++){
			cumul[ir[i]]++; 
		}
	}
	static void get_row_range(mwSize c, mwSize& lr, mwSize& rr){
		lr = jc[c]; rr = jc[c+1];
	}
private:
	// input
	static double* pr;
	static mwIndex* ir;
	static mwIndex* jc;
	static mwSize n;
	static mwSize nz;
	static mwSize k;

	// output
	static double* topk;

	// internal
	static mwSize* rho;
	static mwSize* sorted;
	static mwSize* cumul;

#ifdef TEST
	public: 
		static double* degsize;
#endif
};
double* ComputeTopkGreedy::pr = NULL;
mwIndex* ComputeTopkGreedy::ir = NULL;
mwIndex* ComputeTopkGreedy::jc = NULL;
mwSize ComputeTopkGreedy::n = -1;
mwSize ComputeTopkGreedy::nz = -1;
mwSize ComputeTopkGreedy::k = -1;
double* ComputeTopkGreedy::topk = NULL;
mwSize* ComputeTopkGreedy::rho = NULL;
mwSize* ComputeTopkGreedy::sorted = NULL;
mwSize* ComputeTopkGreedy::cumul = NULL;

#ifdef TEST
double* ComputeTopkGreedy::degsize = NULL;
#endif



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if(nrhs != 2){
		mexErrMsgTxt(ERRMSG_WRONG_INPUT_NUM);
	}
	if(nlhs > 2){
		mexErrMsgTxt(ERRMSG_WRONG_OUTPUT_NUM);
	}
	if( !mxIsSparse(prhs[0]) ){
		mexErrMsgTxt(ERRMSG_WRONG_INPUT_TYPE);
	}
	if( mxGetM(prhs[0]) != mxGetN(prhs[0]) ){
		mexErrMsgTxt(ERRMSG_WRONG_INPUT_TYPE);
	}

	// Prepare variables for the input
	double* pr = mxGetPr(prhs[0]);
	mwIndex* ir = mxGetIr(prhs[0]);
	mwIndex* jc = mxGetJc(prhs[0]);
	mwSize n = mxGetM(prhs[0]);
	mwSize nz = mxGetNzmax(prhs[0]);
	mwSize k = mwSize(mxGetScalar(prhs[1]));

	// Prepare variables for the output
	mxArray* mx_topk = mxCreateDoubleMatrix(k, 1, mxREAL);
	double* top = mxGetPr(mx_topk);

#ifdef TEST
	// Prepare variables for the test output
	mxArray* mx_degsize = mxCreateDoubleMatrix(1,k,mxREAL);
	double* degsize = mxGetPr(mx_degsize);
	ComputeTopkGreedy::degsize = degsize;
#else
	mxArray* mx_degsize = mxCreateDoubleMatrix(0,0,mxREAL);
#endif

	if( n < k ){
		mexErrMsgTxt(ERRMSG_TOO_LARGE_K);
	}

	ComputeTopkGreedy::run(pr, ir, jc, n, nz, top, k);

	if(nlhs > 0) plhs[0] = mx_topk;

#ifdef TEST
	if(nlhs > 1) plhs[1] = mx_degsize;
#endif
}






















