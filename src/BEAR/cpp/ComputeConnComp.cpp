

#include "mex.h"
#include "matrix.h"

#include <queue>

#define FUNC_NAME "ComputeConnComp"
#define ERRMSG_WRONG_INPUT_NUM FUNC_NAME" requires exactly 2 inputs."
#define ERRMSG_WRONG_OUTPUT_NUM FUNC_NAME" produces at most 3 outputs."
#define ERRMSG_WRONG_INPUT_TYPE FUNC_NAME" requires a sparse matrix with the same row and column size as the input."
#define ERRMSG_TOO_MANY_ITERATIONS(method,n) "In "#method" of "FUNC_NAME", the number of iterations exceeds "#n"."


class ConnComp{
public:
	static void run(double* pr, mwIndex* ir, mwIndex* jc, mwSize n, mwSize nz, double* top, mwSize topnum, 
			mwSize* ccnum, double* cclabel, double* cchist, double* ccatt){
		ConnComp::pr = pr;
		ConnComp::ir = ir;
		ConnComp::jc = jc;
		ConnComp::n = n;
		ConnComp::nz = nz;
		ConnComp::top = top;
		ConnComp::topnum = topnum;
		ConnComp::ccnum = ccnum;
		ConnComp::cclabel = cclabel;
		ConnComp::cchist = cchist;
		ConnComp::ccatt = ccatt;

		memset(cclabel, 0, n*sizeof(*cclabel));

		topmark = new mwSize[n];
		memset(topmark, 0, n*sizeof(*topmark));
		for(mwSize i=0; i<topnum; i++){
			topmark[mwSize(top[i])-1] = i + 1;
		}

		run();

		delete [] topmark;
	}
private:
	static void run(){
		double label = 1;
		for(mwSize i=0; i<topnum; i++){
			mwSize c = mwSize(top[i])-1;
			mwSize lr, rr;
			get_row_range(c, lr, rr);

			for(mwSize j=lr; j<rr; j++){
				mwSize u = ir[j];
				if( cclabel[u] == 0 ){
					mwSize num = bfs(u, label);
					if(num > 0){
						cchist[mwSize(label)-1] = double(num);
						label++;
					}
					else{
						continue;
					}
				}
				ccatt[int(cclabel[u])-1] = double(i+1);				
			}
		}
		for(mwSize i=0; i<topnum; i++,label++){
			cclabel[int(top[i])-1] = label;
			cchist[mwSize(label)-1] = 1;
			ccatt[mwSize(label)-1] = label;
		}
		for(mwSize i=0; i<n; i++){
			if( cclabel[i] == 0 ){
				mwSize num = bfs(i, label);
				cchist[mwSize(label)-1] = double(num);
				label++;
				ccatt[int(cclabel[i])] = 0;

				//cclabel[i] = label;
				//cchist[mwSize(label)-1] = 1;
				//label++;
			}
		}
		*ccnum = mwSize(label - 1);
	}

	static mwSize bfs(mwSize id, double label){
		if( topmark[id] > 0 ){
			return 0;
		}

		int num = 1, niter = 0;
		std::queue<mwSize> que;
		que.push(id);
		cclabel[id] = label;
		while( !que.empty() ){
			mwSize u = que.front();
			que.pop();
		
			mwSize lr, rr;
			get_row_range(u, lr, rr);
			for(mwSize i=lr; i<rr; i++){
				if(topmark[ir[i]] == 0 && cclabel[ir[i]] == 0){
					que.push(ir[i]);
					cclabel[ir[i]] = label;
					num++;
				}
			}

			if( niter++ > nz*3 ){
				mexErrMsgTxt(ERRMSG_TOO_MANY_ITERATIONS(bfs, nz));
			}
		}

		return num;
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
	static double* top;
	static mwSize topnum; 

	// output
	static mwSize* ccnum; 
	static double* cclabel; 
	static double* cchist;
	static double* ccatt;

	// internal (should be deleted)
	static mwSize* topmark;
};

double* ConnComp::pr = NULL;
mwIndex* ConnComp::ir = NULL;
mwIndex* ConnComp::jc = NULL;
mwSize ConnComp::n = -1;
mwSize ConnComp::nz = -1;
double* ConnComp::top = NULL;
mwSize ConnComp::topnum = NULL;
mwSize* ConnComp::ccnum = 0;
double* ConnComp::cclabel = NULL;
double* ConnComp::cchist = NULL;
double* ConnComp::ccatt = NULL;
mwSize* ConnComp::topmark = NULL;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	if(nrhs != 2){
		mexErrMsgTxt(ERRMSG_WRONG_INPUT_NUM);
	}
	if(nlhs > 4){
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

	double* top = mxGetPr(prhs[1]);
	mwSize topnum = mxGetM(prhs[1]) * mxGetN(prhs[1]);

	// Prepare variables for the outputs
	mxArray* mx_ccnum = mxCreateDoubleScalar(-1);
	mxArray* mx_cclabel = mxCreateDoubleMatrix(n, 1, mxREAL);
	mxArray* mx_cchist = mxCreateDoubleMatrix(n, 1, mxREAL); 
	mxArray* mx_ccatt = mxCreateDoubleMatrix(n, 1, mxREAL); 
	
	mwSize ccnum = -1;
	double* cclabel = mxGetPr(mx_cclabel);	memset(cclabel, 0, n*sizeof(double));
	double* cchist = mxGetPr(mx_cchist);	memset(cchist, 0, n*sizeof(double));
	double* ccatt = mxGetPr(mx_ccatt);		memset(ccatt, 0, n*sizeof(double));

	try{
		ConnComp::run(
			pr, ir, jc, n, nz, top, topnum,		// input
			&ccnum, cclabel, cchist, ccatt);	// output	
	}catch (std::exception& e){
		mexErrMsgTxt(e.what());
	}

	mxSetM(mx_cchist, ccnum);
	mxSetM(mx_ccatt, ccnum);

	*mxGetPr(mx_ccnum) = double(ccnum);
	if(nlhs > 0) plhs[0] = mx_ccnum;
	if(nlhs > 1) plhs[1] = mx_cclabel;
	if(nlhs > 2) plhs[2] = mx_cchist;
	if(nlhs > 3) plhs[3] = mx_ccatt;
}

