#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <queue>
#include <algorithm>
#include "mex.h"
#include "matrix.h"
#include "lapacke.h"
#include "omp.h"
#include "polymesh_FVM.h"
#include "reconstruction.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [up, um] = vandermonde_from_stencil_order3(vertices, edges, cells, stencil_array)
{
	VerticesFVM vertices(prhs[0]);
	EdgesFVM edges(prhs[1]);
	CellsFVM cells(prhs[2]);
	uint32_t *stencil_array = mxGetUint32s(prhs[3]);
	if (stencil_array == nullptr) {
		mexErrMsgIdAndTxt("MEX:FVM_error", "Can't load input stencil array: type must be uint32");
	}
	size_t m = mxGetM(prhs[3]);
	size_t n = mxGetN(prhs[3]);
	if (n > 1) {
		mexErrMsgIdAndTxt("MEX:FVM_error", "Input stencil must be a column vector");
	}
	if (m < 6) {
		mexErrMsgIdAndTxt("MEX:FVM_error", "Input stencil is too small (less than 6 elements)");
	}

	uint32_t i_center = stencil_array[0];
	double x0 = cells.cx[i_center-1];
	double y0 = cells.cy[i_center-1];
	double h = cells.h[i_center-1];

	std::vector<double> V(6*m);
	for (size_t j = 0; j < m; j++) {
		uint32_t i = stencil_array[j];
		V[    j] = 1.0;
		V[  m+j] = (cells.cx[i-1]-x0)/h;
		V[2*m+j] = (cells.cy[i-1]-y0)/h;
		V[3*m+j] = (cells.camb[i-1 + 0*cells.nc] - 2*cells.cx[i-1]*x0 + x0*x0) / (h*h);
		V[4*m+j] = (cells.camb[i-1 + 1*cells.nc] - cells.cx[i-1]*y0 - cells.cy[i-1]*x0 + x0*y0) / (h*h);
		V[5*m+j] = (cells.camb[i-1 + 2*cells.nc] - 2*cells.cy[i-1]*y0 + y0*y0) / (h*h);
	}

	// argomenti in uscita
    plhs[0] = mxCreateNumericMatrix(m,6,mxDOUBLE_CLASS,mxREAL);
    double *V_matrix = mxGetDoubles(plhs[0]);
    if (V_matrix == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "V_matrix is NULL");
    memcpy(V_matrix, V.data(), V.size()*sizeof(double));
}
