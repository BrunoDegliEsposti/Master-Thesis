#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "mex.h"
#include "matrix.h"
#include "omp.h"
#include "polymesh_FVM.h"
#include "method.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [up, um] = reconstruction_LLS1(vertices, edges, cells, method)
{
	VerticesFVM vertices(prhs[0]);
	EdgesFVM edges(prhs[1]);
	CellsFVM cells(prhs[2]);
	Method method(prhs[3]);

	const mwSize up_dims[] = {edges.ne, cells.nu, edges.nq};
	plhs[0] = mxCreateNumericArray(3, up_dims, mxDOUBLE_CLASS, mxREAL);
	double *up = mxGetDoubles(plhs[0]);

	const mwSize um_dims[] = {edges.ne, cells.nu, edges.nq};
	plhs[1] = mxCreateNumericArray(3, um_dims, mxDOUBLE_CLASS, mxREAL);
	double *um = mxGetDoubles(plhs[1]);

	if (method.order != 1) {
		mexErrMsgIdAndTxt("MEX:FVM_error", "Method.order must be set to 1");
	}

	#pragma omp parallel num_threads(8)
	{

	#pragma omp for
	for (uint32_t i = 1; i <= cells.nc; i++) {
		for (uint32_t j = 0; j < cells.mne; j++) {
			int32_t e = cells.e[i-1 + j*cells.nc];
			if (e > 0) {
				for (uint32_t k = 0; k < edges.nq; k++) {
					for (uint32_t l = 0; l < cells.nu; l++) {
						size_t idx = e-1 + l*edges.ne + k*edges.ne*cells.nu;
						up[idx] = cells.u[i-1 + l*cells.nc];
					}
				}
			} else if (e < 0) {
				for (uint32_t k = 0; k < edges.nq; k++) {
					for (uint32_t l = 0; l < cells.nu; l++) {
						size_t idx = -e-1 + l*edges.ne + k*edges.ne*cells.nu;
						um[idx] = cells.u[i-1 + l*cells.nc];
					}
				}
			}
		}
	}

	// fine del blocco omp parallel
	}
}
