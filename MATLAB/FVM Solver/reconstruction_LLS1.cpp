#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include "mex.h"
#include "matrix.h"
#include "polymesh_FVM.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [up, um] = reconstruction_LLS1(vertices, edges, cells)
{
	Vertices vertices(prhs[0]);
	Edges edges(prhs[1]);
	Cells cells(prhs[2]);

	const mwSize up_dims[] = {edges.ne, cells.nu, edges.nq};
	plhs[0] = mxCreateNumericArray(3, up_dims, mxDOUBLE_CLASS, mxREAL);
	double *up = mxGetDoubles(plhs[0]);

	const mwSize um_dims[] = {edges.ne, cells.nu, edges.nq};
	plhs[1] = mxCreateNumericArray(3, um_dims, mxDOUBLE_CLASS, mxREAL);
	double *um = mxGetDoubles(plhs[1]);

	int32_t e;
	for (uint32_t i = 1; i <= cells.nc; i++) {
		for (uint32_t j = 0; j < cells.mne; j++) {
			e = cells.e[i-1 + j*cells.nc];
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
}
