#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <queue>
#include <algorithm>
#include "mex.h"
#include "matrix.h"
#include "omp.h"
#include "lapacke.h"
#include "polymesh_FVM.h"
#include "method.h"
#include "reconstruction.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [up, um] = reconstruction_T1WENO(vertices, edges, cells, method)
// ricostruzione Type-1 WENO (approccio con penalizzazione)
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

	if (method.order > 3) {
		mexErrMsgIdAndTxt("MEX:FVM_error", "Method.order can be at most 3");
	}

	#pragma omp parallel num_threads(8)
	{
		// variabili di lavoro private di ogni thread
		std::vector<Stencil> stencils;
		std::queue<uint32_t> q;
		std::vector<double> V(6*32);
		std::vector<double> ubar(32);
		std::vector<double> p(6);
		std::vector<double> tau(6);
		
		// per ogni cella, in parallelo
		#pragma omp for
		for (uint32_t i_center = 1; i_center <= cells.nc; i_center++) {
			uint32_t dfb = cells.dfb[i_center-1];
			if (method.order >= 3 && dfb >= 2) {
				reconstruction_T1WENO3(i_center, vertices, edges, cells, up, um, stencils, q, V, ubar, p, tau);
			} else if (method.order >= 2 && dfb >= 1) {
				reconstruction_T1WENO2(i_center, vertices, edges, cells, up, um, stencils, q, V, ubar, p);
			} else {
				reconstruction_LLS1(i_center, vertices, edges, cells, up, um);
			}
		}
	}
}
