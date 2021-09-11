#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <queue>
#include <algorithm>
#include "mex.h"
#include "matrix.h"
#include "polymesh_FVM.h"
#include "reconstruction.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [stencil] = stencil_centered(vertices, edges, cells, i_center, min_stencil_size)
{
	Vertices vertices(prhs[0]);
	Edges edges(prhs[1]);
	Cells cells(prhs[2]);
	uint32_t i_center = (uint32_t)mxGetScalar(prhs[3]);
	uint32_t min_stencil_size = (uint32_t)mxGetScalar(prhs[4]);

	Stencil stencil;
	std::queue<uint32_t> q;
	build_centered_stencil(i_center, min_stencil_size, stencil, q, edges, cells);
	size_t m = stencil.size();

	// argomenti in uscita
    plhs[0] = mxCreateNumericMatrix(m,1,mxUINT32_CLASS,mxREAL);
    uint32_t *stencil_matrix = mxGetUint32s(plhs[0]);
    if (stencil_matrix == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "stencil_matrix is NULL");

	for (uint32_t i = 0; i < m; i++) {
		stencil_matrix[i] = stencil[i].idx;
	}
}
