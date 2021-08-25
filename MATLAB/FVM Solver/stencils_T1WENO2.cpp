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
// function [stencils] = stencils_T1WENO(vertices, edges, cells, i_center)
{
	VerticesFVM vertices(prhs[0]);
	EdgesFVM edges(prhs[1]);
	CellsFVM cells(prhs[2]);
	uint32_t i_center = (uint32_t)mxGetScalar(prhs[3]);

	double x0 = cells.cx[i_center-1];
	double y0 = cells.cy[i_center-1];
	uint32_t ne = cells.ne[i_center-1];

	std::queue<uint32_t> q;
	StencilsVector stencils;
	Stencil stencil_centered;
	build_stencil(i_center, 3, stencil_centered, q, edges, cells);
	stencils.push_back(stencil_centered);
	for (uint32_t j = 0; j < ne; j++) {
		// indice della cella a sinistra rispetto a (x0,y0)
		uint32_t jnext = (j+1) % ne;
		int32_t e_left = cells.e[i_center-1 + jnext*cells.nc];				
		uint32_t i_left = 0;
		if (e_left > 0) {
			i_left = edges.cm[e_left-1];
		} else if (e_left < 0) {
			i_left = edges.cp[-e_left-1];
		}
		if (i_left == 0) continue;

		// indice della cella a destra rispetto a (x0,y0)
		int32_t e_right = cells.e[i_center-1 + j*cells.nc];
		uint32_t i_right = 0;
		if (e_right > 0) {
			i_right = edges.cm[e_right-1];
		} else if (e_right < 0) {
			i_right = edges.cp[-e_right-1];
		}
		if (i_right == 0) continue;

		// costruzione del cono che fa da filtro allo stencil
		double lx = cells.cx[i_left-1];
		double ly = cells.cy[i_left-1];
		double rx = cells.cx[i_right-1];
		double ry = cells.cy[i_right-1];
		Cone cone(x0, y0, lx, ly, rx, ry);
		cone.widen(1e-8);

		// costruzione dello stencil
		auto stencil = build_stencil_from_cone(i_center, 3, cone, q, vertices, edges, cells);
		if (stencil.size() >= 3) {
			stencils.push_back(stencil);
		}
	}

	size_t n = stencils.size();
	size_t m = 0;
	for (const Stencil& stencil: stencils) {
		m = std::max(m,stencil.size());
	}

	// argomenti in uscita
    plhs[0] = mxCreateNumericMatrix(m,n,mxUINT32_CLASS,mxREAL);
    uint32_t *stencils_matrix = mxGetUint32s(plhs[0]);
    if (stencils_matrix == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "stencils_matrix is NULL");

	for (uint32_t j = 0; j < n; j++) {
		for (uint32_t i = 0; i < m; i++) {
			if (stencils[j].size() > i) {
				stencils_matrix[i+j*m] = stencils[j][i];
			}
		}
	}
}
