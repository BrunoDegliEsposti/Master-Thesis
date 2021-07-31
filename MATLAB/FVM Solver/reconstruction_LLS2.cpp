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
#include "polymesh_FVM.h"
#include "reconstruction.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [up, um] = reconstruction_LLS2(vertices, edges, cells)
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

	std::vector<uint32_t> stencil;
	std::queue<uint32_t> q;

	std::vector<double> V(3*16);
	std::vector<double> e1(3);
	std::vector<double> ubar(16);
	std::vector<double> p(3);

	for (uint32_t i_center = 1; i_center <= cells.nc; i_center++) {

		double x0 = cells.cx[i_center-1];
		double y0 = cells.cy[i_center-1];
		double h = cells.h[i_center-1];

		// costruisci lo stencil intorno a i_center
		build_stencil(i_center, 3, stencil, q, edges, cells);
		size_t m = stencil.size();

		// per ogni variabile fisica
		for (uint32_t l = 0; l < cells.nu; l++)	{

			// riempi V una riga alla volta
			V.reserve(m*3);
			for (size_t j = 0; j < m; j++) {
				uint32_t i = stencil[j];
				V[    j] = 1.0;
				V[  m+j] = (cells.cx[i-1]-x0)/h;
				V[2*m+j] = (cells.cy[i-1]-y0)/h;
			}

			// riempi e1
			e1[0] = 1.0;
			for (size_t j = 1; j < 3; j++) {
				e1[j] = 0.0;
			}

			// riempi ubar
			ubar.reserve(m);
			for (size_t j = 0; j < m; j++) {
				uint32_t i = stencil[j];
				ubar[j] = cells.u[i-1 + l*cells.nc];
			}

			// risolvi il sistema ai minimi quadrati V * p = ubar con vincolo p[0] = ubar[0]
			// il contenuto di V, e1 e ubar viene sovrascritto da lapack
			lapack_int info;
			info = least_squares_constrained(V, e1, ubar, p, m, 3);
			if (info != 0) {
				mexErrMsgIdAndTxt("MEX:FVM_error", "Lapack dgels() failed");
			}

			// per ogni spigolo della cella al centro dello stencil
			uint32_t ne = cells.ne[i_center-1];
			for (uint32_t j = 0; j < ne; j++) {
				int32_t e = cells.e[i_center-1 + j*cells.nc];
				if (e > 0) {
					// per ogni nodo di quadratura
					for (uint32_t k = 0; k < edges.nq; k++) {
						// valuta p(csi,eta) nel nodo
						double t = edges.qx[k];
						Vec2 node = edge_lerp(t,vertices,edges,e);
						double csi = (node.x - x0) / h;
						double eta = (node.y - y0) / h;
						size_t idx = e-1 + l*edges.ne + k*edges.ne*cells.nu;
						up[idx] = p[0] + p[1]*csi + p[2]*eta;
					}
				} else if (e < 0) {
					// per ogni nodo di quadratura
					for (uint32_t k = 0; k < edges.nq; k++) {
						// valuta p(csi,eta) nel nodo
						double t = edges.qx[k];
						Vec2 node = edge_lerp(t,vertices,edges,-e);
						double csi = (node.x - x0) / h;
						double eta = (node.y - y0) / h;
						size_t idx = -e-1 + l*edges.ne + k*edges.ne*cells.nu;
						um[idx] = p[0] + p[1]*csi + p[2]*eta;
					}
				}
			}
		}
	}
}
