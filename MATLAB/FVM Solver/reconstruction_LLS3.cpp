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
#include "method.h"
#include "reconstruction.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [up, um] = reconstruction_LLS3(vertices, edges, cells, method)
// ricostruzione Linear Least-Squares di grado 3
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

	if (method.order != 3) {
		mexErrMsgIdAndTxt("MEX:FVM_error", "Method.order must be set to 3");
	}

	char lst = method.least_squares_type;
	if (lst != 'U' && lst != 'C' && lst != 'P') {
		mexErrMsgIdAndTxt("MEX:FVM_error", "Method.least_squares_type must be U/C/P");
	}

	#pragma omp parallel num_threads(8)
	{

	// variabili di lavoro private di ogni thread
	Stencil stencil;
	std::queue<uint32_t> q;
	std::vector<double> V(6*32);
	std::vector<double> V_first_row(6);
	std::vector<double> ubar(32);
	std::vector<double> p(6);

	#pragma omp for
	for (uint32_t i_center = 1; i_center <= cells.nc; i_center++) {

		double x0 = cells.cx[i_center-1];
		double y0 = cells.cy[i_center-1];
		double h = cells.h[i_center-1];

		// costruisci lo stencil intorno a i_center
		bool success;
		uint32_t min_stencil_size = 6;
		if (cells.dfb[i_center-1] == 0) {
			// La dimensione minima dello stencil sulle celle di bordo
			// viene incrementata per questioni di stabilità
			min_stencil_size = 7;
		}
		success = build_centered_stencil(i_center, min_stencil_size, stencil, q, edges, cells);
		if (!success) {
			mexErrMsgIdAndTxt("MEX:FVM_error", "Can't build a large enough stencil");
		}
		size_t m = stencil.size();

		// per ogni variabile fisica
		for (uint32_t l = 0; l < cells.nu; l++)	{

			// riempi V una riga alla volta
			V.reserve(m*6);
			for (size_t j = 0; j < m; j++) {
				uint32_t i = stencil[j].idx;
				V[    j] = 1.0;
				V[  m+j] = (cells.cx[i-1]-x0)/h;
				V[2*m+j] = (cells.cy[i-1]-y0)/h;
				V[3*m+j] = (cells.camb[i-1 + 0*cells.nc] - 2*cells.cx[i-1]*x0 + x0*x0) / (h*h);
				V[4*m+j] = (cells.camb[i-1 + 1*cells.nc] - cells.cx[i-1]*y0 - cells.cy[i-1]*x0 + x0*y0) / (h*h);
				V[5*m+j] = (cells.camb[i-1 + 2*cells.nc] - 2*cells.cy[i-1]*y0 + y0*y0) / (h*h);
			}

			// riempi V_first_row (utile solo se least_squares_type == 'C')
			for (size_t j = 0; j < 6; j++) {
				V_first_row[j] = V[j*m];
			}

			// riempi ubar
			ubar.reserve(m);
			for (size_t j = 0; j < m; j++) {
				uint32_t i = stencil[j].idx;
				ubar[j] = cells.u[i-1 + l*cells.nc];
			}

			// ricostruzione ai minimi quadrati
			if (method.least_squares_type == 'U') {
				// Unconstrained least squares polynomial reconstruction.
				// Risolvi il sistema ai minimi quadrati V * p = ubar.
				// Il contenuto di V e di ubar viene sovrascritto da lapack.
				lapack_int info;
				info = least_squares(V, ubar, p, m, 6);
				if (info != 0) {
					mexErrMsgIdAndTxt("MEX:FVM_error", "Lapack dgels() failed");
				}
			} else if (method.least_squares_type == 'C') {
				// Constrained least squares polynomial reconstruction.
				// Risolvi il sistema ai minimi quadrati V * p = ubar con il
				// vincolo di stabilità V_first_row * p = ubar[0]. Il contenuto
				// di V, V_first_row e ubar viene sovrascritto da lapack.
				lapack_int info;
				info = least_squares_constrained(V, V_first_row, ubar, p, m, 6);
				if (info != 0) {
					mexErrMsgIdAndTxt("MEX:FVM_error", "Lapack dgglse() failed");
				}
			} else {
				// Penalized least squares polynomial reconstruction.
				// Assegna pesi diversi alle varie equazioni che compongono il sistema
				// ai minimi quadrati in base alla distanza dal centro dello stencil.
				for (size_t i = 0; i < m; i++) {
					uint32_t d = stencil[i].distance;
					double weight = pow(10,-(double)d);
					for (size_t j = 0; j < 6; j++) {
						V[i+j*m] *= weight;
					}
					ubar[i] *= weight;
				}
				// Risolvi il sistema ai minimi quadrati V * p = ubar.
				// Il contenuto di V e di ubar viene sovrascritto da lapack.
				lapack_int info;
				info = least_squares(V, ubar, p, m, 6);
				if (info != 0) {
					mexErrMsgIdAndTxt("MEX:FVM_error", "Weighted lapack dgels() failed");
				}
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
						up[idx] = p[0] + p[1]*csi + p[2]*eta + p[3]*csi*csi + p[4]*csi*eta + p[5]*eta*eta;
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
						um[idx] = p[0] + p[1]*csi + p[2]*eta + p[3]*csi*csi + p[4]*csi*eta + p[5]*eta*eta;
					}
				}
			}
		}
	}
	
	// fine del blocco omp parallel
	}
}
