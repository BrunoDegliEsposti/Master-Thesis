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
#include "reconstruction.h"
#include "../Geometry/computational_geometry.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [up, um] = reconstruction_T1WENO3(vertices, edges, cells)
// ricostruzione Type-1 WENO di grado 3 (approccio con penalizzazione)
{
	VerticesFVM vertices(prhs[0]);
	EdgesFVM edges(prhs[1]);
	CellsFVM cells(prhs[2]);

	const mwSize up_dims[] = {edges.ne, cells.nu, edges.nq};
	plhs[0] = mxCreateNumericArray(3, up_dims, mxDOUBLE_CLASS, mxREAL);
	double *up = mxGetDoubles(plhs[0]);

	const mwSize um_dims[] = {edges.ne, cells.nu, edges.nq};
	plhs[1] = mxCreateNumericArray(3, um_dims, mxDOUBLE_CLASS, mxREAL);
	double *um = mxGetDoubles(plhs[1]);

	#pragma omp parallel num_threads(8)
	{

	// variabili di lavoro private di ogni thread
	std::vector<Stencil> stencils;
	std::queue<uint32_t> q;
	std::vector<double> V(6*32);
	std::vector<double> ubar(32);
	std::vector<double> p(6);
	std::vector<double> tau(6);
	
	#pragma omp for
	// per ogni cella, in parallelo
	for (uint32_t i_center = 1; i_center <= cells.nc; i_center++) {

		double x0 = cells.cx[i_center-1];
		double y0 = cells.cy[i_center-1];
		double h = cells.h[i_center-1];
		uint32_t ne = cells.ne[i_center-1];

		// per ogni variabile fisica
		for (uint32_t l = 0; l < cells.nu; l++)	{

			// costruisci più stencil intorno a i_center, di cui uno
			// centrato e altri invece decentrati a forma di cono
			stencils.clear();
			Stencil stencil_centered;
			bool success = build_centered_stencil(i_center, 6, stencil_centered, q, edges, cells);
			if (!success) {
				mexErrMsgIdAndTxt("MEX:FVM_error", "Can't build a large enough centered stencil");
			}
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

				// costruzione dello stencil nella direzione del cono
				Stencil stencil_biased;
				bool success = build_biased_stencil(i_center, 6, stencil_biased, cone, q, vertices, edges, cells);
				if (success) {
					stencils.push_back(stencil_biased);
				}
			}

			double sum_of_omega_tilde = 0.0;
			std::vector<double> p_weno(6, 0.0);

			// itera sugli stencil trovati per costruire p_weno
			for (const Stencil& stencil: stencils) {
				size_t m = stencil.size();

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

				// riempi ubar
				ubar.reserve(m);
				for (size_t j = 0; j < m; j++) {
					uint32_t i = stencil[j].idx;
					ubar[j] = cells.u[i-1 + l*cells.nc];
				}

				// assegna pesi diversi alle varie equazioni che compongono il sistema
				// ai minimi quadrati in base alla distanza dal centro dello stencil.
				for (size_t i = 0; i < m; i++) {
					uint32_t d = stencil[i].distance;
					double weight = pow(10,-(double)d);
					for (size_t j = 0; j < 6; j++) {
						V[i+j*m] *= weight;
					}
					ubar[i] *= weight;
				}
				
				// risolvi il sistema ai minimi quadrati V * p = ubar.
				// il contenuto di V e di ubar viene sovrascritto da lapack
				double condition_number;
				lapack_int info;
				info = least_squares_and_condition_number(V, ubar, p, m, 6, condition_number, tau);
				if (info != 0 || condition_number > 1e3) {
					//mexPrintf("Cell number %u\n", (unsigned int)i_center);
					//for (uint32_t i = 0; i < stencil.size(); i++) {
					//	mexPrintf("%u\n",stencil[i].idx);
					//}
					//mexPrintf("Value of INFO: %u\n", (unsigned int)info);
					//mexPrintf("Value of cond(V): %f\n", condition_number);
					continue;
					//mexErrMsgIdAndTxt("MEX:FVM_error", "Weighted lapack dgels() failed");
				}

				// aggiorna p_weno con una versione pesata del p che abbiamo appena trovato
				double si = smoothness_indicator_WENO3(p, i_center, cells);
				double omega_tilde = pow(1e-5 + si, -4.0);
				sum_of_omega_tilde += omega_tilde;
				for (uint32_t j = 0; j < 6; j++) {
					p_weno[j] += omega_tilde * p[j];
				}
			}

			if (sum_of_omega_tilde == 0.0) {
				mexErrMsgIdAndTxt("MEX:FVM_error", "Not even one stencil was good enough");
			}

			// normalizza p_weno
			for (uint32_t j = 0; j < 6; j++) {
				p_weno[j] /= sum_of_omega_tilde;
			}

			// per ogni spigolo della cella al centro dello stencil
			for (uint32_t j = 0; j < ne; j++) {
				int32_t e = cells.e[i_center-1 + j*cells.nc];
				if (e > 0) {
					// per ogni nodo di quadratura
					for (uint32_t k = 0; k < edges.nq; k++) {
						// valuta p_weno(csi,eta) in tale nodo
						double t = edges.qx[k];
						Vec2 node = edge_lerp(t,vertices,edges,e);
						double csi = (node.x - x0) / h;
						double eta = (node.y - y0) / h;
						size_t idx = e-1 + l*edges.ne + k*edges.ne*cells.nu;
						up[idx] = p_weno[0] + p_weno[1]*csi + p_weno[2]*eta
						        + p_weno[3]*csi*csi + p_weno[4]*csi*eta + p_weno[5]*eta*eta;
					}
				} else if (e < 0) {
					// per ogni nodo di quadratura
					for (uint32_t k = 0; k < edges.nq; k++) {
						// valuta p_weno(csi,eta) in tale nodo
						double t = edges.qx[k];
						Vec2 node = edge_lerp(t,vertices,edges,-e);
						double csi = (node.x - x0) / h;
						double eta = (node.y - y0) / h;
						size_t idx = -e-1 + l*edges.ne + k*edges.ne*cells.nu;
						um[idx] = p_weno[0] + p_weno[1]*csi + p_weno[2]*eta
						        + p_weno[3]*csi*csi + p_weno[4]*csi*eta + p_weno[5]*eta*eta;
					}
				}
			}
		}
	}

	// fine del blocco omp parallel
	}
}
