#pragma once
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
#include "method.h"
#include "../Geometry/computational_geometry.h"

// https://en.cppreference.com/w/cpp/container/vector
// https://en.cppreference.com/w/cpp/container/queue
// https://www.netlib.org/lapack/lug/node27.html
// https://www.netlib.org/lapack/lapacke.html#_calling_code_dgels_code
// https://it.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html

struct StencilCell {
	uint32_t idx;
	uint32_t distance;
	StencilCell(uint32_t idx, uint32_t distance): idx(idx), distance(distance) {}
};

typedef std::vector<StencilCell> Stencil;

bool build_centered_stencil(uint32_t i_center, uint32_t min_stencil_size,
	Stencil &stencil, std::queue<uint32_t> &q, Edges &edges, Cells &cells)
// Crea uno stencil a macchia d'olio a partire dalla cella con indice i_center.
// A ogni iterazione viene aggiunto un intorno a forma di anello finché non viene
// raggiunta min_stencil_size. In caso di successo la funzione ritorna true,
// altrimenti false. Le distanze sono definite come in teoria dei grafi.
{
	stencil.clear();
	uint32_t distance = 0;
	stencil.push_back(StencilCell(i_center,distance));
	
	// clear queue
	while (!q.empty()) {
		q.pop();
	}
	q.push(i_center);

	while (stencil.size() < min_stencil_size) {
		uint32_t n = q.size();
		if (n == 0) {
			// sono finite le celle da poter aggiungere allo stencil
			return false;
		}
		distance += 1;
		while (n--) {
			uint32_t i = q.front();
			uint32_t ne = cells.ne[i-1];
			for (uint32_t j = 0; j < ne; j++) {
				int32_t e = cells.e[i-1 + j*cells.nc];
				uint32_t i_new = 0;
				if (e > 0) {
					i_new = edges.cm[e-1];
				} else if (e < 0) {
					i_new = edges.cp[-e-1];
				}
				if (i_new == 0) continue;
				// linear search on the stencil to find i_new
				bool found = false;
				for (uint32_t k = 0; k < stencil.size(); k++) {
					if (stencil[k].idx == i_new) {
						found = true;
						break;
					}
				}
				if (!found) {
					stencil.push_back(StencilCell(i_new,distance));
					q.push(i_new);
				}
			}
			q.pop();
		}
	}
	return true;
}

bool build_onesided_stencil(uint32_t i_center, uint32_t min_stencil_size,
	Stencil &stencil, Cone &cone, std::queue<uint32_t> &q,
	Vertices &vertices, Edges &edges, Cells &cells)
// Crea uno stencil nella direzione di cone a partire dalla cella con indice i_center.
// A ogni iterazione viene aggiunto un intorno a forma di settore di anello
// finché non viene raggiunta min_stencil_size. In caso di successo la funzione
// ritorna true, altrimenti false. Le distanze sono definite come in teoria dei grafi.
{
	stencil.clear();
	uint32_t distance = 0;
	stencil.push_back(StencilCell(i_center,distance));

	// clear queue
	while (!q.empty()) {
		q.pop();
	}
	q.push(i_center);

	while (stencil.size() < min_stencil_size) {
		uint32_t n = q.size();
		if (n == 0) {
			// sono finite le celle da poter aggiungere allo stencil
			return false;
		}
		distance += 1;
		while (n--) {
			uint32_t i = q.front();
			uint32_t ne = cells.ne[i-1];
			for (uint32_t j = 0; j < ne; j++) {
				int32_t e = cells.e[i-1 + j*cells.nc];
				uint32_t i_new = 0;
				if (e > 0) {
					i_new = edges.cm[e-1];
				} else if (e < 0) {
					i_new = edges.cp[-e-1];
				}
				if (i_new == 0) continue;
				Polygon p(cells.mne);
				p.load(vertices,edges,cells,i_new);
				if (!cone.overlaps_with_convex_polygon(p)) continue;
				// linear search on the stencil to find i_new
				bool found = false;
				for (uint32_t k = 0; k < stencil.size(); k++) {
					if (stencil[k].idx == i_new) {
						found = true;
						break;
					}
				}
				if (!found) {
					stencil.push_back(StencilCell(i_new,distance));
					q.push(i_new);
				}
			}
			q.pop();
		}
	}
	return true;
}

lapack_int least_squares(std::vector<double> &V, std::vector<double> &ubar,
	std::vector<double> &p, lapack_int m, lapack_int n)
{
	lapack_int nrhs = 1;
	double *a = V.data();
	lapack_int lda = m;
	double *b = ubar.data();
	lapack_int ldb = m;
	lapack_int info;
	info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', m, n, nrhs, a, lda, b, ldb);
	for (lapack_int i = 0; i < n; i++) {
		p[i] = b[i];
	}
	return info;
}

lapack_int least_squares_constrained(std::vector<double> &V, std::vector<double> &V_first_row,
	std::vector<double> &ubar, std::vector<double> &p, lapack_int m, lapack_int n)
{
	double *a = V.data();
	lapack_int lda = m;
	double *b = V_first_row.data();
	lapack_int ldb = 1;
	double *c = ubar.data();
	double ubar_center = ubar[0];
	double *d = &ubar_center;
	double *x = p.data();
	return LAPACKE_dgglse(LAPACK_COL_MAJOR, m, n, 1, a, lda, b, ldb, c, d, x);
}

lapack_int least_squares_and_condition_number(std::vector<double> &V, std::vector<double> &ubar,
	std::vector<double> &p, lapack_int m, lapack_int n, double &condition_number, std::vector<double> &tau)
{
	// compute QR factorization (https://www.netlib.org/lapack/lug/node40.html)
	double *qr = V.data();
	lapack_int lda = m;
	lapack_int info;
	info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, n, qr, lda, tau.data());
	if (info != 0) return info;

	// approximate the condition number of V using R
	double rcond;
	info = LAPACKE_dtrcon(LAPACK_COL_MAJOR, 'I', 'U', 'N', n, qr, lda, &rcond);
	if (info == 0) {
		condition_number = 1/rcond;
	} else {
		return info;
	}

	// compute c = Q^T ubar
	double *c = ubar.data();
	info = LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'T', m, 1, n, qr, lda, tau.data(), c, m);
	if (info != 0) return info;

	// solve the linear system R p = c(1:n)
	info = LAPACKE_dtrtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', n, 1, qr, lda, c, n);
	if (info == 0) {
		memcpy(p.data(), c, n*sizeof(double));
	}
	return info;
}

double smoothness_indicator_WENO2(std::vector<double> p, uint32_t i, CellsFVM &cells)
// calcola l'indicatore di smoothness per il polinomio locale p(\csi,\eta) di grado 1
{
	double p1 = p[1];
	double p2 = p[2];
	return p1*p1 + p2*p2;
}

double smoothness_indicator_WENO3(std::vector<double> p, uint32_t i, CellsFVM &cells)
// calcola l'indicatore di smoothness per il polinomio locale p(\csi,\eta) di grado 2
{
	double p1 = p[1];
	double p2 = p[2];
	double p3 = p[3];
	double p4 = p[4];
	double p5 = p[5];
	double h = cells.h[i-1];
	double x0 = cells.cx[i-1];
	double y0 = cells.cy[i-1];
	double xx_mean = cells.camb[i-1 + 0*cells.nc];
	double yy_mean = cells.camb[i-1 + 2*cells.nc];
	double csicsi_mean = (xx_mean - x0*x0) / (h*h);
	double etaeta_mean = (yy_mean - y0*y0) / (h*h);
	double s = p1*p1 + p2*p2 + 4*p3*p3 + p4*p4 + 4*p5*p5;
	s += (4*p3*p3 + p4*p4) * csicsi_mean;
	s += (4*p5*p5 + p4*p4) * etaeta_mean;
	return s;
}

void reconstruction_LLS1(uint32_t i_center, VerticesFVM &vertices,
	EdgesFVM &edges, CellsFVM &cells, double *up, double *um)
{
	for (uint32_t j = 0; j < cells.mne; j++) {
		int32_t e = cells.e[i_center-1 + j*cells.nc];
		if (e > 0) {
			for (uint32_t k = 0; k < edges.nq; k++) {
				for (uint32_t l = 0; l < cells.nu; l++) {
					size_t idx = e-1 + l*edges.ne + k*edges.ne*cells.nu;
					up[idx] = cells.u[i_center-1 + l*cells.nc];
				}
			}
		} else if (e < 0) {
			for (uint32_t k = 0; k < edges.nq; k++) {
				for (uint32_t l = 0; l < cells.nu; l++) {
					size_t idx = -e-1 + l*edges.ne + k*edges.ne*cells.nu;
					um[idx] = cells.u[i_center-1 + l*cells.nc];
				}
			}
		}
	}
}

void reconstruction_T1WENO2(
	// variabili in ingresso
	uint32_t i_center, VerticesFVM &vertices, EdgesFVM &edges, CellsFVM &cells, Method &method,
	// variabili in uscita
	double *up, double *um,
	// variabili di lavoro
	std::vector<Stencil> &stencils,	std::queue<uint32_t> &q, std::vector<double> &V,
	std::vector<double> &ubar, std::vector<double> &p)
{
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
		bool success = build_centered_stencil(i_center, 3, stencil_centered, q, edges, cells);
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
			Stencil stencil_onesided;
			bool success = build_onesided_stencil(i_center, 3, stencil_onesided, cone, q, vertices, edges, cells);
			if (success) {
				stencils.push_back(stencil_onesided);
			}
		}

		double sum_of_omega_tilde = 0.0;
		std::vector<double> p_weno(3, 0.0);

		// itera sugli stencil trovati per costruire p_weno
		for (const Stencil& stencil: stencils) {
			size_t m = stencil.size();

			// riempi V una riga alla volta
			V.reserve(m*3);
			for (size_t j = 0; j < m; j++) {
				uint32_t i = stencil[j].idx;
				V[    j] = 1.0;
				V[  m+j] = (cells.cx[i-1]-x0)/h;
				V[2*m+j] = (cells.cy[i-1]-y0)/h;
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
				for (size_t j = 0; j < 3; j++) {
					V[i+j*m] *= weight;
				}
				ubar[i] *= weight;
			}
			
			// risolvi il sistema ai minimi quadrati V * p = ubar.
			// il contenuto di V e di ubar viene sovrascritto da lapack
			lapack_int info;
			info = least_squares(V, ubar, p, m, 3);
			if (info != 0) {
				mexErrMsgIdAndTxt("MEX:FVM_error", "Weighted lapack dgels() failed");
			}

			// aggiorna p_weno con una versione pesata del p che abbiamo appena trovato
			double si = smoothness_indicator_WENO2(p, i_center, cells);
			double omega_tilde = pow(method.WENO_epsilon + si, -method.WENO_power);
			sum_of_omega_tilde += omega_tilde;
			for (uint32_t j = 0; j < 3; j++) {
				p_weno[j] += omega_tilde * p[j];
			}
		}

		// normalizza p_weno
		for (uint32_t j = 0; j < 3; j++) {
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
					up[idx] = p_weno[0] + p_weno[1]*csi + p_weno[2]*eta;
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
					um[idx] = p_weno[0] + p_weno[1]*csi + p_weno[2]*eta;
				}
			}
		}
	}
}

void reconstruction_T1WENO3(
	// variabili in ingresso
	uint32_t i_center, VerticesFVM &vertices, EdgesFVM &edges, CellsFVM &cells, Method &method,
	// variabili in uscita
	double *up, double *um,
	// variabili di lavoro
	std::vector<Stencil> &stencils,	std::queue<uint32_t> &q, std::vector<double> &V,
	std::vector<double> &ubar, std::vector<double> &p, std::vector<double> &tau)
{
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
			Stencil stencil_onesided;
			bool success = build_onesided_stencil(i_center, 6, stencil_onesided, cone, q, vertices, edges, cells);
			if (success) {
				stencils.push_back(stencil_onesided);
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
			double omega_tilde = pow(method.WENO_epsilon + si, -method.WENO_power);
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