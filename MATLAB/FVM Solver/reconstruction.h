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
#include "../Geometry/computational_geometry.h"

// https://en.cppreference.com/w/cpp/container/vector
// https://en.cppreference.com/w/cpp/container/queue
// https://www.netlib.org/lapack/lug/node27.html
// https://www.netlib.org/lapack/lapacke.html#_calling_code_dgels_code
// https://it.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html

typedef std::vector<uint32_t> Stencil;
typedef std::vector<Stencil> StencilsVector;

void build_stencil(uint32_t i_center, uint32_t min_stencil_size,
	std::vector<uint32_t> &stencil, std::queue<uint32_t> &q, EdgesFVM &edges, CellsFVM &cells)
{
	stencil.clear();
	stencil.push_back(i_center);
	// clear queue
	while (!q.empty()) {
		q.pop();
	}
	q.push(i_center);

	while (stencil.size() < min_stencil_size) {
		uint32_t n = q.size();
		if (n == 0) {
			mexErrMsgIdAndTxt("MEX:FVM_error", "Can't build a large enough stencil");
		}
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
				auto search = std::find(stencil.begin(), stencil.end(), i_new);
				if (search == std::end(stencil)) {
					stencil.push_back(i_new);
					q.push(i_new);
				}
			}
			q.pop();
		}
	}
}

Stencil build_stencil_from_cone(uint32_t i_center, uint32_t min_stencil_size, Cone &cone,
	std::queue<uint32_t> &q, VerticesFVM &vertices, EdgesFVM &edges, CellsFVM &cells)
{
	Stencil stencil = {i_center};
	// clear queue
	while (!q.empty()) {
		q.pop();
	}
	q.push(i_center);

	while (stencil.size() < min_stencil_size) {
		uint32_t n = q.size();
		if (n == 0) {
			// non è stato possibile costruire uno stencil grande a sufficienza;
			// ritorna comunque quello che c'è. stencil.size() sarà minore di min_stencil_size.
			return stencil;
		}
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
				auto search = std::find(stencil.begin(), stencil.end(), i_new);
				if (search == std::end(stencil)) {
					stencil.push_back(i_new);
					q.push(i_new);
				}
			}
			q.pop();
		}
	}
	return stencil;
}

lapack_int least_squares(std::vector<double> &V, std::vector<double> &p,
	std::vector<double> &ubar, lapack_int m, lapack_int n)
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

lapack_int least_squares_constrained(std::vector<double> &V, std::vector<double> &e1,
	std::vector<double> &ubar, std::vector<double> &p, lapack_int m, lapack_int n)
{
	double *a = V.data();
	lapack_int lda = m;
	double *b = e1.data();
	lapack_int ldb = 1;
	double *c = ubar.data();
	double ubar_center = ubar[0];
	double *d = &ubar_center;
	double *x = p.data();
	return LAPACKE_dgglse(LAPACK_COL_MAJOR, m, n, 1, a, lda, b, ldb, c, d, x);
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
