#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include "mex.h"
#include "matrix.h"
#include "polymesh.h"

struct Polygon {
	uint32_t ne;
	std::vector<double> vx;
	std::vector<double> vy;

	Polygon(uint32_t mne): vx(mne,0.0), vy(mne,0.0) {}

	void load(Vertices &vertices, Edges &edges, Cells &cells, uint32_t i)
	{
		ne = cells.ne[i-1];
		for (uint32_t j = 0; j < ne; j++) {
			int32_t e = cells.e[i-1 + j*cells.nc];
			if (e > 0) {
				uint32_t v = edges.v1[e-1];
				vx[j] = vertices.x[v-1];
				vy[j] = vertices.y[v-1];
			} else if (e < 0) {
				uint32_t v = edges.v2[-e-1];
				vx[j] = vertices.x[v-1];
				vy[j] = vertices.y[v-1];
			}
		}
	}
};

inline bool ray_segment_overlap(double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4);

bool segment_overlaps_with_convex_polygon(double x1, double y1, double x2, double y2, Polygon &p);

inline bool point_in_convex_polygon(double x, double y, Polygon &p);

inline bool segments_overlap(double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4);

inline bool ccw(double x1, double y1, double x2, double y2, double x3, double y3);

struct Cone {
// Il vertice del cono è (cx,cy). Il cono è compreso tra la semiretta sinistra
// avente direzione (lx-cx,ly-cy) e la semiretta destra avente direzione (rx-cx,ry-cy).
	double cx, cy;
	double lx, ly;
	double rx, ry;

	Cone(double cx, double cy, double lx, double ly, double rx, double ry):
		cx(cx), cy(cy), lx(lx), ly(ly), rx(rx), ry(ry)
	{
		if (!ccw(lx,ly,cx,cy,rx,ry)) {
			mexErrMsgIdAndTxt("MEX:Geometry_error", "Bad arguments passed to cone constructor");
		}
	}

	void widen(double theta)
	// Allarga l'apertura del cono di theta radianti (utile per assicurarsi che
	// punti sul bordo del cono superino sempre il test di appartenenza)
	{
		double costheta = cos(theta);
		double sintheta = sin(theta);
		// ccw rotation
		double lx_new = cx + costheta * (lx-cx) - sintheta * (ly-cy);
		double ly_new = cy + sintheta * (lx-cx) + costheta * (ly-cy);
		// cw rotation
		double rx_new = cx + costheta * (rx-cx) + sintheta * (ry-cy);
		double ry_new = cy - sintheta * (rx-cx) + costheta * (ry-cy);
		// assignments
		lx = lx_new; ly = ly_new; rx = rx_new; ry = ry_new;
		if (!ccw(lx,ly,cx,cy,rx,ry)) {
			mexErrMsgIdAndTxt("MEX:Geometry_error", "Cone is not in ccw order anymore");
		}
	}

	bool overlaps_with_point(double x, double y)
	{
		return ccw(lx,ly,cx,cy,x,y) && ccw(x,y,cx,cy,rx,ry);
	}

	bool overlaps_with_convex_polygon(const Polygon &p)
	{
		// testa l'appartenenza dei vertici del poligono
		for (uint32_t i = 0; i < p.ne; i++) {
			if (overlaps_with_point(p.vx[i],p.vy[i])) return true;
		}
		// testa le intersezioni con i lati del poligono
		for (uint32_t i = 0; i < p.ne; i++) {
			uint32_t inext = (i+1) % p.ne;
			if (ray_segment_overlap(cx,cy,lx,ly,p.vx[i],p.vy[i],p.vx[inext],p.vy[inext])) return true;
			if (ray_segment_overlap(cx,cy,rx,ry,p.vx[i],p.vy[i],p.vx[inext],p.vy[inext])) return true;
		}
		return false;
	}
};

inline bool ccw(double x1, double y1, double x2, double y2, double x3, double y3)
// controlla se i tre punti (x1,y1), (x2,y2), (x3,y3) si avvolgono in senso antiorario
{
	return (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1) >= 0;
}

inline bool segments_overlap(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
// controlla l'intersezione tra il segmento con estremi (x1,y1), (x2,y2) e il segmento con estremi (x3,y3), (x4,y4)
{
	double A_11 = x2 - x1; double A_12 = x3 - x4;
	double A_21 = y2 - y1; double A_22 = y3 - y4;
	double b_1  = x3 - x1; double b_2  = y3 - y1;

	double det = A_11 * A_22 - A_21 * A_12;
	if (det == 0.0) {
		return false;
	}
	double det1 = b_1 * A_22 - b_2 * A_12;
	double det2 = A_11 * b_2 - A_21 * b_1;
	double x_1 = det1/det;
	double x_2 = det2/det;
	return 0.0 <= x_1 && x_1 <= 1.0 && 0.0 <= x_2 && x_2 <= 1.0;
}

inline bool ray_segment_overlap(double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4)
// controlla l'intersezione tra la semiretta con estremo (x1,y1) e
// direzione (x2-x1,y2-y1) e il segmento con estremi (x3,y3) e (x4,y4).
{
	double A_11 = x2 - x1; double A_12 = x3 - x4;
	double A_21 = y2 - y1; double A_22 = y3 - y4;
	double b_1  = x3 - x1; double b_2  = y3 - y1;

	double det = A_11 * A_22 - A_21 * A_12;
	if (det == 0.0) {
		return false;
	}
	double det1 = b_1 * A_22 - b_2 * A_12;
	double det2 = A_11 * b_2 - A_21 * b_1;
	double x_1 = det1/det;
	double x_2 = det2/det;
	return 0.0 <= x_1 && 0.0 <= x_2 && x_2 <= 1.0;
}

inline bool point_in_convex_polygon(double x, double y, const Polygon &p)
// controlla se il punto (x,y) si trova all'interno del poligono convesso p,
// i cui vertici sono memorizzati in senso antiorario
{
	for (uint32_t i = 0; i < p.ne-1; i++) {
		if (!ccw(x,y,p.vx[i],p.vy[i],p.vx[i+1],p.vy[i+1])) {
			return false;
		}
	}
	if (!ccw(x,y,p.vx[p.ne-1],p.vy[p.ne-1],p.vx[0],p.vy[0])) {
		return false;
	}
	return true;
}

bool segment_overlaps_with_convex_polygon(double x1, double y1, double x2, double y2, const Polygon &p)
// controlla se il segmento con estremi (x1,y1), (x2,y2) interseca il poligono p
{
	if (point_in_convex_polygon(x1,y1,p)) return true;
	if (point_in_convex_polygon(x2,y2,p)) return true;
	for (uint32_t i = 0; i < p.ne-1; i++) {
		if (segments_overlap(x1,y1,x2,y2,p.vx[i],p.vy[i],p.vx[i+1],p.vy[i+1])) return true;
	}
	if (segments_overlap(x1,y1,x2,y2,p.vx[p.ne-1],p.vy[p.ne-1],p.vx[0],p.vy[0])) return true;
	return false;
}
