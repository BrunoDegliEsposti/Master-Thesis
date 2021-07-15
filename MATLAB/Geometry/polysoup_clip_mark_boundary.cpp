#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <queue>
#include <unordered_set>
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

struct PairHash {
    size_t operator()(std::pair<uint32_t,uint32_t> xy) const noexcept {
        return (size_t)xy.first << 32 | (size_t)xy.second;
    }
};

bool segment_overlaps_with_convex_polygon(double x1, double y1, double x2, double y2, Polygon &p);

inline bool point_in_convex_polygon(double x, double y, Polygon &p);

inline bool segments_overlap(double x1, double y1, double x2, double y2,
	double x3, double y3, double x4, double y4);

inline bool ccw(double x1, double y1, double x2, double y2, double x3, double y3);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [marker] = polysoup_clip_mark_boundary(vertices,edges,cells,crx,cry);
{
	// argomenti in ingresso
    Vertices vertices(prhs[0]);
    Edges edges(prhs[1]);
    Cells cells(prhs[2]);
    double *crx = mxGetDoubles(prhs[3]);
    if (crx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "crx is NULL");
    double *cry = mxGetDoubles(prhs[4]);
    if (cry == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cry is NULL");
    size_t ncr = mxGetNumberOfElements(prhs[3]) - 1;	// crx[0] == crx[ncr]

    // argomenti in uscita
    plhs[0] = mxCreateNumericMatrix(cells.nc,1,mxUINT8_CLASS,mxREAL);
    uint8_t *marker = mxGetUint8s(plhs[0]);
    if (marker == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "marker is NULL");

    // trova la prima coppia segmento-poligono che si interseca
    Polygon p(cells.mne);
    uint32_t i_first = 0;
    for (uint32_t i = 1; i <= cells.nc; i++) {
    	p.load(vertices,edges,cells,i);
    	if (segment_overlaps_with_convex_polygon(crx[0],cry[0],crx[1],cry[1],p)) {
    		i_first = i;
    		break;
    	}
    }
    if (i_first == 0) mexErrMsgIdAndTxt("MEX:input_error", "clipping region is not a subset of polymesh");

    // inizializza le strutture dati di supporto
    std::queue<std::pair<uint32_t,uint32_t>> q;
    std::unordered_set<std::pair<uint32_t,uint32_t>,PairHash> us;
    q.push(std::make_pair(i_first,0));
    us.insert(std::make_pair(i_first,0));
    marker[i_first-1] = 1;

    // marca tutti i poligoni che intersecano il bordo della clipping region
    // sfruttando le informazioni di adiacenza sia dei poligoni che del bordo
    while (!q.empty()) {
    	uint32_t i = q.front().first;		// indice del poligono (1-based)
    	uint32_t k = q.front().second;		// indice del segmento (0-based)
    	
    	uint32_t i_new = i;
    	uint32_t k_new = k;
    	p.load(vertices,edges,cells,i_new);
    	for (int32_t l : {-1,1}) {
	    	k_new = (k+l) % ncr;
	    	if (!segment_overlaps_with_convex_polygon(
	    		crx[k_new],cry[k_new],crx[k_new+1],cry[k_new+1],p)) continue;
	    	auto search = us.find(std::make_pair(i_new,k_new));
	    	if (search == us.end()) {
	    		q.push(std::make_pair(i_new,k_new));
	    		us.insert(std::make_pair(i_new,k_new));
	    		marker[i_new-1] = 1;
	    	}
	    }

    	for (uint32_t j = 0; j < cells.ne[i-1]; j++) {
    		int32_t e = cells.e[i-1+j*cells.nc];
    		if (e > 0) {
    			i_new = edges.cm[e-1];
    		} else if (e < 0) {
    			i_new = edges.cp[-e-1];
    		} else {
    			mexErrMsgIdAndTxt("MEX:mesh_error", "cells.ne is wrong");
    		}
    		if (i_new == 0) {
    			// lo spigolo era di bordo
    			continue;
    		}
    		p.load(vertices,edges,cells,i_new);
    		for (int32_t l : {-1,0,1}) {
    			k_new = (k+l) % ncr;
    			if (!segment_overlaps_with_convex_polygon(
    				crx[k_new],cry[k_new],crx[k_new+1],cry[k_new+1],p)) continue;
    			auto search = us.find(std::make_pair(i_new,k_new));
    			if (search == us.end()) {
    				q.push(std::make_pair(i_new,k_new));
    				us.insert(std::make_pair(i_new,k_new));
    				marker[i_new-1] = 1;
    			}
    		}
    	}
    	q.pop();
    }
}

bool segment_overlaps_with_convex_polygon(double x1, double y1, double x2, double y2, Polygon &p)
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

inline bool point_in_convex_polygon(double x, double y, Polygon &p)
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
