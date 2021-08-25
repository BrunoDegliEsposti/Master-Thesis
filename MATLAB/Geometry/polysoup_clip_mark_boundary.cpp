#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <queue>
#include <unordered_set>
#include "mex.h"
#include "matrix.h"
#include "polymesh.h"
#include "computational_geometry.h"

struct PairHash {
    size_t operator()(std::pair<uint32_t,uint32_t> xy) const noexcept {
        return (size_t)xy.first << 32 | (size_t)xy.second;
    }
};

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
