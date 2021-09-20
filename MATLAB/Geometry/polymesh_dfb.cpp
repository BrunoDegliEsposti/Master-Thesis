#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <queue>
#include "mex.h"
#include "matrix.h"
#include "polymesh.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [dfb] = polymesh_dfb(vertices,edges,cells);
{
	// argomenti in ingresso
    Vertices vertices(prhs[0]);
    Edges edges(prhs[1]);
    Cells cells(prhs[2]);

    // argomenti in uscita
    uint32_t nc = cells.nc;
    plhs[0] = mxCreateNumericMatrix(nc, 1, mxUINT32_CLASS, mxREAL);
    uint32_t *dfb = mxGetUint32s(plhs[0]);

    // inizializza dfb con distanze massime
    for (uint32_t i = 1; i <= nc; i++) {
    	dfb[i-1] = UINT32_MAX;
    }

    // la distanza dal bordo viene calcolata con un algoritmo
    // flood-fill a partire dal bordo. L'algoritmo usa una coda:
    std::queue<uint32_t> q;

    // inizializza la coda iterando sugli spigoli di bordo.
    // eventuali duplicati nella coda non sono un problema.
    for (uint32_t j = edges.nie+1; j <= edges.ne; j++) {
    	uint32_t i = edges.cp[j-1];
    	dfb[i-1] = 0;
    	q.push(i);
    }

    // algoritmo flood fill
    while (!q.empty()) {
    	uint32_t i = q.front();
    	q.pop();
    	uint32_t distance_i = dfb[i-1];
		uint32_t ne = cells.ne[i-1];
		// per ogni cella adiacente c
		uint32_t c;
		for (uint32_t j = 0; j < ne; j++) {
			int32_t e = cells.e[i-1 + j*cells.nc];
			if (e > 0) {
				c = edges.cm[e-1];
			} else if (e < 0) {
				c = edges.cp[-e-1];
			} else {
				continue;
			}
			if (c == 0) continue;
			uint32_t distance_c = dfb[c-1];
			if (distance_i + 1 < distance_c) {
				dfb[c-1] = distance_i + 1;
				q.push(c);
			}
		}
    }
}
