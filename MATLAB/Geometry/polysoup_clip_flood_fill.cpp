#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <queue>
#include "mex.h"
#include "matrix.h"
#include "polymesh.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// function [marker] = polysoup_clip_flood_fill(vertices,edges,cells,marker,starting_cell);
{
	// argomenti in ingresso
    Vertices vertices(prhs[0]);
    Edges edges(prhs[1]);
    Cells cells(prhs[2]);
    uint32_t starting_cell = (uint32_t)mxGetScalar(prhs[4]);

    // argomenti in uscita
    plhs[0] = mxDuplicateArray(prhs[3]);
    uint8_t *marker = mxGetUint8s(plhs[0]);

    // algoritmo flood fill
    uint8_t value = marker[starting_cell-1];
    std::queue<uint32_t> q;
    q.push(starting_cell);
    while (!q.empty()) {
    	uint32_t i = q.front();
    	q.pop();
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
			// nessuna cella adiacente (siamo sul bordo)
			if (c == 0) continue;
			// se la cella adiacente è già marcata, saltala
			if (marker[c-1] != 0) continue;
			// altrimenti, marcala e aggiungila alla coda
			marker[c-1] = value;
			q.push(c);
		}
    }
}