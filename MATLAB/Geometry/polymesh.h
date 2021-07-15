#pragma once
#include <cstdint>
#include "mex.h"
#include "matrix.h"

struct Vertices {
	uint32_t nv;
	double *x;
	double *y;
	Vertices(const mxArray *in)
	{
		nv = (uint32_t)mxGetScalar(mxGetField(in,0,"nv"));
		x = mxGetDoubles(mxGetField(in,0,"x"));
		y = mxGetDoubles(mxGetField(in,0,"y"));
		if (x == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "vertices.x is NULL");
		if (y == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "vertices.y is NULL");
	}
};

struct Edges {
	uint32_t nie;
	uint32_t nbe;
	uint32_t ne;
	uint32_t *type;
    uint32_t *v1;
    uint32_t *v2;
    double *length;
    double *nx;
    double *ny;
    uint32_t *cp;
    uint32_t *cm;
	Edges(const mxArray *in)
	{
		nie = (uint32_t)mxGetScalar(mxGetField(in,0,"nie"));
		nbe = (uint32_t)mxGetScalar(mxGetField(in,0,"nbe"));
		ne  = (uint32_t)mxGetScalar(mxGetField(in,0,"ne"));
		type = mxGetUint32s(mxGetField(in,0,"type"));
		v1 = mxGetUint32s(mxGetField(in,0,"v1"));
		v2 = mxGetUint32s(mxGetField(in,0,"v2"));
		length = mxGetDoubles(mxGetField(in,0,"length"));
		nx = mxGetDoubles(mxGetField(in,0,"nx"));
		ny = mxGetDoubles(mxGetField(in,0,"ny"));
		cp = mxGetUint32s(mxGetField(in,0,"cp"));
		cm = mxGetUint32s(mxGetField(in,0,"cm"));
		if (type == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.type is NULL");
		if (v1 == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.v1 is NULL");
		if (v2 == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.v2 is NULL");
		if (length == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.length is NULL");
		if (nx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.nx is NULL");
		if (ny == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.ny is NULL");
		if (cp == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.cp is NULL");
		if (cm == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.cm is NULL");
	}
};

struct Cells {
    uint32_t nc;
    double *cx;
    double *cy;
    double *area;
    double *perimeter;
    double *h;
    uint32_t mne;
    uint32_t *ne;
    int32_t *e;
    uint32_t *nac;
    uint32_t *ac;
    Cells(const mxArray *in)
    {
    	nc = (uint32_t)mxGetScalar(mxGetField(in,0,"nc"));
    	cx = mxGetDoubles(mxGetField(in,0,"cx"));
    	cy = mxGetDoubles(mxGetField(in,0,"cy"));
    	area = mxGetDoubles(mxGetField(in,0,"area"));
    	perimeter = mxGetDoubles(mxGetField(in,0,"perimeter"));
    	h = mxGetDoubles(mxGetField(in,0,"h"));
    	mne = (uint32_t)mxGetScalar(mxGetField(in,0,"mne"));
    	ne = mxGetUint32s(mxGetField(in,0,"ne"));
    	e = mxGetInt32s(mxGetField(in,0,"e"));
    	nac = mxGetUint32s(mxGetField(in,0,"nac"));
    	ac = mxGetUint32s(mxGetField(in,0,"ac"));
    	if (cx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.cx = NULL");
    	if (cy == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.cy = NULL");
    	if (area == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.area = NULL");
    	if (perimeter == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.perimeter = NULL");
    	if (h == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.h = NULL");
    	if (ne == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.ne = NULL");
    	if (e == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.e = NULL");
    	if (nac == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.nac = NULL");
    	if (ac == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.ac = NULL");
    }
};
