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
	// Geometric attributes
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

    // FVM attributes
    uint32_t nq;
    double *qx;
    double *qw;
    double *up;
    double *um;
    double *tnf;
    double *mws;

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
		nq = (uint32_t)mxGetScalar(mxGetField(in,0,"nq"));
		qx = mxGetDoubles(mxGetField(in,0,"qx"));
		qw = mxGetDoubles(mxGetField(in,0,"qw"));
		up = mxGetDoubles(mxGetField(in,0,"up"));
		um = mxGetDoubles(mxGetField(in,0,"um"));
		tnf = mxGetDoubles(mxGetField(in,0,"tnf"));
		mws = mxGetDoubles(mxGetField(in,0,"mws"));
		if (qx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.qx is NULL");
		if (qw == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.qw is NULL");
		// if (up == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.up is NULL");
		// if (um == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.um is NULL");
		// if (tnf == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.tnf is NULL");
		// if (mws == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.mws is NULL");
	}
};

struct Cells {
	// Geometric attributes
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

    // FVM attributes
    uint32_t nu;
    double *u;
    double *mws;
    double *camb;

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
    	nu = (uint32_t)mxGetScalar(mxGetField(in,0,"nu"));
    	u = mxGetDoubles(mxGetField(in,0,"u"));
    	mws = mxGetDoubles(mxGetField(in,0,"mws"));
    	camb = mxGetDoubles(mxGetField(in,0,"camb"));
    }
};

struct Vec2 {
	double x;
	double y;	
};

Vec2 edge_lerp(double t, Vertices &vertices, Edges &edges, uint32_t j)
{
	uint32_t v1 = edges.v1[j-1];
	uint32_t v2 = edges.v2[j-1];
	double v1x = vertices.x[v1-1];
	double v2x = vertices.x[v2-1];
	double v1y = vertices.y[v1-1];
	double v2y = vertices.y[v2-1];
	return {(1.0-t)*v1x + t*v2x,
		    (1.0-t)*v1y + t*v2y};
}
