#pragma once
#include <cstdint>
#include "mex.h"
#include "matrix.h"
#include "../Geometry/polymesh.h"

struct VerticesFVM: Vertices {
	VerticesFVM(const mxArray *in): Vertices(in) {}
};

struct EdgesFVM: Edges {

    // FVM attributes
    uint32_t nq;
    double *qx;
    double *qw;
    double *up;
    double *um;
    double *tnf;
    double *mws;

	EdgesFVM(const mxArray *in): Edges(in)
	{
		mxArray *nq_field = mxGetField(in,0,"nq");
		if (nq_field == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges has no field nq");
		nq = (uint32_t)mxGetScalar(nq_field);
		qx = mxGetDoubles(mxGetField(in,0,"qx"));
		if (qx == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.qx is NULL");
		qw = mxGetDoubles(mxGetField(in,0,"qw"));
		if (qw == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "edges.qw is NULL");
		up = mxGetDoubles(mxGetField(in,0,"up"));
		um = mxGetDoubles(mxGetField(in,0,"um"));
		tnf = mxGetDoubles(mxGetField(in,0,"tnf"));
		mws = mxGetDoubles(mxGetField(in,0,"mws"));
	}
};

struct CellsFVM: Cells {

    // FVM attributes
    uint32_t nu;
    double *u;
    double *mws;
    double *camb;

    CellsFVM(const mxArray *in): Cells(in)
    {
    	mxArray *nu_field = mxGetField(in,0,"nu");
    	if (nu_field == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells has no field nu");
    	nu = (uint32_t)mxGetScalar(nu_field);
    	u = mxGetDoubles(mxGetField(in,0,"u"));
    	if (u == nullptr) mexErrMsgIdAndTxt("MEX:nullptr", "cells.u is NULL");
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
