/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Kernel 1: CRS construction
// Simple two-pass CRS construction using Active Messages

#include "common.h"
#include "csr_reference.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <search.h>

int64_t nverts_known = 0;
int *degrees;
int64_t *column;
float *weights;
extern oned_csr_graph g; //from bfs_reference for isisolated function

//this function is needed for roots generation
int isisolated(int64_t v) {
	return (g.rowstarts[v]==g.rowstarts[v+1]);
}

void halfedgehndl(int data)
{  degrees[data]++; }

void fulledgehndl(int *data) {
	int vloc = *data;
	int64_t gtgt = *((int64_t*)(data+1));
	SETCOLUMN(degrees[vloc]++,gtgt);
#ifdef SSSP
	float w = ((float*)data)[3];
	weights[degrees[vloc]-1]=w;
#endif
}

void send_half_edge (int64_t src,int64_t tgt) {
	int vloc=src;
    halfedgehndl(vloc);
	if(tgt>=nverts_known) nverts_known=tgt+1;
}
#ifdef SSSP
void send_full_edge (int64_t src,int64_t tgt,float w) {
	int vloc[4];
	vloc[0]=src;
	memcpy(vloc+1,&tgt,8);
	memcpy(vloc+3,&w,4);
    fulledgehndl(vloc);
}
#else
void send_full_edge (int64_t src,int64_t tgt) {
	int vloc[3];
	vloc[0]=src;
	memcpy(vloc+1,&tgt,8);
    fulledgehndl(vloc);
}
#endif

void convert_graph_to_oned_csr(const tuple_graph* const tg, oned_csr_graph* const g) {
	g->tg = tg;

	size_t i,j,k;

	int64_t nvert=tg->nglobaledges/2;
	nvert/=1;
	nvert+=1;
	degrees=xcalloc(nvert,sizeof(int));

	int numiters=ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg);
	// First pass : calculate degrees of each vertex
	ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
		ptrdiff_t j;
		for (j = 0; j < bufsize; ++j) {
			int64_t v0 = get_v0_from_edge(&buf[j]);
			int64_t v1 = get_v1_from_edge(&buf[j]);
			if(v0==v1) continue;
			send_half_edge(v0, v1);
			send_half_edge(v1, v0);
		}
	} ITERATE_TUPLE_GRAPH_END;

	int64_t nglobalverts = 0;
	nglobalverts=nverts_known+1;
	g->nglobalverts = nglobalverts;
	size_t nlocalverts = nglobalverts;
	g->nlocalverts = nlocalverts;

	//graph stats printing
#ifdef DEBUGSTATS
	long maxdeg=0,isolated=0,totaledges=0,originaledges;
	long maxlocaledges,minlocaledges;
	for(i=0;i<g->nlocalverts;i++) {
		long deg = degrees[i];
		totaledges+=deg;
		if(maxdeg<deg) maxdeg=deg;
		if(!deg) isolated++;
	}
	originaledges=totaledges;
	maxlocaledges=totaledges;
	minlocaledges=totaledges;
	long averageedges = totaledges/num_pes();
	double disbalance = (double)(maxlocaledges-minlocaledges)/(double)averageedges * 100.0;
	if(!my_pe()) printf("\n maxdeg %lld verts %lld, isolated %lld edges %lld\n\t A max %ld min %ld ave %ld delta %ld percent %3.2f\n ",
			maxdeg,g->nglobalverts,isolated,totaledges,maxlocaledges,minlocaledges,averageedges,maxlocaledges-minlocaledges,disbalance);

	// finished stats printing

	g->notisolated=g->nglobalverts-isolated;
#endif
	unsigned int *rowstarts = xmalloc((nlocalverts + 1) * sizeof(int));
	g->rowstarts = (int *) rowstarts;

	rowstarts[0] = 0;
	for (i = 0; i < nlocalverts; ++i) {
		rowstarts[i + 1] = rowstarts[i] + (i >= nlocalverts ? 0 : degrees[i]);
		degrees[i] = rowstarts[i];
	}

	size_t nlocaledges = rowstarts[nlocalverts];
	g->nlocaledges = nlocaledges;

	int64_t colalloc = BYTES_PER_VERTEX*nlocaledges;
	colalloc += (4095);
	colalloc /= 4096;
	colalloc *= 4096;
	column = xmalloc(colalloc);
#ifdef SSSP
	weights = xmalloc(4*nlocaledges);
	g->weights = weights;
#endif
	//long allocatededges=colalloc;
	g->column = column;

	//Next pass , actual data transfer: placing edges to its places in column and hcolumn
	ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
		ptrdiff_t j;
		for (j = 0; j < bufsize; ++j) {
			int64_t v0 = get_v0_from_edge(&buf[j]);
			int64_t v1 = get_v1_from_edge(&buf[j]);
			if(v0==v1) continue;
#ifdef SSSP
			send_full_edge(v0, v1,wbuf[j]);
			send_full_edge(v1, v0,wbuf[j]);
#else
			send_full_edge(v0, v1);
			send_full_edge(v1, v0);
#endif
		}
	} ITERATE_TUPLE_GRAPH_END;

	free(degrees);
}

void free_oned_csr_graph(oned_csr_graph* const g) {
	if (g->rowstarts != NULL) {free(g->rowstarts); g->rowstarts = NULL;}
	if (g->column != NULL) {free(g->column); g->column = NULL;}
}
