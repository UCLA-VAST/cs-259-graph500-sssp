/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:                Anton Korzh anton@korzh.us
                                Graph500 Steering Committee
                                http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Validation kernel for both BFS and SSSP
//  1. building CRS graph representation from tuple graph,
//  2. creating distances map (for BFS)
//  3. checking all edges to follow the triangle rule
//  4. checking all claimed vertices actually reached with correct distance
//  5. validating number of edges visited

#include "common.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <search.h>
#include <float.h>

#ifdef REUSE_CSR_FOR_VALIDATION
#include "csr_reference.h"
extern int64_t* column;
extern int* rowstarts;
#ifdef SSSP
extern float* weights;
#endif
#else
#define SETCOLUMN(a,b) vcolumn[a]=b
#define COLUMN(i) vcolumn[i]
#endif
int firstvalidationrun=1;
//int failedttovalidate=0;
int validatingbfs=0;

unsigned int *vdegrees,*vrowstarts;
int64_t *vcolumn;
#ifdef SSSP
float* vweights;
#endif
int64_t *globpred,nedges_traversed;
float *globdist,prevlevel;
int64_t val_errors=0;
int failedtovalidate=0;
int64_t newvisits;
int * confirmed=NULL;
int64_t maxvertex;

void frompredhndl(int* data) {
	int vfrom = *data;
	int64_t predfrom = vfrom;
	int vloc = *(data+1);

	if(globpred[vloc] == predfrom && globdist[vloc]==FLT_MAX) globdist[vloc]=prevlevel+1.0,newvisits++;
}

void send_frompred (int vfrom,int64_t src) {
	int vloc[2]={vfrom, src};
    frompredhndl(vloc);
}

void vhalfedgehndl(int data)
{  vdegrees[data]++; }

void send_half (int64_t src) {
	int vloc=src;
    vhalfedgehndl(vloc);
}

void vfulledgehndl(int * data) {
	int vloc = *data;
	int64_t gtgt = *((int64_t*)(data+1));
	int next = vdegrees[vloc]++;
	SETCOLUMN(next,gtgt);
#ifdef SSSP
	vweights[next] = ((float*)data)[3];
#endif
}

void vsend_full_edge (int64_t src,int64_t tgt,float w) {
	int vloc[4];
	vloc[0]=src;
	memcpy(vloc+1,&tgt,8);
#ifdef SSSP
	memcpy(vloc+3,&w,4);
    vfulledgehndl(vloc);
#else
    vfulledgehndl(vloc);
#endif
}

typedef struct edgedist {
	unsigned int vfrom;
	unsigned int vloc;
	int64_t predfrom;
	float distfrom;
#ifdef SSSP
	float w;
#endif
} edgedist;

#define DUMPERROR(text) { printf("Validation Error: %s, edge %llu %llu weight %f pred0 %llu pred1 %llu dist0 %f dist1 %f\n",text,v0,v1,w,predv0,predv1,distv0,distv1); val_errors++; return; }

//main validation handler: tracks all edges and at delivery has both vertex preds and distances to be checked
void edgepreddisthndl(edgedist* m) {
	unsigned int v1loc = m->vloc;

	int64_t v1 = v1loc;
	int64_t v0 = m->vfrom;
	int64_t predv1 = globpred[v1loc];
	int64_t predv0 = m->predfrom;
#ifdef SSSP
	float w = (validatingbfs?1.0:m->w);
#else
	float w = 1.0;
#endif
	float distv0 = m->distfrom;
	float distv1 = globdist[v1loc];
	if(predv0==-1 && predv1==-1) return; else if(v0<v1) nedges_traversed++;

	if((predv0==-1 && predv1!=-1) || (predv1==-1 && predv0!=-1)) DUMPERROR("edge connecting visited and unvisited vertices");
	if(predv1==v0 && distv1 == distv0+w) confirmed[v1loc]=1; //confirm pred/dist as existing edge

	if(distv0+w < distv1 || distv1+w<distv0) DUMPERROR("triangle rule violated");
}

void sendedgepreddist(unsigned vloc,unsigned int vedge) {
	edgedist m = {vloc,COLUMN(vedge),globpred[vloc],globdist[vloc]
#ifdef SSSP
		,vweights[vedge]
#endif
	};
    edgepreddisthndl(&m);
}

void makedepthmapforbfs(const size_t nlocalverts,const int64_t root,int64_t * const pred,float* dist) {

	int i,j;
	for(i=0;i<nlocalverts;i++) {
		dist[i]=FLT_MAX; //at the end there should be no FLT_MAX left
		if(pred[i]==-1) dist[i]=-1.0;
		if(pred[i]==root) dist[i]=1.0;
	}

	//fix root distance
	dist[root]=0.0;

	newvisits=1;
	prevlevel=0.0;

	while(newvisits!=0) {
		newvisits=0;
		prevlevel+=1.0;

		for(i=0;i<nlocalverts;i++)
			if(dist[i]==prevlevel)
				for(j=vrowstarts[i];j<vrowstarts[i+1];j++)
					send_frompred(i,COLUMN(j));
	}
}

int validate_result(int isbfs,const tuple_graph* const tg, const size_t nlocalverts, const int64_t root, int64_t* const pred, float* dist,int64_t *nedges_in) {
	validatingbfs=isbfs;

	//if(failedtovalidate) return 0; //failed to allocate lots of memory for validation: skipping all validation now

	size_t i,j,k;
	if(firstvalidationrun) {
		firstvalidationrun=0;
		confirmed = xmalloc(nlocalverts*sizeof(int));
#ifdef REUSE_CSR_FOR_VALIDATION
vrowstarts=(unsigned int *)rowstarts;
#ifdef SSSP
vweights=weights;
#endif
#else
		vdegrees=xcalloc(nlocalverts,sizeof(int));

		int numiters=ITERATE_TUPLE_GRAPH_BLOCK_COUNT(tg);
		// First pass : calculate degrees of each vertex
		ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
			ptrdiff_t j;
			for (j = 0; j < bufsize; ++j) {
				int64_t v0 = get_v0_from_edge(&buf[j]);
				int64_t v1 = get_v1_from_edge(&buf[j]);
				if(v0==v1) continue;
				send_half(v0);
				send_half(v1);
			}
		} ITERATE_TUPLE_GRAPH_END;

		vrowstarts = xmalloc((nlocalverts + 1) * sizeof(int));

		vrowstarts[0] = 0;
		for (i = 0; i < nlocalverts; ++i) {
			vrowstarts[i + 1] = vrowstarts[i] + (i >= nlocalverts ? 0 : vdegrees[i]);
			vdegrees[i] = vrowstarts[i];
		}

		vcolumn = xmalloc(8*vrowstarts[nlocalverts]);
#ifdef SSSP
		vweights = xmalloc(4*vrowstarts[nlocalverts]);
#endif
		//Second pass , actual data transfer: placing edges to its places in vcolumn
		ITERATE_TUPLE_GRAPH_BEGIN(tg, buf, bufsize,wbuf) {
			ptrdiff_t j;
			for (j = 0; j < bufsize; ++j) {
				int64_t v0 = get_v0_from_edge(&buf[j]);
				int64_t v1 = get_v1_from_edge(&buf[j]);
				if(v0==v1) continue;
#ifdef SSSP
				vsend_full_edge(v0, v1,wbuf[j]);
				vsend_full_edge(v1, v0,wbuf[j]);
#else
				vsend_full_edge(v0, v1,1.0);
				vsend_full_edge(v1, v0,1.0);
#endif
			}
		} ITERATE_TUPLE_GRAPH_END;

		free(vdegrees); vdegrees=NULL;
#endif
	for (i = 0; i < nlocalverts; ++i)
		if(i>maxvertex) maxvertex = i;

	} //only run at first validation

	//Actual validation here:
	globpred=pred;
	globdist=dist;
	for (i = 0; i < nlocalverts; ++i)
		confirmed[i]=0;

	for (i = 0; i < nlocalverts; ++i)
		if((pred[i]!=-1 && pred[i]<0) || pred[i]>maxvertex) 
		printf("Validation Error: predecessor %llu of vertex %lu is out of range\n",pred[i],i),val_errors++;
	if(val_errors>0) return 0;

	if(validatingbfs) 
	makedepthmapforbfs(nlocalverts,root,pred,dist);

	for (i = 0; i < nlocalverts; ++i) {
		if(dist[i]!=-1.0 && dist[i]<0.0) 
		printf("Validation Error: distance/depth %3.2f of vertex %lu is out of range\n",dist[i],i),val_errors++;
		if((pred[i]==-1 && dist[i]!=-1.0) || (pred[i]!=-1 && dist[i]==-1.0)) 
		printf("Validation Error: vertex %lu has inconsistent predecessor %llu and distance %f\n",i,pred[i],dist[i]),val_errors++;
	}

    int vloc = root;
    if(pred[vloc]!=root || dist[vloc]!=0.0) 
    printf("Validation Error: root vertex %llu has predecessor %llu and distance %f\n",root,pred[vloc],dist[vloc]),val_errors++;
    else confirmed[vloc]=1;

	nedges_traversed=0;

	for (i = 0; i < nlocalverts; ++i)
		for(j = vrowstarts[i];j<vrowstarts[i+1];j++)
			sendedgepreddist(i,j);

	for (i = 0; i < nlocalverts; ++i)
		if(confirmed[i]==0 && pred[i]!=-1) 
		printf("Validation Error: path to vertex %lu not confirmed from predecessor %llu with distance %f\n",i,pred[i],dist[i]),val_errors++;

	if(nedges_in!=NULL && *nedges_in!=nedges_traversed) 
	printf("Validation Error: wrong nedge_traversed %llu (correct number is %llu)\n",*nedges_in,nedges_traversed),val_errors++;
	if(val_errors>0) return 0; else return 1;
}
