/* Copyright (c) 2011-2017 Graph500 Steering Committee
   All rights reserved.
   Developed by:		Anton Korzh anton@korzh.us
				Graph500 Steering Committee
				http://www.graph500.org
   New code under University of Illinois/NCSA Open Source License
   see license.txt or https://opensource.org/licenses/NCSA
*/

// Graph500: Kernel 2: BFS
// Simple level-synchronized BFS with visits as Active Messages

#include "common.h"
#include "csr_reference.h"
#include "bitmap_reference.h"
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>

// two arrays holding visited vertices for current and next level
// we swap pointers each time
int *q1,*q2;
int qc,q2c; //pointer to first free element

//VISITED bitmap parameters
unsigned long *visited;
int64_t visited_size;

//global variables of CSR graph to be used inside of AM-handlers
int64_t *column;
int64_t *pred_glob;
int * rowstarts;

oned_csr_graph g;

typedef struct visitmsg {
	int vloc;
	int vfrom;
} visitmsg;

//AM-handler for check&visit
void visithndl(const visitmsg *m) {
	if (!TEST_VISITEDLOC(m->vloc)) {
		SET_VISITEDLOC(m->vloc);
		q2[q2c++] = m->vloc;
		pred_glob[m->vloc] = m->vfrom;
	}
}

inline void send_visit(int64_t glob, int from) {
	visitmsg m = {glob,from};
    visithndl(&m);
}

void make_graph_data_structure(const tuple_graph* const tg) {
	int i,j,k;
	convert_graph_to_oned_csr(tg, &g);
	column=g.column;
	rowstarts=g.rowstarts;

	visited_size = (g.nlocalverts + ulong_bits - 1) / ulong_bits;
	q1 = xmalloc(g.nlocalverts*sizeof(int)); //100% of vertexes
	q2 = xmalloc(g.nlocalverts*sizeof(int));
	for(i=0;i<g.nlocalverts;i++) q1[i]=0,q2[i]=0; //touch memory
	visited = xmalloc(visited_size*sizeof(unsigned long));
}

void run_bfs(int64_t root, int64_t* pred) {
	int64_t nvisited;
	long sum;
	unsigned int i,j,k,lvl=1;
	pred_glob=pred;

	CLEAN_VISITED();

	qc=0; sum=1; q2c=0;

	nvisited=1;
    pred[root]=root;
    SET_VISITED(root);
    q1[0]=root;
    qc=1;

	// While there are vertices in current level
	while(sum) {
		//for all vertices in current level send visit AMs to all neighbours
		for(i=0;i<qc;i++)
			for(j=rowstarts[q1[i]];j<rowstarts[q1[i]+1];j++)
				send_visit(COLUMN(j),q1[i]);

		qc=q2c;int *tmp=q1;q1=q2;q2=tmp;
		sum=qc;

		nvisited+=sum;

		q2c=0;
	}

}

//we need edge count to calculate teps. Validation will check if this count is correct
void get_edge_count_for_teps(int64_t* edge_visit_count) {
	long i,j;
	long edge_count=0;
	for(i=0;i<g.nlocalverts;i++)
		if(pred_glob[i]!=-1) {
			for(j=rowstarts[i];j<rowstarts[i+1];j++)
				if(COLUMN(j)<=i)
					edge_count++;

		}

	*edge_visit_count=edge_count;
}

void clean_pred(int64_t* pred) {
	int i;
	for(i=0;i<g.nlocalverts;i++) pred[i]=-1;
}
void free_graph_data_structure(void) {
	int i; 
	free_oned_csr_graph(&g);
	free(q1); free(q2); free(visited);
}

size_t get_nlocalverts_for_pred(void) {
	return g.nlocalverts;
}
