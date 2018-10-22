MPI_CFLAGS := $(shell mpicc -showme:compiler)
MPI_LDFLAGS := $(shell mpicc -showme:linker)

CFLAGS = -Drestrict=__restrict__ -O3 -DGRAPH_GENERATOR_MPI -DREUSE_CSR_FOR_VALIDATION -I../aml $(MPI_CFLAGS)
LDFLAGS = -lpthread $(MPI_LDFLAGS)

CC = gcc

all: graph500_sssp

GENERATOR_SOURCES = ../generator/graph_generator.c ../generator/make_graph.c ../generator/splittable_mrg.c ../generator/utils.c
SOURCES = main.c utils.c validate.c ../aml/aml.c
HEADERS = common.h csr_reference.h bitmap_reference.h

graph500_sssp: bfs_reference.c $(SOURCES) $(HEADERS) $(GENERATOR_SOURCES) csr_reference.c sssp_reference.c
	$(CC) $(CFLAGS) -DSSSP -DSKIP_BFS $(LDFLAGS) -o graph500_sssp bfs_reference.c sssp_reference.c csr_reference.c $(SOURCES) $(GENERATOR_SOURCES) -lm

clean:
	-rm -f graph500_*
