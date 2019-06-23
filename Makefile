CFLAGS = -Wall
FFLAGS = -Wall
LDFLAGS = -O3

APPS = ddt_bw col_matrix_lat_indexed col_matrix_lat_vector

all: ${APPS}

clean:
	rm -f *.o ${APPS}

ddt_bw: ddt_bw.f90
	mpif90 ${FFLAGS} ${LDFLAGS} ddt_bw.f90 -o ddt_bw

col_matrix_lat_indexed: col_matrix_lat_indexed.c
	mpicc ${CFLAGS} ${LDFLAGS} col_matrix_lat_indexed.c -o col_matrix_lat_indexed

col_matrix_lat_vector: col_matrix_lat_vector.c
	mpicc ${CFLAGS} ${LDFLAGS} col_matrix_lat_vector.c -o col_matrix_lat_vector

