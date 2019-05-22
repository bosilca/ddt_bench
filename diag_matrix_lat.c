#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

#define TAG_INDEXED 999

#if 0
#define MAX_ELEM  10
#define NB_WARMUP 0
#define NB_LOOPS  1
#else

#define MAX_ELEM  4096
#define NB_WARMUP 1000
#define NB_LOOPS  10000

#endif

#define OPT_INT    0
#define OPT_DOUBLE 1

/*
 * code that sends a matrix diagonal and measures the time spent
 *
 * The part that is sent is represented below:
 *
 *                               
 [00][00] ...................................................... 
 ........ [01][01] .............................................
 ................. [02][02] ....................................
...
 .................................... [97][97] .................
 ............................................. [98][98] ........
 ...................................................... [99][99]
 */


double snd_matrix_double[MAX_ELEM][MAX_ELEM];
double rcv_matrix_double[MAX_ELEM][MAX_ELEM];
double ref_rcv_matrix_double[MAX_ELEM][MAX_ELEM];

int snd_matrix_int[MAX_ELEM][MAX_ELEM];
int rcv_matrix_int[MAX_ELEM][MAX_ELEM];
int ref_rcv_matrix_int[MAX_ELEM][MAX_ELEM];

void init_matrices(int opt, int diag)
{
    int i, j;
    int min_idx, max_idx;

    if (diag >= 0) {
        min_idx = 0;
        max_idx = MAX_ELEM - diag;
    } else {
        min_idx = -diag;
        max_idx = MAX_ELEM;
    }

    /*
     * Init the matrices
     */
    if (opt == OPT_DOUBLE) {
        for (i = 0; i < MAX_ELEM; i++) {
            for (j = 0; j < MAX_ELEM; j++) {
                snd_matrix_double[i][j] = 1.0 * i * MAX_ELEM + j;
                rcv_matrix_double[i][j] = -1.0;
                ref_rcv_matrix_double[i][j] = -1.0;
            }
        }
        for (i = min_idx; i < max_idx; i++) {
            ref_rcv_matrix_double[i][i + diag] = 1.0 * i * MAX_ELEM + i + diag;
        }
    } else {
        for (i = 0; i < MAX_ELEM; i++) {
            for (j = 0; j < MAX_ELEM; j++) {
                snd_matrix_int[i][j] = i * MAX_ELEM + j;
                rcv_matrix_int[i][j] = -1;
                ref_rcv_matrix_int[i][j] = -1;
            }
        }
        for (i = min_idx; i < max_idx; i++) {
            ref_rcv_matrix_int[i][i + diag] = i * MAX_ELEM + i + diag;
        }
    }
}


void warmup(int my_rank, MPI_Datatype dtt, int opt)
{
    MPI_Status status;
    int i, j;
    int errors;
    void *snd_matrix, *rcv_matrix;

    if (opt == OPT_DOUBLE) {
        snd_matrix = snd_matrix_double;
	rcv_matrix = rcv_matrix_double;
    } else {
        snd_matrix = snd_matrix_int;
	rcv_matrix = rcv_matrix_int;
    }

    if (!my_rank) {
        for (i = 0; i < NB_WARMUP; i++) {
            MPI_Send(snd_matrix, 1, dtt, 1, TAG_INDEXED, MPI_COMM_WORLD);
            MPI_Recv(rcv_matrix, 1, dtt, 1, TAG_INDEXED, MPI_COMM_WORLD, &status);
        }
    } else {
        for (i = 0; i < NB_WARMUP; i++) {
            MPI_Recv(rcv_matrix, 1, dtt, 0, TAG_INDEXED, MPI_COMM_WORLD, &status);
            MPI_Send(snd_matrix, 1, dtt, 0, TAG_INDEXED, MPI_COMM_WORLD);
        }

        /* Check for errors */
        errors = 0;
        for (i = 0; i < MAX_ELEM; i++) {
            for (j = 0; j < MAX_ELEM; j++) {
                if (opt == OPT_DOUBLE) {
                    if (rcv_matrix_double[i][j] != ref_rcv_matrix_double[i][j]) {
                            printf("!!!!!!!!! rcv_matrix[%d][%d] : expected %f GOT %f\n", i, j, ref_rcv_matrix_double[i][j], rcv_matrix_double[i][j]);
                        errors++;
                    }
                } else {
                    if (rcv_matrix_int[i][j] != ref_rcv_matrix_int[i][j]) {
                        printf("!!!!!!!!! rcv_matrix[%d][%d] : expected %d GOT %d\n", i, j, ref_rcv_matrix_int[i][j], rcv_matrix_int[i][j]);
                    }
                }
            }
        }
        if (errors) {
            printf("FOUND %d errors", errors);
            printf("    )-:\n");
        }
#if 0
        else        { printf("w000t!!!\n");   }
#endif
    }
}


int main(int argc, char **argv)
{
    double t_beg, t_end;
    double latency;
    int disps[MAX_ELEM], blklens[MAX_ELEM];
    int i, diag, opt;
    int min_idx, max_idx;
    int size, rank;
    char *type;
    void *snd_matrix, *rcv_matrix;
    MPI_Datatype newtype;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size != 2) {
        if (!rank) {
            printf("shoulf be run with 2 ranks!!!!\n");
        }
        MPI_Finalize();
        exit(1);
    }
    if (argc != 3) {
        if (!rank) {
            printf("Syntax: %s <diagonal number {double|int}>\n", argv[0]);
            printf("           0: actual diagonal\n");
            printf("           x (x>0): diagonal from [0][x] to [%d - x][%d]\n", MAX_ELEM - 1, MAX_ELEM - 1);
            printf("           x (x<0): diagonal from [x][0] to [%d][%d - x]\n", MAX_ELEM - 1, MAX_ELEM - 1);
        }
        MPI_Finalize();
        exit(1);
    }
    diag = atoi(argv[1]);
    if (diag <= -MAX_ELEM || diag >= MAX_ELEM) {
        if (!rank) {
            printf("%s: <diagonal number> should be > %d and < %d\n", argv[0], -MAX_ELEM, MAX_ELEM);
        }
        MPI_Finalize();
        exit(1);
    }
    type = argv[2];
    if (strcmp(type, "double")) {
        if (strcmp(type, "int")) {
            if (!rank) {
                printf("Syntax: %s <diagonal number> {double|int}\n", argv[0]);
            }
            MPI_Finalize();
            exit(1);
        } else {
            opt = OPT_INT;
            snd_matrix = snd_matrix_int;
	    rcv_matrix = rcv_matrix_int;
        }
    } else {
        opt = OPT_DOUBLE;
        snd_matrix = snd_matrix_double;
	rcv_matrix = rcv_matrix_double;
    }

    /* Compute start and size of each row.
     * Actually size is constant since we are transferring a diagonal.
     *
     *     | Starting from this diagonal, we have 0 displacements and blklens
     *     |          at the end of the arrays
     *     v 
     *   X D X X X X X X X
     *   d X D X X X X X X
     *   X d X D X X X X X
     *   X X d X D X X X X
     *   X X X d X D X X X
     *   X X X X d X D X X
     *   X X X X X d X D X
     *   X X X X X X d X D
     *   X X X X X X X d X
     *                 ^
     *                 | Starting from this diagonal, we have 0 displacements and
     *                 |          blklens at the beginning of the arrays
     */

    for (i = 0 ; i < MAX_ELEM; i++) {
        disps[i] = 0;
        blklens[i] = 0;
    }
    if (diag >= 0) {
        min_idx = 0;
        max_idx = MAX_ELEM - diag;
    } else {
        min_idx = -diag;
        max_idx = MAX_ELEM;

    }
    for (i = min_idx; i < max_idx; i++) {
        disps[i] = i * MAX_ELEM + i + diag;
        blklens[i] = 1;
    }

    /*
     * count    = # of blocks
     * array_of_blocklen = # of elems of oldtype in each block
     * array_of_displacements   = displacement of each block (in # of elems) wrt beginning
     */
    if (opt == OPT_INT) {
        MPI_Type_indexed(MAX_ELEM, blklens, disps, MPI_INT, &newtype);
    } else {
        MPI_Type_indexed(MAX_ELEM, blklens, disps, MPI_DOUBLE, &newtype);
    }
    MPI_Type_commit(&newtype);

    init_matrices(opt, diag);

    warmup(rank, newtype, opt);

    if (rank == 0) {

        t_beg = MPI_Wtime();
        for (i = 0; i < NB_LOOPS; i++) {
            MPI_Send(snd_matrix, 1, newtype, 1, TAG_INDEXED, MPI_COMM_WORLD);
            MPI_Recv(rcv_matrix, 1, newtype, 1, TAG_INDEXED, MPI_COMM_WORLD, &status);
        }
        t_end = MPI_Wtime();

    } else {

        for (i = 0; i < NB_LOOPS; i++) {
            MPI_Recv(rcv_matrix, 1, newtype, 0, TAG_INDEXED, MPI_COMM_WORLD, &status);
            MPI_Send(snd_matrix, 1, newtype, 0, TAG_INDEXED, MPI_COMM_WORLD);
        }

    }

    if (rank == 0) {
        latency = (t_end - t_beg) * 1e6 / (2.0 * NB_LOOPS);

        fprintf(stdout, "(%s) LATENCY MATRIX %d x %d %ss (usecs)\n",
                argv[0], MAX_ELEM, MAX_ELEM,
                (opt == OPT_INT) ? "int" : "double");
        fprintf(stdout, "diag=%d -----------------> %f usecs\n", diag, latency);
        fflush(stdout);
    }


    MPI_Finalize();
    exit(0);
}


