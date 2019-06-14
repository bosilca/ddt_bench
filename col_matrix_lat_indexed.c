#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <mpi.h>

#define TAG_INDEXED 999

#if 0
#define MAX_ELEM 10
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
 * code that sends a matrix column
 * and measures the time spent
 *
 * The part that is sent is represented below:
 *
 *                              V
 [00][00] ................. [00][CC] .................. [00][99] 
 [01][00] ................. [01][CC] .................. [01][99]
 [02][00] ................. [02][CC] .................. [02][99]
...
 [97][00] ................. [97][CC] .................. [97][99]
 [98][00] ................. [98][CC] .................. [98][99]
 [99][00] ................. [99][CC] .................. [99][99]
 */


double snd_matrix_double[MAX_ELEM][MAX_ELEM];
double rcv_matrix_double[MAX_ELEM][MAX_ELEM];
double ref_rcv_matrix_double[MAX_ELEM][MAX_ELEM];

int snd_matrix_int[MAX_ELEM][MAX_ELEM];
int rcv_matrix_int[MAX_ELEM][MAX_ELEM];
int ref_rcv_matrix_int[MAX_ELEM][MAX_ELEM];

void init_matrices(int opt, int col, int blocklen)
{
    int i, j;

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
        for (i = 0; i < MAX_ELEM; i++) {
            for( j = 0; j < blocklen; j++ ) {
                ref_rcv_matrix_double[i][col+j] = 1.0 * i * MAX_ELEM + col + j;
            }
        }
    } else {
        for (i = 0; i < MAX_ELEM; i++) {
            for (j = 0; j < MAX_ELEM; j++) {
                snd_matrix_int[i][j] = i * MAX_ELEM + j;
                rcv_matrix_int[i][j] = -1;
                ref_rcv_matrix_int[i][j] = -1;
            }
        }
        for (i = 0; i < MAX_ELEM; i++) {
            for( j = 0; j < blocklen; j++ ) {
                ref_rcv_matrix_int[i][col+j] = i * MAX_ELEM + col + j;
            }
        }
    }
}


void warmup(int my_rank, int col, int blocklen, MPI_Datatype dtt, int opt)
{
    MPI_Status status;
    int i, j;
    int errors;
    void *ptr_send, *ptr_recv;

    if (!NB_WARMUP) return;

    if (opt == OPT_DOUBLE) {
        ptr_send = snd_matrix_double;
        ptr_recv = rcv_matrix_double;
    } else {
        ptr_send = snd_matrix_int;
        ptr_recv = rcv_matrix_int;
    }

    if (!my_rank) {
        for (i = 0; i < NB_WARMUP; i++) {
            MPI_Send(ptr_send, 1, dtt, 1, TAG_INDEXED, MPI_COMM_WORLD);
            MPI_Recv(ptr_recv, 1, dtt, 1, TAG_INDEXED, MPI_COMM_WORLD, &status);
        }
    } else {
        for (i = 0; i < NB_WARMUP; i++) {
            MPI_Recv(ptr_recv, 1, dtt, 0, TAG_INDEXED, MPI_COMM_WORLD, &status);
            MPI_Send(ptr_send, 1, dtt, 0, TAG_INDEXED, MPI_COMM_WORLD);
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
    int i, col, opt, blocklen = 1;
    int size, rank;
    char *type;
    void *ptr_send, *ptr_recv;
    MPI_Datatype newtype;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    //sleep(20);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size != 2) {
        if (!rank) {
            printf("shoulf be run with 2 ranks!!!!\n");
        }
        MPI_Finalize();
        exit(1);
    }
    if (argc < 3) {
        if (!rank) {
            printf("Syntax: %s <column number> {double|int} [blocklen=1]\n", argv[0]);
        }
        MPI_Finalize();
        exit(1);
    }
    col = atoi(argv[1]);
    if (col < 0 || col >= MAX_ELEM) {
        if (!rank) {
            printf("%s: <column number> should be >= 0 and < %d\n", argv[0], MAX_ELEM);
        }
        MPI_Finalize();
        exit(1);
    }
    type = argv[2];
    if (strcmp(type, "double")) {
        if (strcmp(type, "int")) {
            if (!rank) {
                printf("Syntax: %s <column number> {double|int} [blocklen=1]\n", argv[0]);
            }
            MPI_Finalize();
            exit(1);
        } else {
            opt = OPT_INT;
            ptr_send = snd_matrix_int;
            ptr_recv = rcv_matrix_int;
        }
    } else {
        opt = OPT_DOUBLE;
        ptr_send = snd_matrix_double;
        ptr_recv = rcv_matrix_double;
    }
    if( argc > 3 ) {
        blocklen = atoi(argv[3]);
        if( (blocklen < 1) || ((col + blocklen) > MAX_ELEM) ) {
            if(!rank) {
                printf("When col = %d blocklen must be [1 .. %d[\n", col, MAX_ELEM-col);
                MPI_Finalize();
                exit(1);
            }
        }
    }

    /* Compute start and size of each row.
     * Actually size is constant since we are transferring a column. */
    for (i = 0; i < MAX_ELEM; i++) {
        disps[i] = i * MAX_ELEM + col;
        blklens[i] = blocklen;
    }

    /*
     * count    = # of blocks
     * array_of_blocklen = # of elems of oldtype in each block
     * array_of_displacements   = displacement of each block (in # of elems) wrt beginning
     */
    MPI_Type_indexed(MAX_ELEM, blklens, disps,
                     (opt == OPT_INT) ? MPI_INT : MPI_DOUBLE, &newtype);
    MPI_Type_commit(&newtype);

    init_matrices(opt, col, blocklen);

    warmup(rank, col, blocklen, newtype, opt);

    if (rank == 0) {

        t_beg = MPI_Wtime();
        for (i = 0; i < NB_LOOPS; i++) {
            MPI_Send(ptr_send, 1, newtype, 1, TAG_INDEXED, MPI_COMM_WORLD);
            MPI_Recv(ptr_recv, 1, newtype, 1, TAG_INDEXED, MPI_COMM_WORLD, &status);
        }
        t_end = MPI_Wtime();

    } else {

        for (i = 0; i < NB_LOOPS; i++) {
            MPI_Recv(ptr_recv, 1, newtype, 0, TAG_INDEXED, MPI_COMM_WORLD, &status);
            MPI_Send(ptr_send, 1, newtype, 0, TAG_INDEXED, MPI_COMM_WORLD);
        }

    }

    if (rank == 0) {
        latency = (t_end - t_beg) * 1e6 / (2.0 * NB_LOOPS);

        fprintf(stdout, "(%s) LATENCY MATRIX %d x %d %ss (usecs)\n",
                argv[0], MAX_ELEM, MAX_ELEM,
                (opt == OPT_INT) ? "int" : "double");
        fprintf(stdout, "col=%d blocklen=%d -----------------> %f usecs\n", col, blocklen, latency);
        fflush(stdout);
    }


    MPI_Finalize();
    exit(0);
}


