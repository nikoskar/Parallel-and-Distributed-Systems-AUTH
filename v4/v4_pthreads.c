#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include "mmio.h"
#include <pthread.h>
#define NUM_THREADS 12

/* Parallel implementation of v4 using POSIX Threads */

/* define the struct that is going to be used in pthread_create */
typedef struct {
  uint32_t i_start;      /* Starting loop index for the current thread */
  uint32_t i_end;        /* Ending loop index for the current threads  */
  uint32_t *csc_col;     /* CSC Column array for the multiplication    */
  uint32_t *csc_row;     /* CSC row array for the multiplication       */
  uint32_t *masked_vals; /* masked matrix result array                 */
} ARRINFO;

/* A .* (A*A) implementation routine */
void *sparseSelfProduct(void *args);

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);

int main(int argc, char *argv[])
{
    FILE *f;
    double *val;
    int ret_code;
    uint32_t *I, *J;
    MM_typecode matcode;
    uint32_t M, N, nz, sum;
    long elapsed_sec, elapsed_nsec;
    struct timespec ts_start, ts_end;
		struct timeval ts_start, ts_end;

    if (argc < 2)
	  {
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	  }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }


    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    I = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));

    if (!mm_is_pattern(matcode))
    {
      for (uint32_t i=0; i<nz; i++)
      {
          fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
          I[i]--;
          J[i]--;
          I[nz + i] = J[i];
          J[nz + i] = I[i];
      }
    }
    else
    {
      for (uint32_t i=0; i<nz; i++)
      {
          fscanf(f, "%d %d\n", &I[i], &J[i]);
          val[i]=1;
          I[i]--;
          J[i]--;
          I[nz + i] = J[i];
          J[nz + i] = I[i];
      }
    }

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);

    uint32_t * csc_row = (uint32_t *)malloc(2 * nz     * sizeof(uint32_t));
    uint32_t * csc_col = (uint32_t *)malloc((N + 1) * sizeof(uint32_t));

    if (I[0] > J[0]) coo2csc(csc_row, csc_col, J, I, 2 * nz, N, 0);
    else coo2csc(csc_row, csc_col, I, J, 2 * nz, N, 0);

    uint32_t *masked_vals = malloc(2 * nz * sizeof(uint32_t));
    for(uint32_t i = 0; i < 2 * nz; i++) {
      masked_vals[i] = 0;
    }

    uint32_t *e_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *c3_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *result = malloc(N * sizeof(uint32_t));

    for(uint32_t i = 0; i < N; i++) {
      e_vector[i] = 1 ;
      result[i] = 0;
      c3_vector[i] = 0;
    }

    printf("\n*** Triangle counting has now started ***\n");

    /* In this case we need to pareallelize the algorithm with pthreads. In this
     * case the code structure depends heavily upon the number of threads. We
     * fire up threads by invokind funcitons so in this case we need to port the
     * triangle counter into a function. In order to parallelize the execution
     * we will split the for loop (practically the calculation of a specific
     * amount of non zero elements of A.*(A*A)) into smaller loops - each executed
     * by a separate thread.
     */



    /* Create threads to perform the triangle search. Also make them joinable  */
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    int rc;;
    ARRINFO arr_info[NUM_THREADS];

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
    //gettimeofday(&ts_start,NULL);

    for(uint32_t i = 0; i < NUM_THREADS ; i++) {

      arr_info[i].csc_col = csc_col;
      arr_info[i].csc_row = csc_row;
      arr_info[i].i_start = i * N / NUM_THREADS;
      arr_info[i].i_end = i * N / NUM_THREADS + N / NUM_THREADS;
      arr_info[i].masked_vals = masked_vals;
      //printf("\n***************test**********************\n");
      //printf("in main with index start: %d, end:  %d\n",arr_info[i].i_start, arr_info[i].i_end );

      rc = pthread_create(&threads[i], &attr, sparseSelfProduct, &arr_info[i]);
      if(rc) {
        printf("\nERROR: Something went wrong while creating thread %d\n", i);
        printf("\nThe error code returned by pthread_create is %d", rc);
        exit(-1);
      }
    }

    uint32_t final_left = N - (N / NUM_THREADS) * NUM_THREADS;

    for(uint32_t t = 0; t < NUM_THREADS; t++) {
       rc = pthread_join(threads[t], NULL);
       if (rc) {
          printf("ERROR; return code from pthread_join() is %d\n", rc);
          exit(-1);
      }
    }

    ARRINFO arr_info_single;
    if (final_left != 0) {
      arr_info_single.csc_col = csc_col;
      arr_info_single.csc_row = csc_row;
      arr_info_single.i_start = N - final_left;
      arr_info_single.i_end = N;
      arr_info_single.masked_vals = masked_vals;
      //arr_info[i].id = i;
      printf("\n***************test**********************\n");
      //printf("in main with index start: %d, end:  %d\n",arr_info[i].i_start, arr_info[i].i_end );

      rc = pthread_create(&threads[0], &attr, sparseSelfProduct, &arr_info_single);
      if(rc) {
        printf("\nERROR: Something went wrong while creating  the final thread\n");
        printf("\nThe error code returned by pthread_create is %d", rc);
        exit(-1);
      }
    }


    pthread_attr_destroy(&attr);

    for(uint32_t i = 0; i < N; i++)
    {
      uint32_t start = csc_col[i];
      uint32_t end = csc_col[i + 1];

      for(uint32_t j = 0; j < end - start; j++) {
        result[i] += masked_vals[start + j];

      }
    }

    uint32_t total = 0;
    for(uint32_t i = 0; i < N; i++) {
      c3_vector[i] = result[i]/2;
      total += c3_vector[i];
    }

    printf("\n*** Total triangles: %d ***\n: ", total/3);



    // gettimeofday(&ts_end,NULL);
    // double elapsed = (ts_end.tv_sec + (double)ts_end.tv_usec / 1000000) - (ts_start.tv_sec + (double)ts_start.tv_usec / 1000000);
		// printf("\nOverall elapsed time: %f\n", elapsed);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    elapsed_sec = ts_end.tv_sec - ts_start.tv_sec;

    if ((ts_end.tv_nsec - ts_start.tv_nsec) < 0) {
       elapsed_nsec = 1000000000 + ts_end.tv_nsec - ts_start.tv_nsec;
    } else {
      elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
    }

    printf("\n   Overall elapsed time: %f <---\n", elapsed_sec + (double)elapsed_nsec/1000000000);


    free(I);
    free(J);
    free(result);
    free(e_vector);
    free(c3_vector);
    free(csc_col);
    free(csc_row);
    free(masked_vals);
    pthread_exit(0);

  	return 0;
}

void *sparseSelfProduct(void *args) {

  uint32_t sum;
  ARRINFO *info = args;
  uint32_t *csc_col = info->csc_col;
  uint32_t *csc_row = info->csc_row;
  uint32_t *masked_vals = info->masked_vals;
  uint32_t start = info->i_start;
  uint32_t end = info->i_end;

  //printf("in thread: %d with index start: %d, end:  %d\n",id, start, end );
  for(int i = start; i < end; i++) {

    uint32_t starti = csc_col[i];
    uint32_t endi = csc_col[i + 1];

    for(int j = 0; j < endi - starti; j++) {

      uint32_t startj = csc_col[csc_row[ csc_col[i] + j ]];
      uint32_t endj = csc_col[ csc_row[ csc_col[i] + j ] + 1];

      uint32_t size_row = endi - starti;
      uint32_t size_col = endj - startj;

      uint32_t *r = malloc(size_row * sizeof(uint32_t));
      uint32_t *c = malloc(size_col * sizeof(uint32_t));

      for(int pos = 0; pos < size_row; pos++) {
        r[pos] = csc_row[starti + pos];
      }

      for(int pos = 0; pos < size_col; pos++) {
        c[pos] = csc_row[startj + pos];
      }

      int ptr_row = 0;
      int ptr_col = 0;
      sum = 0;

      while( ptr_row != size_row && ptr_col != size_col ) {
        if (r[ptr_row] > c[ptr_col]) {
          ptr_col++;
        } else if (r[ptr_row] < c[ptr_col]) {
          ptr_row++;
        } else {
          ptr_row++;
          ptr_col++;
          sum++;
        }
      }
      if(sum) {
        masked_vals[starti + j] = sum;
      }
      free(r);free(c);
    }

  }
  pthread_exit(NULL);
}

void coo2csc(
  uint32_t       * const row,
  uint32_t       * const col,
  uint32_t const * const row_coo,
  uint32_t const * const col_coo,
  uint32_t const         nnz,
  uint32_t const         n,
  uint32_t const         isOneBased
) {

  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  for (uint32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;

  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }

  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }

}
