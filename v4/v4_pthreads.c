#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include "mmio.h"
#include <pthread.h>

#define NUM_THREADS 12

typedef struct {
  uint32_t i_start;  // starting index of loop
  uint32_t i_end;    // final index of loop
  uint32_t *csc_col; // CSC column array
  uint32_t *csc_row; // CSC row array
  uint32_t *masked_vals; //masked matrix values
} ARRINFO;

void *mulA(void *args) {

  uint32_t sum;
  ARRINFO *info = args;
  uint32_t *csc_col = info->csc_col;
  uint32_t *csc_row = info->csc_row;
  uint32_t *masked_vals = info->masked_vals;

  for(int i = info->i_start; i < info->i_end; i++) {

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
  pthread_exit(0);

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

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    uint32_t *I, *J;
    double *val;

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

    int sum;
    int count = 0;
    int flag = 1;
		struct timeval ts_start, ts_end;

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

    /* define the struct that will be passed in the threads. It hold all of the
     * info related to the array
    */
    ARRINFO arr_info;
    pthread_t threads[NUM_THREADS];


    /* Create threads to perform the triangle search  */
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    uint32_t group = N / NUM_THREADS;
    int rc;


    gettimeofday(&ts_start,NULL);

    for(uint32_t i = 0; i < NUM_THREADS; i++) {
      ARRINFO arr_info;
      arr_info.i_start = i * group;
      arr_info.i_end = i * group + group;
      arr_info.csc_col = csc_col;
      arr_info.csc_row = csc_row;
      arr_info.masked_vals = masked_vals;
      rc = pthread_create(&threads[i], &attr, mulA, &arr_info);
      if(rc) {
        printf("\nERROR: Something went wrong while creating a thread\n");
        printf("\nThe error code returned by pthread_create is %d", rc);
      }
    }

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



    gettimeofday(&ts_end,NULL);
    double elapsed = (ts_end.tv_sec + (double)ts_end.tv_usec / 1000000) - (ts_start.tv_sec + (double)ts_start.tv_usec / 1000000);

		printf("\nOverall elapsed time: %f\n", elapsed);

    free(I);
    free(J);
    free(result);
    free(e_vector);
    free(c3_vector);
    free(csc_col);
    free(csc_row);
    free(masked_vals);

  	return 0;
}
