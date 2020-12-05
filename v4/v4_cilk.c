#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include "mmio.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <pthread.h>


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

    int numWorkers = __cilkrts_get_nworkers();
    printf("\n*** Currently using %d workers ***.\n",numWorkers);


    printf("\n*** Triangle counting has now started ***\n");

    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);

    gettimeofday(&ts_start,NULL);


    for (uint32_t i = 0; i < N; i++) {

      uint32_t starti = csc_col[i];
      uint32_t endi = csc_col[i + 1];

      for (uint32_t j = 0; j < endi - starti; j++) {

        uint32_t startj = csc_col[csc_row[ csc_col[i] + j ]];
        uint32_t endj = csc_col[ csc_row[ csc_col[i] + j ] + 1];

        uint32_t size_row = endi - starti;
        uint32_t size_col = endj - startj;

        uint32_t *r = malloc(size_row * sizeof(uint32_t));
        uint32_t *c = malloc(size_col * sizeof(uint32_t));


        for(uint32_t pos = 0; pos < size_row; pos++) {
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


    cilk_for(uint32_t i = 0; i < N; i++)
    {
      uint32_t start = csc_col[i];
      uint32_t end = csc_col[i + 1];

      for(uint32_t j = 0; j < end - start; j++) {
        result[i] += masked_vals[start + j];
      }
    }

    uint32_t total = 0;
    cilk_for(uint32_t i = 0; i < N; i++) {
      c3_vector[i] = result[i]/2;
      /* If we add the mutex here we increase the time by ~ 0.1 sec*/
      //pthread_mutex_lock(&mutex);
      total += c3_vector[i];
      //pthread_mutex_unlock(&mutex);

    }


    gettimeofday(&ts_end,NULL);

    printf("\n*** Total triangles: %d ***\n: ", total/3);



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
