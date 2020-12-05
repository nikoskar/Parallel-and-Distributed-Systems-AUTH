#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include "mmio.h"

/* The following program reads an undirected, unweighted random graph (Matrix Market file format), 
 * represented by its adjacency matrix (sparse) and counts the triangles incident with each node 
 * using the formula ((A .* (A*A)) * e) / 2. 
*/

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
    
    FILE *f;
    double *val;
    int M, N, nz;
    int ret_code;
    uint32_t *I, *J;
    MM_typecode matcode;
    int sum, count = 0, flag = 1;
    long elapsed_sec, elapsed_nsec;
		struct timespec ts_start, ts_end;
		//struct timeval ts_start, ts_end;

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

    /* Reserve memory for the CSC arrays */
    uint32_t * csc_row = (uint32_t *)malloc(2 * nz     * sizeof(uint32_t));
    uint32_t * csc_col = (uint32_t *)malloc((N + 1) * sizeof(uint32_t));

    /* Figure out which part of the array you fetch(upper, lower) */
    if (I[0] > J[0]) coo2csc(csc_row, csc_col, J, I, 2 * nz, N, 0);
    else coo2csc(csc_row, csc_col, I, J, 2 * nz, N, 0);

    /* Initialize the masked sparse matrix product array (result of the product) */
    uint32_t *masked_vals = malloc(2 * nz * sizeof(uint32_t));
    for(uint32_t i = 0; i < 2 * nz; i++) {
      masked_vals[i] = 0;
    }

    /* Reserve memory and initialize the 3 following arrays: 
     * e_vector --> filled with ones
     * c3_vector --> triangles vector
     * result --> the result of the 3 multiplications
     */
    uint32_t *e_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *c3_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *result = malloc(N * sizeof(uint32_t));

    for(uint32_t i = 0; i < N; i++) {
      e_vector[i] = 1 ;
      result[i] = 0;
      c3_vector[i] = 0;
    }

    printf("\n*** Triangle counting has now started ***\n");

    //gettimeofday(&ts_start,NULL);
    clock_gettime(CLOCK_MONOTONIC, &ts_start);


    for(int i = 0; i < N; i++) {

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

    /* multiply with the ones vector using only the amount 
     * of non zero elements
     */
    for(uint32_t i = 0; i < N; i++)
    {
      uint32_t start = csc_col[i];
      uint32_t end = csc_col[i + 1];

      for(uint32_t j = 0; j < end - start; j++) {
        // the vector can safely be ignored (all elements = 1)
        result[i] += masked_vals[start + j];
      }
    }

    /* multiply by 2 and get the final result */
    uint32_t total = 0;
    for(uint32_t i = 0; i < N; i++) {
      c3_vector[i] = result[i]/2;
      total += c3_vector[i];
    }

    //gettimeofday(&ts_end,NULL);
    //double elapsed = (ts_end.tv_sec + (double)ts_end.tv_usec / 1000000) - (ts_start.tv_sec + (double)ts_start.tv_usec / 1000000);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    elapsed_sec = ts_end.tv_sec - ts_start.tv_sec;

    if ((ts_end.tv_nsec - ts_start.tv_nsec) < 0) {
    	 elapsed_nsec = 1000000000 + ts_end.tv_nsec - ts_start.tv_nsec;
  	} else {
  	  elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
  	}
        
    printf("\n   Total triangles: %d <---\n: ", total/3);
    printf("\n   Overall elapsed time: %f <---\n", elapsed_sec + (double)elapsed_nsec/1000000000);

		//printf("\nOverall elapsed time: %f\n", elapsed);

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
