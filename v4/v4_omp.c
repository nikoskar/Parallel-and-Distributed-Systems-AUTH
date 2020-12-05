#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mmio.h"
#include <omp.h>
#include <sys/time.h>

/* Parallel implementation of v4 using OpenMP */

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
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    uint32_t *I, *J;
    double *val;
    uint32_t sum;
    long elapsed_sec, elapsed_nsec;
    struct timespec ts_start, ts_end;

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

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */
    I = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    if (!mm_is_pattern(matcode))
    {
      for (uint32_t i=0; i<nz; i++)
      {
          fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
          I[i]--;  /* adjust from 1-based to 0-based */
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
          I[i]--;  /* adjust from 1-based to 0-based */
          J[i]--;
          I[nz + i] = J[i];
          J[nz + i] = I[i];
      }
    }

    /************************/
    /* now write out matrix */
    /************************/
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);


    // for (uint32_t i=0; i< 2 * nz; i++)
    //     fprintf(stdout, "%d %d \n", I[i]+1, J[i]+1);

    uint32_t * csc_row = (uint32_t *)malloc(2 * nz     * sizeof(uint32_t));
    uint32_t * csc_col = (uint32_t *)malloc((N + 1) * sizeof(uint32_t));

    // Figure out which part of the matrix you fetch from the mtx file
    // upper if I[0] < J[0] and lower if I[0] > J[0] and then invoke coo2csc
    // acoordingly
    if (I[0] > J[0]) coo2csc(csc_row, csc_col, J, I, 2 * nz, N, 0);
    else coo2csc(csc_row, csc_col, I, J, 2 * nz, N, 0);



    uint32_t *masked_vals = malloc(2 * nz * sizeof(uint32_t));
    for(uint32_t i = 0; i < 2 * nz; i++) {
      masked_vals[i] = 0;
    }

    /* ones vector, c3 vector and the final triangle vector */
    uint32_t *e_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *c3_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *result = malloc(N * sizeof(uint32_t));

    for(uint32_t i = 0; i < N; i++) {
      e_vector[i] = 1 ; // guarantees no zeros
      result[i] = 0;
      c3_vector[i] = 0;
    }

    printf("\n*** Triangle counting has now started ***\n");
		/* Start counting triangles */
		clock_gettime(CLOCK_MONOTONIC, &ts_start);

  #pragma omp parallel
  {
    /*
     we have default(shared) but default(none) is usefull if we want to debug
     We will use static scheduling. If we could predetermine and predict it then
     it would propably be a good option to try static schedule too. Also, we need
     to add the sum as a private variable. We must be careful because valiables listed in the private clause
     are not initialized unless we make sure for they do. We dont need to include the
     shared variables (csc_col for example) because they are by default shared among the
     threads.
     */
    //omp_set_nested(1);
    #pragma omp for schedule(dynamic) private(sum)

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
            //#pragma omp critical // we want only one thread at a time to access the variable below
            //pragma opm atomic
            masked_vals[starti + j] = sum;
          }
          free(r);free(c);
        }

      } // omp barriers are implied at the end of each loop so we dont have to hardtype it.

  }



    for(uint32_t i = 0; i < N; i++)
    {
      uint32_t start = csc_col[i];
      uint32_t end = csc_col[i + 1];

      for(uint32_t j = 0; j < end - start; j++) {
        //e_vector's values are all 1 so it can be ignored for the multiplication

        result[i] += masked_vals[start + j];

      }
    }


    uint32_t total = 0;
      //#pragma omp parallel
      //#pragma omp for schedule(static) reduction(+:total)
      for(uint32_t i = 0; i < N; i++) {
        //#pragma omp critical
        c3_vector[i] = result[i]/2;
        total += c3_vector[i];
      }

    printf("\n*** Total triangles: %d ***\n: ", total/3);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);


		elapsed_sec = ts_end.tv_sec - ts_start.tv_sec;

		if ((ts_end.tv_nsec - ts_start.tv_nsec) < 0) {
		  elapsed_nsec = 1000000000 + ts_end.tv_nsec - ts_start.tv_nsec;
		} else {
		  elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
		}

		printf("\nOverall elapsed time: %f\n", elapsed_sec + (double)elapsed_nsec/1000000000);

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

void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

  // ----- cannot assume that input is already 0!
  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (uint32_t l = 0; l < nnz; l++)
    col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }

}
