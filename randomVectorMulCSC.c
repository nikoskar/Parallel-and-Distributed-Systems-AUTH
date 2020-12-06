#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mmio.h"

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

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    uint32_t *I, *J;
    double *val;
    uint32_t total = 0;

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
    }
    }


    if (f !=stdin) fclose(f);

    for(uint32_t i = 0; i < nz; i++) {
        I[nz + i] = J[i];
        J[nz + i] = I[i];
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

    /*
     * Create a random vector of dimension Nx1 (N = number of rows/cols of A)
     * Its cells contain values ranging from 1 to 5 (randomly filled)
     *
    */
    srand(time(NULL));
    uint32_t *vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *result = malloc(N * sizeof(uint32_t));

    for(uint32_t i = 0; i < N; i++) {
      vector[i] = 1 + rand() % 5; // guarantees no zeros
      result[i] = 0;
    }

    for(uint32_t i = 0; i < N; i++)
    {
      uint32_t start = csc_col[i];
      uint32_t end = csc_col[i + 1];
      uint32_t sum = 0;

      //printf("i = %d, start = %d, end = %d\n", i, start, end);

      for(uint32_t j = 0; j < end - start; j++) {
        sum = sum + vector[csc_row[j]];
        result[i] += vector[start + j];
      }
      //printf(" the total sum is %d\n: ", sum);

    }
	return 0; 
  }
