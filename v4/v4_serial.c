#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mmio.h"

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
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

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, *I, *J;
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

    I = (int *) malloc(2 * nz * sizeof(int));
    J = (int *) malloc(2 * nz * sizeof(int));
    val = (double *) malloc(2 * nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    /* Replace missing val column with 1s and change the fscanf to match patter matrices*/

    if (!mm_is_pattern(matcode))
    {
	    for (i=0; i<nz; i++)
	    {
	        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
	        I[i]--;  /* adjust from 1-based to 0-based */
	        J[i]--;
	    }
    }
    else
    {
	    for (i=0; i<nz; i++)
	    {
	        fscanf(f, "%d %d\n", &I[i], &J[i]);
					I[nz + i] = J[i];
					J[nz + i] = I[i];
					val[nz + i] = val [i] = 1;
	        I[i]--;  /* adjust from 1-based to 0-based */
	        J[i]--;
					I[nz + i]--;
					J[nz + i]--;
	    }
    }



    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    // for (i=0; i< 2 * nz; i++)
    //     fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);


		int tempi, tempj;
		for(int i = 0; i < 2 * nz; i++)  {
			for (int j = i + 1; j < 2 * nz; j++){
		 		if (J[i] > J[j]) {
					tempj = J[i];
					tempi = I[i];

					J[i] = J[j];
					I[i] = I[j];

					J[j] = tempj;
					I[j] = tempi;
				}
			}
		}


    uint32_t * csc_row = (uint32_t *)malloc(2 * nz     * sizeof(uint32_t));
    uint32_t * csc_col = (uint32_t *)malloc((N + 1) * sizeof(uint32_t));

    coo2csc(csc_row, csc_col, I, J, 2 * nz, N, 0);


		for(int i = 0; i < N; i++) {
			int start = csc_col[i];
			int end = csc_col[i+1];
			int count = 0;
			int *temp = malloc( (end - start) * sizeof(int));

			for (int j = start; j < end; j++) {
				temp[count] = csc_row[j];
				count++;
			}
			qsort(temp, end-start, sizeof(int), compare);
			count = 0;
			for (int j = start; j < end; j++) {
				csc_row[j] = temp[count];
				count++;
			}
		}


    // printf("csc_row: \n\n");
    // for(int i = 0; i < 2* nz; i++) {
    //   printf("%d, ", csc_row[i]);
    // }
    //
    // printf("csc_col: \n\n");
    // for(int i = 0; i < N + 1; i++) {
    //   printf("%d, ", csc_col[i]);
    // }
		// printf("\n");


    /*
     * Create a random vector of dimension Nx1 (N = number of rows/cols of A)
     * Its cells contain values ranging from 1 to 5 (randomly filled)
     *
    */
    // srand(time(NULL));
    // uint32_t *vector = (malloc)(N * sizeof(uint32_t));
    // uint32_t *result = malloc(N * sizeof(uint32_t));
    //
    // for(uint32_t i = 0; i < N; i++) {
    //   vector[i] = 1 + rand() % 5; // guarantees no zeros
    //   result[i] = 0;
    // }
    //
    // for(uint32_t i = 0; i < N; i++)
    // {
    //   uint32_t start = csc_col[i];
    //   uint32_t end = csc_col[i + 1];
    //   uint32_t sum = 0;
    //
    //   //printf("i = %d, start = %d, end = %d\n", i, start, end);
    //
    //   for(uint32_t j = start; j < end; j++) {
    //     sum = sum + vector[csc_row[j]];
    //   }
    //   result[i] = sum;
    //   //printf(" the total sum is %d\n: ", sum);
    //
    // }
    //
    //
    // printf("\n");





















    uint32_t *masked_vals = malloc(2 * nz * sizeof(uint32_t));

    //optimize later
    for(uint32_t i = 0; i < 2 * nz; i++) {
      masked_vals[i] = 0;
    }

    //int count = 0, count1 = 0;
    int sum, count = 0;
    int flag = 1;

    // FIll the array who is the result of the multiplication
    for(int i = 0; i < N; i++) {

      /* For any row i, I only want to calculate the columns that are not
       * 0 in A, since the Hadamard product will erase them anyways.
       * The first step is filtering out the lines im not gonna use
      */
      uint32_t starti = csc_col[i];
      uint32_t endi = csc_col[i + 1];

      /* Create a pointer to the list of non zero indices for line i */
      uint32_t *pti = & csc_row[starti];
      uint32_t *pti_temp = & csc_row[starti];
      printf("\neimaste sti grammi %d, me prwto mi mideniko stoixeio to %d\n", i, *pti);


      // now loop through only the useful lines
      for(int j = starti; j < endi; j++) {
        pti = pti_temp;

        /* I am now parsing only the correct columens: csc_row[j]
         * For each and every one of these, i need to get the indices of
         * their non zero elements. In order to get the index of the current
         * column, I need to use csc_row[j]
         * --The following is specific for our test mtx--
         * For example for i = o ==> j = 0:4 . This index corresponds
         * to the columns csc_row[0] = 1, csc_row[1] = 2 etc
        */



        uint32_t startj = csc_col[ csc_row[j] ];
        uint32_t endj = csc_col[ csc_row[j] + 1];

        /* Create a pointer to the list of non zero indices for line i */
        uint32_t *ptj = & csc_row[startj];

        printf("\n-- stili : %d me prwto mi mideniko %d \n", csc_row[j], *ptj);


        int counti = 0, countj = 0;
        sum = 0;
        while( (counti < endi - starti) && (countj < endj - startj) ) {

          if (*ptj < *pti) {
            ptj++;
            countj++;
          } else if (*pti < *ptj) {
            pti++;
            counti++;
          } else {
            counti++;
            countj++;
            ptj++;
            pti++;
            sum = sum + 1;
            printf("xxxxxxxxxxxxxxxxxxxx\n");
          }


          // if (*p > p2) {
          //   p2 = *(ptj++);
          //   countj++;
          // } else if (p1 < p2) {
          //   p1 = *(pti++);
          //   counti++;
          // } else {
          //   countj++;
          //   counti++;
          //   p2 = *(ptj++);
          //   p1 = *(pti++);
          //   sum = sum + 1;
          //   printf("oooooooooooooooooooooooooo\n");
          //   printf("%d\n",sum);
          // }
        }
        int y = sum;
        printf("y is %d\n", y);

        masked_vals[count] = y;
        count++;


      }

    }


    /* vector array */
    uint32_t *e_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *c3_vector = (malloc)(N * sizeof(uint32_t));
    uint32_t *result = malloc(N * sizeof(uint32_t));

    for(uint32_t i = 0; i < N; i++) {
      e_vector[i] = 1 ; // guarantees no zeros
      result[i] = 0;
    }

    for(uint32_t i = 0; i < N; i++)
    {
      uint32_t start = csc_col[i];
      uint32_t end = csc_col[i + 1];
      uint32_t s = 0;

      //printf("i = %d, start = %d, end = %d\n", i, start, end);

      for(uint32_t j = start; j < end; j++) {
        //s = s + vector[ csc_row[j] ] ;
        // e_vector is a (ones) X N matrix so we dont need to add in the multiplication
        s = s + masked_vals[csc_row[j]];
      }
      result[i] = s;
      c3_vector[i] = s/2;
      printf(" triangles per node:%d\n: ", c3_vector[i]);

    }



    // printf("\n\n\n\ncsc_row: \n");
    // for(int i = 0; i < 2* nz; i++) {
    //   printf("%d, ", csc_row[i]);
    // }
    // printf("\n");
    //
    // printf("csc_col: \n");
    // for(int i = 0; i < N + 1; i++) {
    //   printf("%d, ", csc_col[i]);
    // }
    // printf("\n");
    //
    // printf("masked_vals: \n");
    // for(int i = 0; i < 2 * nz; i++) {
    //   printf("%d, ", masked_vals[i]);
    // }
    // printf("\n");


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


/*
uint32_t *c3 = malloc(N*sizeof(uint32_t));
for (uint32_t i = 0; i < N; i++) {
  c3[i] = 0;
}
for(uint32_t j = 0; j < N - 1; j++) {
  uint32_t start = csc_col[j];
  uint32_t end = csc_col[j + 1];
  uint32_t idx = 0;
  uint32_t size = end - start;
  uint32_t *log = malloc(size*sizeof(uint32_t));
  //printf("\n\n");
  printf("j: %d has these zeros: ", j);
  for(uint32_t pos = start; pos < end; pos++) {
    log[idx] = csc_row[pos];
    printf("%d, ", log[idx]);
    idx++;
  }
  //printf("\n\n");
  for(uint32_t i = start; i < end; i++) {
    uint32_t start1 = csc_col[csc_row[i]];
    uint32_t end1 = csc_col[csc_row[i] + 1];
    //printf("in i = %d with start = %d and end = %d", i, start1, end1);
    //printf("\n ---i = %d", csc_row[i]);
    for(uint32_t k = start1; k < end1; k++) {
      //printf("\n ------k = %d (j = %d, i = %d)", csc_row[k], j, csc_row[i]);
      if(csc_row[k]) {
        for(uint32_t x = 0; x < size; x++)
        {
          if(csc_row[k] == log[x]) {
            printf("\n ------k = %d (j = %d, i = %d)", csc_row[k], j, csc_row[i]);
            c3[csc_row[k]]++;
            c3[j]++;
            c3[csc_row[i]]++;
            printf("\ntriangle found\n");
          }
        }
      }
      // for(uint32_t x = 0; x < size; x++)
      // {
      // 	if(csc_row[k] == log[x]) {
      // 		c3[csc_row[k]]++;
      // 		c3[csc_row[j]]++;
      // 		c3[csc_row[i]]++;
      // 		//printf("triangle found\n");
      // 	}
      // }
    }
    printf("\n");
  }
  printf("\n");
  free(log);
}
for (uint32_t i = 0; i < N; i++) {
  printf("node %d has %d triangles\n", i, c3[i]);
}
free(c3);
*/
