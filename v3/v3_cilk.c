#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mmio.h"
#include <pthread.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

/* Parallel Implementation of V3 using OpenCilk */

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
    uint32_t M, N, nz;
    uint32_t i, *I, *J;
    double *val;
		int ret_code;
		MM_typecode matcode;
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

    I = (uint32_t *) malloc(nz * sizeof(uint32_t));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));


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
          val[i]=1;
          I[i]--;  /* adjust from 1-based to 0-based */
          J[i]--;
      }
    }



    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    // for (i=0; i<nz; i++)
    //     fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);

/*****************************************************************************/

    uint32_t * csc_row = (uint32_t *)malloc(nz     * sizeof(uint32_t));
    uint32_t * csc_col = (uint32_t *)malloc((N + 1) * sizeof(uint32_t));

    coo2csc(csc_row, csc_col, I, J, nz, N, 0);

    // printf("csc_row: \n\n");
    // for(int i = 0; i < nz; i++) {
    //   printf("%d, ", csc_row[i]);
    // }
		//
    // printf("csc_col: \n\n");
    // for(int i = 0; i < N + 1; i++) {
    //   printf("%d, ", csc_col[i]);
    // }

		uint32_t *c3 = malloc(N*sizeof(uint32_t));
		for (uint32_t i = 0; i < N; i++) {
			c3[i] = 0;
		}

    __cilkrts_set_param("nworkers","8");
		int numWorkers = __cilkrts_get_nworkers();
    printf("number of workers to %d.\n",numWorkers);

		pthread_mutex_t mutex1;
		pthread_mutex_init(&mutex1, NULL);

    printf("\n***Now we start counting triangles***\n");

	
    /*
    	if you are running this locally you should be able to use clock_gettime without issues. In case you face problems
	with that, we have added gettimeofday. Just uncomment ---> //struct timeval ts_start, ts_end; <--- and the 
	gettimeofday commands. Finally comment out any clock_gettime implementation
    
    */
    //gettimeofday(&ts_start,NULL);
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    cilk_for(uint32_t j = 0; j < N - 1; j++) {

      uint32_t start = csc_col[j];
      uint32_t end = csc_col[j + 1];
			uint32_t idx = 0;
			uint32_t size = end - start;

			uint32_t *log = malloc(size*sizeof(uint32_t));

			//printf("\n\n");

			for(uint32_t pos = start; pos < end; pos++) {
				log[idx] = csc_row[pos];
				//printf("%d, ", log[idx]);
				idx++;
			}
			//printf("\n\n");

      cilk_for(uint32_t i = start; i < end; i++) {

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
								pthread_mutex_lock(&mutex1);
								c3[csc_row[k]]++;
								c3[j]++;
								c3[csc_row[i]]++;
								pthread_mutex_unlock(&mutex1);
							}
						}
					}

        }
      }
			free(log);
    }

		int sum = 0;
		for (uint32_t i = 0; i < N; i++) {
			//printf("node %d has %d triangles\n", i, c3[i]);
			sum += c3[i];
		}


    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    elapsed_sec = ts_end.tv_sec - ts_start.tv_sec;

    if ((ts_end.tv_nsec - ts_start.tv_nsec) < 0) {
    	 elapsed_nsec = 1000000000 + ts_end.tv_nsec - ts_start.tv_nsec;
  	} else {
  	  elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
  	}

    printf("\n   Overall elapsed time: %f <---\n", elapsed_sec + (double)elapsed_nsec/1000000000);
    printf("\n   The total amount of triangles is: %d <---\n\n", sum/3);

		// gettimeofday(&ts_end,NULL);
		// double elapsed = (ts_end.tv_sec + (double)ts_end.tv_usec / 1000000) - (ts_start.tv_sec + (double)ts_start.tv_usec / 1000000);

    free(I);
    free(J);
		free(c3);
    free(val);
    free(csc_col);
    free(csc_row);

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
