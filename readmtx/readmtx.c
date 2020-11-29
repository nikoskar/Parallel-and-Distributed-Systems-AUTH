
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "mmio.h"

#define foreach(item, array) \
    for(int keep = 1, \
            count = 0,\
            size = sizeof (array) / sizeof *(array); \
        keep && count != size; \
        keep = !keep, count++) \
      for(item = (array) + count; keep; keep = !keep)


int *printsliced (uint32_t start,
                  uint32_t end,
                  uint32_t i,
                  uint32_t *csc_col,
                  uint32_t *csc_row,
                  uint32_t N) {


  int size = end - start;
  int *sliced = malloc(size*sizeof(int));
  int *none = malloc(1*sizeof(int));
  none[0] = 0;
  int flag = 0;



  for (int x = 0; x < size; x++) {
    sliced[x] = 0;
  }


  int id = 0;
  for(uint32_t x = start; x < end; x++) {
    if(csc_row[x] > i && csc_row[x] < N-1)
    {
      sliced[id] = csc_row[x];
      id++;
      flag = 1;
    }
  }

  int *final_sliced = malloc(size*sizeof(int));
  int final_size = id;
  for (int x = 0; x < size; x++) {
    final_sliced[x] = 0;
  }

  for(int x = 0; x < id; x++) {
    final_sliced[x] = sliced[x];
  }
  // printf("\nrow : %d, SLICED FINAL IS: ", i);
  //
  // for(int x = 0; x < final_size; x++) {
  //   printf("%d, ", final_sliced[x]);
  // }
  if(flag) {
    return final_sliced;
  } else {
    return none;
  }

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
    uint32_t M, N, nz;
    uint32_t i, *I, *J;
    double *val;

    if (argc < 2)
	  {
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	  } else {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /* find out size of sparse matrix .... */
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);

    /* reseve memory for matrices */
    I = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(2*nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));

    if (!mm_is_pattern(matcode))
    {
      for (i=0; i<2*nz; i++) {
          fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
          I[nz+i] = J[i];
          J[nz+i] = I[i];
          I[i]--;  /* adjust from 1-based to 0-based */
          J[i]--;
          J[nz+i]--;
          I[nz+i]--;
      }
    }
    else{
      for (i=0; i<nz; i++)
      {
          fscanf(f, "%d %d\n", &I[i], &J[i]);
          I[nz+i] = J[i];
          J[nz+i] = I[i];
          val[i]=1;
          I[i]--;  /* adjust from 1-based to 0-based */
          J[i]--;
          J[nz+i]--;
          I[nz+i]--;
      }
    }

    if (f !=stdin) fclose(f);

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);


    uint32_t isOneBased = 0;
    uint32_t *csc_row = (uint32_t *)malloc(2*nz*sizeof(uint32_t));
    uint32_t *csc_col = (uint32_t *)malloc((N + 1)*sizeof(uint32_t));

    coo2csc(csc_row, csc_col, I, J, nz*2, N, isOneBased);

    // printf("print csc_row(must me %d elements):\n", 2*nz);
    // for (int i = 0; i < 2*nz; i++) {
    //   printf("%d :  %d\n", i, csc_row[i]);
    // }
    //
    // printf("print csc_col(must me %d elements):\n", N+1);
    // for (int i = 0; i < N + 1; i++) {
    //   printf("%d :  %d\n", i, csc_col[i]);
    // }
    printf("\n\n");

    /* TRIANGLE COUNTER AREA  */
    long elapsed_sec, elapsed_nsec;
  	struct timespec ts_start, ts_end;

    uint32_t *c3 = (uint32_t*)malloc(N*sizeof(uint32_t));
    for(int i = 0; i < N; i++) {
      c3[i] = 0;
    }

    printf("\n\n");
    printf("size of array is %d x %d", N,N);
    printf("started triangle counting\n");
    clock_gettime(CLOCK_MONOTONIC, &ts_start);


    for (int i = 0; i < N-2; i++) {

      uint32_t start = csc_col[i];
      uint32_t end = csc_col[i+1];


      int *sliced = printsliced(start, end, i, csc_col, csc_row, N);
      int size = sizeof(sliced)/sizeof(sliced[0]);
      printf("\n\n");

      printf("currently in loop i: %d", i);

      // check for the zero case
      if (sliced[0] > i) {
        foreach(int *j, sliced) {

          uint32_t start2 = csc_col[*j];
          uint32_t end2 = csc_col[*j+1];

          int *sliced2 = printsliced(start2, end2, *j, csc_col, csc_row, N+1);
          int size2 = sizeof(sliced2)/sizeof(sliced2[0]);
          printf("\n\n");

          printf("******currently in loop j: %d", *j);
          // foreach(int *j, sliced2) {
          //   printf(" %d, ", *j);
          // }
          printf("\n\n");

          foreach (int *k, sliced2) {
            printf("******--------------currently in loop k: %d", *k);
            printf("\n\n");

            c3[*k]++;

          }


          free(sliced2);

        }
      }




      free(sliced);

    }

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    printf("finished triangle counting\n");


    printf("\n\n");
    elapsed_sec = ts_end.tv_sec - ts_start.tv_sec;

    if ((ts_end.tv_nsec - ts_start.tv_nsec) < 0) {
      elapsed_nsec = 1000000000 + ts_end.tv_nsec - ts_start.tv_nsec;
    } else {
      elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
    }
    printf("Overall elapsed time: %f\n", elapsed_sec + (double)elapsed_nsec/1000000000);

    for (int i = 0; i < N; i++) {
  		printf(" Node %d has %d incident triangles\n", i, c3[i]);
  	}

	return 0;
}







/*

























































*/
