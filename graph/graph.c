#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include "mmio.h"
#define NUM_NODES 1000

/*
 * we defined NUM_NODES in order to test the code
 * The code can me easily modified to accept NUM_NODES as user input
 * This is why we use malloc
*/

int** buildGraph() {

	srand(time(NULL));
	int **adj = malloc(sizeof(int*)*NUM_NODES);

	for(size_t i = 0; i < NUM_NODES; i++) {
		adj[i] = malloc(NUM_NODES*sizeof(int)); //*(adj+i)
	}

	for(size_t i = 0; i < NUM_NODES; i++) {
		for(size_t j = 0; j < NUM_NODES; j++) {
			if (i == j) adj[i][j] = 0;
			else adj[i][j] = adj[j][i] = rand() % 2;
		}
	}

	return adj;
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


// uint32_t checkifone(int x, int y, uint32_t *row, uint32_t *col) {
// 	uint32_t start = col[y];
// 	uint32_t end = col[y+1];
// 	uint32_t size = end - start;
// 	uint32_t *sliced = (uint32_t *)malloc(size * sizeof(uint32_t));
// 	uint32_t count = 0;
// 	for (int i = start; i < end; i++) {
// 		sliced[count] = row[i];
// 		count++;
// 	}
//
// 	for(int i = 0; i < size; i++) {
// 		if (x == sliced[i]) {
// 			return 1;
// 		}
// 	}
//
// 	return 0;
// }


uint32_t checkifone(int x, int y, uint32_t *row, uint32_t *col) {
	uint32_t nz = col[j+1] + col[j];
	uint32_t counter = 0;
	do {
		if (row[col[y] + counter] == x) return 1;
		if (row[col[y] + counter] > x) return 0;
		counter++;
	} while (nz > x)
	return 0;
}


int main(int argc, char *argv[]) {

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

	if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
			exit(1);

	I = (uint32_t *) malloc(nz * sizeof(uint32_t));
	J = (uint32_t *) malloc(nz * sizeof(uint32_t));

	for (i=0; i<nz; i++)
	{
			fscanf(f, "%d %d \n", &I[i], &J[i]);
			I[i]--;  /* adjust from 1-based to 0-based */
			J[i]--;
	}

	mm_write_banner(stdout, matcode);
	mm_write_mtx_crd_size(stdout, M, N, nz);
	// for (i=0; i<nz; i++)
	//     fprintf(stdout, "%d %d \n", I[i]+1, J[i]+1);

	// NOW CONVERT COO TO CSC
	uint32_t isOneBased = 0;
	uint32_t *csc_row = (uint32_t *)malloc(nz     * sizeof(uint32_t));
	uint32_t *csc_col = (uint32_t *)malloc((M + 1) * sizeof(uint32_t));
	coo2csc(csc_row, csc_col, I, J, nz, M, isOneBased);



	int *c3 = (int *)malloc(sizeof(int)*M);
	for (int i = 0; i < M; i++) {
		c3[i] = 0;
	}
	int **arr = buildGraph();
	long elapsed_sec, elapsed_nsec, elapsed_sec1, elapsed_nsec1;
	struct timespec ts_start, ts_end, ts_start1, ts_end1;

	printf("\n########## ADJ MATRIX ##########\n\n");

	// for (int i = 0; i < NUM_NODES; i++) {
	// 	for (int j = 0; j < NUM_NODES; j++) {
	// 		printf(" %d",arr[i][j]);
	// 	}
	// 	printf("\n");
	// }

	printf ("\n################################\n\n");

	// OPTIMIZED VERSION
	clock_gettime(CLOCK_MONOTONIC, &ts_start);

	for (int i = 1; i < M - 2; i++) {
		for (int j = i + 1; j < M - 1; j++) {
			for (int k = j + 1; k < M; k++) {
				uint32_t temp1, temp2, temp3;
				temp1 = checkifone(i, j, csc_row, csc_col);
				temp2 = checkifone(j, k, csc_row, csc_col);
				temp3 = checkifone(k, i, csc_row, csc_col);
				//arr[i][j] == arr[j][k] == arr[k][i] == 1
				if (temp1 & temp2 & temp3) {
					c3[i]++;
					c3[j]++;
					c3[k]++;
				}
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &ts_end);

	// SLOW VERSION
	// clock_gettime(CLOCK_MONOTONIC, &ts_start1);
	// for (int i = 0; i < M; i++) {
	// 	for (int j = 0; j < M; j++) {
	// 		for (int k = 0; k < M; k++) {
	// 			//if ( arr[i][j] == arr[j][k] == arr[k][i] == 1) {
	// 			if (arr[i][j] == arr[j][k] && arr[j][k] == arr[k][i] && arr[k][i] == 1) {
	// 				c3[i]++;
	// 				c3[j]++;
	// 				c3[k]++;
	// 			}
	// 		}
	// 	}
	// }
  // clock_gettime(CLOCK_MONOTONIC, &ts_end1);


	elapsed_sec = ts_end.tv_sec - ts_start.tv_sec;
	//elapsed_sec1 = ts_end1.tv_sec - ts_start1.tv_sec;
	//elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;

	if ((ts_end.tv_nsec - ts_start.tv_nsec) < 0) {
  	elapsed_nsec = 1000000000 + ts_end.tv_nsec - ts_start.tv_nsec;
  } else {
    elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;
  }

	// if ((ts_end1.tv_nsec - ts_start1.tv_nsec) < 0) {
	// 	elapsed_nsec1 = 1000000000 + ts_end1.tv_nsec - ts_start1.tv_nsec;
	// } else {
	// 	elapsed_nsec1 = ts_end1.tv_nsec - ts_start1.tv_nsec;
	// }

	// for (int i = 0; i < NUM_NODES; i++) {
	// 	printf(" Node %d has %d incident triangles\n", i, c3[i]);
	// }

	printf ("\n################################\n\n");

	// printf("SLOW VERSION\n");
	// printf("Elapsed time in seconds: %ld\n", elapsed_sec1);
	// printf("Elapsed time in nanoseconds: %ld\n", elapsed_nsec1);
	// printf("Overall elapsed time: %f\n", elapsed_sec1 + (double)elapsed_nsec1/1000000000);

	printf("\n\nOPTIMIZED VERSION\n");
	printf("Elapsed time in seconds: %ld\n", elapsed_sec);
	printf("Elapsed time in nanoseconds: %ld\n", elapsed_nsec);
	printf("Overall elapsed time: %f\n", elapsed_sec + (double)elapsed_nsec/1000000000);

	printf ("\n################################\n\n");

	free(arr);
}
