#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#define NUM_NODES 5

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

int main(int argc, char *argv[]) {

	int c3[NUM_NODES] = {0};
	int **arr = buildGraph();
	long elapsed_sec, elapsed_nsec;
	struct timespec ts_start, ts_end;
	
	printf("\n ########## ADJ MATRIX ##########\n");

	for (int i = 0; i < NUM_NODES; i++) {
		for (int j = 0; j < NUM_NODES; j++) {
			printf(" %d",arr[i][j]);
		}
		printf("\n");
	}
	printf (" ################################\n\n");

	clock_gettime(CLOCK_MONOTONIC, &ts_start);

	// START COUNTING THE TRIANGLES
	for (int i = 0; i < NUM_NODES; i++) {
		for (int j = 0; j < NUM_NODES; j++) {
			if(i >= j) {
				continue;
			}
			for (int k = 0; k < NUM_NODES; k++) {
				if(j >= k) {
					continue;
				}
				//if ( arr[i][j] == arr[j][k] == arr[k][i] == 1) {
				if (arr[i][j] == arr[j][k] && arr[j][k] == arr[k][i] && arr[k][i] == 1) {
					c3[i]++; 
					c3[j]++; 
					c3[k]++;
				}
			}
		}
	}	
	
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	elapsed_sec = ts_end.tv_sec - ts_start.tv_sec;
	elapsed_nsec = ts_end.tv_nsec - ts_start.tv_nsec;

	printf("Elapsed time in seconds: %ld\n", elapsed_sec);
	printf("Elapsed time in nanoseconds: %ld\n", elapsed_nsec);
	printf (" \n################################\n\n");
	
	for (int i = 0; i < NUM_NODES; i++) {
		printf(" Node %d has %d incident triangles\n", i, c3[i]);
	}

	printf (" \n################################\n\n");

	free(arr);
}
