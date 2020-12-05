#include <time.h>
#include <stdio.h>
#include <stdlib.h>

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
