#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_ROWS 5

int main() {
	
	int nz;
	int count = 0;
	srand(time(NULL));
	int **arr = malloc(sizeof(int*)*NUM_ROWS);

	for(size_t i = 0; i < NUM_ROWS; i++)	
		arr[i] = malloc(sizeof(int)*NUM_ROWS);
	
	// generate the following arr just to test the code works
	for (size_t i = 0; i < NUM_ROWS; i++)
		for (size_t j = 0; j < NUM_ROWS; j++)
			if ( i == j) arr[i][j] = i+j+1;
			else if ( (i == 0 || i == 4) && (j == 0 || j == 3)) arr[i][j] = rand() % 10;
			else arr[i][j] = 0;
	
	// allocate space for the 3 arrays (row, col, value) - worst case scenario
	int *row = (int*)malloc(sizeof(int)*NUM_ROWS*NUM_ROWS);
	int *col = (int*)malloc(sizeof(int)*NUM_ROWS*NUM_ROWS);
	int *values = (int*)malloc(sizeof(int)*NUM_ROWS*NUM_ROWS);

	for (size_t i = 0; i < NUM_ROWS; i++) {
		for (size_t j = 0; j < NUM_ROWS; j++) {
			printf("%d ",arr[i][j]);
		
			if ( arr[i][j] != 0 ) {
				row[count] = i;
				col[count] = j;
				values[count] = arr[i][j];
				count++;
			}
		}	
		printf("\n");
	}
	printf("%d\n",count);


	for (size_t i = 0; i < count; i++) {
			printf("row: %d, column: %d, value: %d\n", row[i], col[i], values[i]);
	}
	
	
	// COO to dense makes no sense. So lets turn arr into a CSC













	
	free(arr)
	free(row);	
	free(col);
	free(values);
	
			
}
