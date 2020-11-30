#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[]) {


	int arr[4][4] = {
		{0, 0, 0, 0},
		{5, 8, 0, 0},
		{0, 0, 3, 0},
		{0, 6, 0, 0}
	};
	
	int v[4] = {5, 8, 3, 6};
	int cols[4] = {0, 1, 2, 1};
	int rows[4] = {1, 1, 2, 3};

	int lenrindex = 5; // = 4 + 1 // m + 1 h orisma apo
	int *temp = (int*)malloc(lenrindex*sizeof(int));

	// need to find the number of non zeros in every line	
	int *nnz = (int*)malloc(4); // number of rows
	
	for(int i = 0; i < 4; i++) {
		nnz[i] = 0;
	}

	int flag = rows[0], count;
	count = 1;
	for(int i = 1; i < 4; i++) {
		if(flag == rows[i]) {
			count++;
		} else { 
			nnz[flag] = count;
			flag = rows[i];
			count = 1;
		}
	}

	nnz[3] = count;
	
	int A[5];
	A[0] = 0;
	for(size_t i = 1; i < lenrindex; i++) {
		A[i] = A[i-1] + nnz[i-1];
	}
	
	printf("\n\n");
	for(int i = 0; i < 4; i++)
	{		
		printf("hello: %d\n",nnz[i]);
	}
}

