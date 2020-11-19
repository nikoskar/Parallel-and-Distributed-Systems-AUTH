#include <stdio.h>
#include <stdlib.h>

int coo2csr(double **sparse, int num_rows ,size_t num_nnz, int *rows, int *cols, double* vals) {

	int count, flag;

	int *row_index = (int*)malloc((num_rows+1)*sizeof(int));
	int *nnz = (int*)malloc(num_rows*sizeof(int));

	//can go into thread 1
	for (size_t i = 0; i < num_rows; i++) {
		nnz[i] = 0;
	}

	//can go into thread 2
	flag = rows[0]; count = 1;
	for(size_t i = 1; i < num_nnz; i++) {
		if (flag == rows[i]) {
			count++;
		} else {
			nnz[flag] = count;
			flag = rows[i];
			count = 1;
		}
	}
	nnz[num_nnz - 1] = count;

	row_index[0] = 0;
	for (size_t i = 1; i < num_rows + 1; i++) {
		row_index[i] = row_index[i-1] + nnz[i-1];
	}
}

