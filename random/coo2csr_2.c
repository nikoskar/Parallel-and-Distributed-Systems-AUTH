#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void coo2csr(
	uint32_t 	   num, // number of rows/columns
	uint32_t 	   num_nnz, // number of non zero elements
	uint32_t 	   *const csr_row, // array to fill
	uint32_t 	   *const csr_col,
	uint32_t const *const rows, // COO rows
	uint32_t const *const cols, // COO columns
	uint32_t const *const vals // COO values

){
	uint32_t count, flag;
	uint32_t *nnz = (uint32_t*)malloc(num*sizeof(uint32_t));

	for (uint32_t i = 0; i < num; i++) {
		nnz[i] = 0;
	}

	flag = rows[0]; count = 1;
	for(uint32_t i = 1; i < num; i++) {
		if (flag == rows[i]) {
			count++;
		} else {
			nnz[flag] = count;
			flag = rows[i];
			count = 1;
		}
	}
	nnz[num - 1] = count;

	csr_row[0] = 0;
	for (uint32_t i = 1; i < num + 1; i++) {
		csr_row[i] = csr_row[i-1] + nnz[i-1];
	}
	free(nmz);
}
