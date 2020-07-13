#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define A_MATRIX_HEIGHT 1000
#define A_MATRIX_WIDTH 1000
#define B_MATRIX_HEIGHT 1000
#define B_MATRIX_WIDTH A_MATRIX_WIDTH
#define C_MATRIX_HEIGHT A_MATRIX_HEIGHT
#define C_MATRIX_WIDTH B_MATRIX_HEIGHT

void print_matrix(int *matrix, size_t h, size_t w)
{
	for(size_t i = 0; i < h; ++i) {
		for(size_t j = 0; j < w; ++j) {
			printf("%d ", matrix[i * w + j]);
		}
		printf("\n");
	}
}

int main(int argc, char **argv)
{
	int num_of_procs, *a, *b, *c, i, j, k;
	double start_wtime = omp_get_wtime(), end_wtime = 0.0;
	a = (int*)malloc(sizeof(int)*A_MATRIX_HEIGHT*A_MATRIX_WIDTH);
	b = (int*)malloc(sizeof(int)*B_MATRIX_HEIGHT*B_MATRIX_WIDTH);
	c = (int*)malloc(sizeof(int)*C_MATRIX_HEIGHT*C_MATRIX_WIDTH);
	for(i = 0; i < A_MATRIX_HEIGHT; ++i)
		for(j = 0; j < A_MATRIX_WIDTH; ++j)
			a[i * A_MATRIX_WIDTH + j] = i * A_MATRIX_WIDTH + j;
	for(i = 0; i < B_MATRIX_HEIGHT; ++i)
		for(j = 0; j < B_MATRIX_WIDTH; ++j)
			b[i * B_MATRIX_WIDTH + j] = i * B_MATRIX_WIDTH + j;
#pragma omp parallel for private(i, j, k)
	for(i = 0; i < A_MATRIX_WIDTH; ++i) {
		for(j = 0; j < B_MATRIX_WIDTH; ++j) {
			int temp = 0;
			num_of_procs = omp_get_num_threads();
			for(k = 0; k < A_MATRIX_WIDTH; ++k)
				temp += a[i*A_MATRIX_WIDTH + k] * b[j*B_MATRIX_WIDTH + k];
			c[i*C_MATRIX_WIDTH + j] = temp;
		}
	}
	free(a);
	free(b);
	free(c);
	end_wtime = omp_get_wtime();
	printf("%d %d %d %d %d %f\n", num_of_procs,
		A_MATRIX_WIDTH, A_MATRIX_HEIGHT, B_MATRIX_WIDTH, B_MATRIX_HEIGHT, end_wtime - start_wtime);
	return 0;
}
