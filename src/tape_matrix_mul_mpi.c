#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

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
	int my_rank, num_of_procs, *a, *b, *c;
	int *a_sendcounts, *a_displs, *b_sendcounts, *b_displs, *c_sendcounts, *c_displs;
	double start_wtime = 0.0, end_wtime = 0.0, *all_times;
	MPI_Init(&argc, &argv);
	start_wtime = MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &num_of_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int a_div = A_MATRIX_HEIGHT / num_of_procs, a_mod = A_MATRIX_HEIGHT % num_of_procs;
	int b_div = B_MATRIX_HEIGHT / num_of_procs, b_mod = B_MATRIX_HEIGHT % num_of_procs;
	if (my_rank == 0) {
		a = (int*)malloc(sizeof(int)*A_MATRIX_HEIGHT*A_MATRIX_WIDTH);
		b = (int*)malloc(sizeof(int)*B_MATRIX_HEIGHT*B_MATRIX_WIDTH);
		c = (int*)malloc(sizeof(int)*C_MATRIX_HEIGHT*C_MATRIX_WIDTH);
		for(size_t i = 0; i < A_MATRIX_HEIGHT; ++i)
			for(size_t j = 0; j < A_MATRIX_WIDTH; ++j)
				a[i * A_MATRIX_WIDTH + j] = i * A_MATRIX_WIDTH + j;
		for(size_t i = 0; i < B_MATRIX_HEIGHT; ++i)
			for(size_t j = 0; j < B_MATRIX_WIDTH; ++j)
				b[i * B_MATRIX_WIDTH + j] = i * B_MATRIX_WIDTH + j;
		a_sendcounts = (int*)malloc(sizeof(int)*num_of_procs);
		b_sendcounts = (int*)malloc(sizeof(int)*num_of_procs);
		c_sendcounts = (int*)malloc(sizeof(int)*num_of_procs);
		a_displs = (int*)malloc(sizeof(int)*num_of_procs);
		b_displs = (int*)malloc(sizeof(int)*num_of_procs);
		c_displs = (int*)malloc(sizeof(int)*num_of_procs);
		all_times = (double*)malloc(sizeof(double)*num_of_procs);
		for(size_t i = 0; i < num_of_procs; ++i) {
			a_sendcounts[i] = (a_div + (i < a_mod ? 1 : 0)) * A_MATRIX_WIDTH;
			b_sendcounts[i] = (b_div + (i < b_mod ? 1 : 0)) * B_MATRIX_WIDTH;
			c_sendcounts[i] = (a_div + (i < a_mod ? 1 : 0)) * C_MATRIX_WIDTH;
		}
		a_displs[0] = 0;
		b_displs[0] = 0;
		c_displs[0] = 0;
		for(size_t i = 1; i < num_of_procs; ++i) {
			a_displs[i] = a_displs[i - 1] + a_sendcounts[i - 1];
			b_displs[i] = b_displs[i - 1] + b_sendcounts[i - 1];
			c_displs[i] = c_displs[i - 1] + c_sendcounts[i - 1];
		}
	}
	int a_block_size = a_div + (a_mod != 0 ? 1 : 0);
	int b_block_size = b_div + (b_mod != 0 ? 1 : 0);
	int *a_recv = (int*)malloc(sizeof(int)*a_block_size*A_MATRIX_WIDTH);
	int *b_recv = (int*)malloc(sizeof(int)*b_block_size*B_MATRIX_WIDTH);
	int *c_recv = (int*)malloc(sizeof(int)*a_block_size*C_MATRIX_WIDTH);
	int a_block_recv = a_div + (my_rank < a_mod ? 1 : 0); 
	int b_block_recv = b_div + (my_rank < b_mod ? 1 : 0);
	MPI_Scatterv(a, a_sendcounts, a_displs, MPI_INT, a_recv, a_block_recv*A_MATRIX_WIDTH, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(b, b_sendcounts, b_displs, MPI_INT, b_recv, b_block_recv*B_MATRIX_WIDTH, MPI_INT, 0, MPI_COMM_WORLD);
	int temp, j_start = (b_block_size - 1) * my_rank + (my_rank < b_mod ? my_rank : b_mod);
	int next_rank = my_rank == 0 ? num_of_procs - 1 : my_rank - 1;
	int prev_rank = my_rank == num_of_procs - 1 ? 0 : my_rank + 1; 
	MPI_Status status;
	for(size_t n = 0; n < num_of_procs; ++n) {
		for(size_t i = 0; i < a_block_recv; ++i) {
			for(size_t j = 0; j < b_block_recv; ++j) {
				temp = 0;
				for(size_t k = 0; k < A_MATRIX_WIDTH; ++k)
					temp += a_recv[i * A_MATRIX_WIDTH + k] * b_recv[j * B_MATRIX_WIDTH + k];
				c_recv[i * C_MATRIX_WIDTH + j_start + j] = temp;
			}
		}
		if (n != num_of_procs - 1) {
			j_start = n + my_rank == num_of_procs - 1 ? 0 : j_start + b_block_recv;
			MPI_Sendrecv_replace(b_recv, b_block_size*B_MATRIX_WIDTH, MPI_INT, next_rank, 0, prev_rank, 0, MPI_COMM_WORLD, &status);
			MPI_Sendrecv_replace(&b_block_recv, 1, MPI_INT, next_rank, 1, prev_rank, 1, MPI_COMM_WORLD, &status);
		}
	}
	end_wtime = MPI_Wtime() - start_wtime;
	MPI_Gatherv(c_recv, a_block_recv*C_MATRIX_WIDTH, MPI_INT, c, c_sendcounts, c_displs, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&end_wtime, 1, MPI_DOUBLE, all_times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (my_rank == 0) {
		double max_time = all_times[0];
		for(size_t i = 1; i < num_of_procs; ++i) {
			if (all_times[i] > max_time) {
				max_time = all_times[i];
			}
		}
//		printf("%lf ", max_time);
		printf("%d %d %d %d %d %f\n", num_of_procs,
			A_MATRIX_WIDTH, A_MATRIX_HEIGHT, B_MATRIX_WIDTH, B_MATRIX_HEIGHT, max_time);
		free(a);
		free(b);
		free(c);
		free(a_sendcounts);
		free(b_sendcounts);
		free(c_sendcounts);
		free(a_displs);
		free(b_displs);
		free(c_displs);
	}
	free(a_recv);
	free(b_recv);
	free(c_recv);
	MPI_Finalize();
	return 0;
}
