CXXFLAGS += -Wall -g -fopenmp

SRC_DIR = ./src
BIN_DIR = ./bin

all: directories omp mpi

omp: $(BIN_DIR)/tape_matrix_mul_omp

mpi: $(BIN_DIR)/tape_matrix_mul_mpi

$(BIN_DIR)/tape_matrix_mul_omp: directories $(SRC_DIR)/tape_matrix_mul_omp.c
	$(CXX) $(CXXFLAGS) $(lastword $^) -o $@
	#mpixlc_r -qsmp=omp $(lastword $^) -o $@

$(BIN_DIR)/tape_matrix_mul_mpi: directories $(SRC_DIR)/tape_matrix_mul_mpi.c
	mpiCC $(lastword $^) -o $@
	#mpixlc_r -qsmp=omp $(lastword $^) -o $@

directories:
	mkdir -p $(BIN_DIR)

clean:
	rm -rvf $(BIN_DIR)
	
