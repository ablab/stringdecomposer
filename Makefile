export ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BUILD_DIR = $(ROOT_DIR)/build
BIN_DIR = $(BUILD_DIR)/bin
SRC_DIR = $(ROOT_DIR)/src


all:
	mkdir -p $(BIN_DIR)
	${CXX} -o $(BIN_DIR)/dp $(SRC_DIR)/main.cpp -fopenmp --std=c++11 -O2 -Wall -Wextra -pedantic -Wshadow -Wfloat-equal -fsanitize=address

test_launch: all
	python string_decomposer.py ./test_data/read.fa ./test_data/DXZ1_star_monomers.fa -o test_data
	grep -q "Thank you for using StringDecomposer!" test_data/string_decomposer.log
	diff -q ./test_data/final_decomposition_fc89af8.tsv ./test_data/final_decomposition.tsv

clean:
	-rm -rf $(BUILD_DIR)
	-rm -rf test_data/final_decomposition{_alt.tsv,_raw.tsv,.tsv} test_data/string_decomposer.log
