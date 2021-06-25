export ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BUILD_DIR = $(ROOT_DIR)/build
BIN_DIR = $(BUILD_DIR)/bin
SRC_DIR = $(ROOT_DIR)/src


build:
	mkdir -p $(BIN_DIR)
	${CXX} -o $(BIN_DIR)/dp $(SRC_DIR)/main.cpp -fopenmp --std=c++11 -O2 -Wall -Wextra -pedantic -Wshadow -Wfloat-equal -fsanitize=address

test_launch: build
	python string_decomposer.py ./test_data/read.fa ./test_data/DXZ1_star_monomers.fa -o test_data
	grep -q "Thank you for using StringDecomposer!" test_data/string_decomposer.log
	diff -q ./test_data/final_decomposition_fc89af8.tsv ./test_data/final_decomposition.tsv

install: build
	python setup.py install --record install_footprint.txt

clean:
	-rm -rf $(BUILD_DIR)
	-rm -rf test_data/final_decomposition{_alt.tsv,_raw.tsv,.tsv} test_data/string_decomposer.log
	@if [ -f install_footprint.txt ]; then\
		echo "removing install footprint from install_footprint.txt";\
		cat install_footprint.txt | xargs rm -rf;\
		rm -rf install_footprint.txt;\
	fi
