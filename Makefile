export ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BUILD_DIR = $(ROOT_DIR)/build
BIN_DIR = $(BUILD_DIR)/bin
SRC_DIR = $(ROOT_DIR)/src

TEST_QUERY = "test_data/read.fa"
TEST_MONOMERS = "test_data/DXZ1_star_monomers.fa"
TEST_OUTDIR = "test_data"
TEST_REFERENCE = "test_data/final_decomposition_fc89af8.tsv"

build:
	mkdir -p $(BIN_DIR)
	${CXX} -o $(BIN_DIR)/dp $(SRC_DIR)/main.cpp -fopenmp --std=c++11 -O2 -Wall -Wextra -pedantic -Wshadow -Wfloat-equal -fsanitize=address

test_launch: build
	python stringdecomposer.py $(TEST_QUERY) $(TEST_MONOMERS) -o $(TEST_OUTDIR)
	grep -q "Thank you for using StringDecomposer!" $(TEST_OUTDIR)/stringdecomposer.log
	diff -q $(TEST_REFERENCE) $(TEST_OUTDIR)/final_decomposition.tsv

install: build
	python setup.py install --record install_footprint.txt

test_launch_install: install
	stringdecomposer $(TEST_QUERY) $(TEST_MONOMERS) -o $(TEST_OUTDIR)
	grep -q "Thank you for using StringDecomposer!" $(TEST_OUTDIR)/stringdecomposer.log
	diff -q $(TEST_REFERENCE) $(TEST_OUTDIR)/final_decomposition.tsv


clean:
	-rm -rf $(BUILD_DIR)
	-rm -rf test_data/final_decomposition{_alt.tsv,_raw.tsv,.tsv} test_data/stringdecomposer.log
	@if [ -f install_footprint.txt ]; then\
		echo "removing install footprint from install_footprint.txt";\
		cat install_footprint.txt | xargs rm -rf;\
		rm -rf install_footprint.txt;\
	fi
	-rm -rf StringDecomposer.egg-info dist
	python setup.py clean
