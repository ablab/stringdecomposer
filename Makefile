export ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
SD_DIR = $(ROOT_DIR)/stringdecomposer
BUILD_DIR = $(SD_DIR)/build
BIN_DIR = $(BUILD_DIR)/bin
SRC_DIR = $(SD_DIR)/src

TEST_QUERY = $(SD_DIR)/test_data/read.fa
TEST_MONOMERS = $(SD_DIR)/test_data/DXZ1_star_monomers.fa
TEST_OUTDIR = $(SD_DIR)/test_data
TEST_REFERENCE = $(SD_DIR)/test_data/final_decomposition_fc89af8.tsv

build:
	mkdir -p $(BIN_DIR)
	${CXX} -o $(BIN_DIR)/dp $(SRC_DIR)/main.cpp $(SRC_DIR)/edlib.cpp -fopenmp --std=c++11 -O2 -Wall -Wextra -pedantic -Wshadow -Wfloat-equal

test_launch: build
	bin/stringdecomposer $(TEST_QUERY) $(TEST_MONOMERS) -o $(TEST_OUTDIR) --second-best
	grep -q "Thank you for using StringDecomposer!" $(TEST_OUTDIR)/stringdecomposer.log
	diff -q $(TEST_REFERENCE) $(TEST_OUTDIR)/final_decomposition.tsv

install: build
	python setup.py install --record install_footprint.txt

test_launch_install: install
	stringdecomposer $(TEST_QUERY) $(TEST_MONOMERS) -o $(TEST_OUTDIR) --second-best
	grep -q "Thank you for using StringDecomposer!" $(TEST_OUTDIR)/stringdecomposer.log
	diff -q $(TEST_REFERENCE) $(TEST_OUTDIR)/final_decomposition.tsv

clean:
	-rm -rf $(BUILD_DIR)
	-rm -rf $(SD_DIR)/test_data/final_decomposition_alt.tsv
	-rm -rf $(SD_DIR)/test_data/final_decomposition_raw.tsv
	-rm -rf $(SD_DIR)/test_data/final_decomposition.tsv
	-rm -rf $(SD_DIR)/test_data/stringdecomposer.log
	-rm -rf StringDecomposer.egg-info dist build

uninstall:
	@if [ -f install_footprint.txt ]; then\
		echo "removing install footprint from install_footprint.txt";\
		cat install_footprint.txt | xargs rm -rf;\
		rm -rf install_footprint.txt;\
	fi
