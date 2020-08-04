export ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
SD_DIR = $(ROOT_DIR)/sd
BIN_DIR = ${SD_DIR}/bin
SRC_DIR = ${SD_DIR}/src


all: sd
	mkdir -p $(BIN_DIR)
	${CXX} -o $(BIN_DIR)/dp $(SRC_DIR)/main.cpp -fopenmp --std=c++11

clean:
	-rm -rf $(BIN_DIR)
	-rm -rf build
