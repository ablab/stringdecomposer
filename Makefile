export ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BIN_DIR = $(ROOT_DIR)/build/bin
SRC_DIR = $(ROOT_DIR)/src


all:
	mkdir -p $(BIN_DIR)
	${CXX} -o $(BIN_DIR)/dp $(SRC_DIR)/main.cpp -fopenmp --std=c++11

clean:
	-rm -rf $(BIN_DIR)
