ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BIN_DIR = ${ROOT_DIR}/bin

all:
	mkdir -p $(BIN_DIR)
	g++ -o $(BIN_DIR)/dp src/main.cpp -fopenmp --std=c++11

clean:
	-rm -r ${BIN_DIR}
