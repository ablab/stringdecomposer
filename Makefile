export ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BUILD_DIR = $(ROOT_DIR)/build
BIN_DIR = $(BUILD_DIR)/bin
SRC_DIR = $(ROOT_DIR)/src


all:
	mkdir -p $(BIN_DIR)
	${CXX} -o $(BIN_DIR)/dp $(SRC_DIR)/main.cpp -fopenmp --std=c++11 -O2 -Wall -Wextra -pedantic -Wshadow -Wfloat-equal -fsanitize=address

clean:
	-rm -rf $(BUILD_DIR)
