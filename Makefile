CC = g++
OPT = -std=c++11 -O3 -pthread -fopenmp -g
INC = -I src/. -I src/snappy/build/.
LIB = src/snappy/build/libsnappy.a -lz -lm -lboost_unit_test_framework

all: examples

docs: .doxyconf
	doxygen .doxyconf

examples: pagerank

pagerank: examples/pagerank.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

tests: test/test_pagerank.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o tests $^ $(LIB)