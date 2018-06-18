CC = g++
OPT = -std=c++11 -O3 -pthread -fopenmp
INC = -I src/. -I src/snappy/build/.
LIB = src/snappy/build/libsnappy.a -lz -lm

all: examples

docs: .doxyconf
	doxygen .doxyconf

examples: pagerank

pagerank: examples/pagerank.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)
