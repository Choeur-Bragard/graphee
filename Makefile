CC = g++
OPT = -std=c++11 -O3 -pthread
INC = -I src/.
LIB = src/snappy/build/libsnappy.a -lz -lm

all: examples

docs: .doxyconf
	doxygen .doxyconf

examples: pagerank

pagerank: examples/pagerank.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)
