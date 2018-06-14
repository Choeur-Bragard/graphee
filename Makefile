CC = g++
OPT = -std=c++11 -O3 -pthread
#OPT = -std=c++11 -g -pthread
INC = -I src/.
LIB = src/snappy/build/libsnappy.a -lz -lm

all: examples

docs: .doxyconf
	doxygen .doxyconf

examples: split matvecprod pagerank edgelist

split: examples/split.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

matvecprod: examples/matvecprod.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

pagerank: examples/pagerank.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

edgelist: examples/edgelist.cpp src/edgelist.cpp src/utils.cpp
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)
