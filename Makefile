CC = g++
#OPT = -std=c++11 -O3 -pthread
OPT = -std=c++11 -g -pthread
INC = -I src/.
LIB = src/snappy/build/libsnappy.a -lgzstream -lz -lm

examples: split load vector matvecprod pagerank

pagerank: examples/pagerank.cc src/utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

matvecprod: examples/matvecprod.cc src/utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

vector: examples/vector.cc src/utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

split: examples/split.cc src/utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

load: examples/load.cc src/utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)
