CC = g++
OPT = -std=c++11 -O3 -pthread
INC = -I src/.
LIB = -lgzstream -lsnappy -lz -lm

examples: split load vector

vector: examples/vector.cc src/gpe_utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

split: examples/split.cc src/gpe_utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)

load: examples/load.cc src/gpe_utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)
