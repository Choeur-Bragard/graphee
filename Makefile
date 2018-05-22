CC = g++
OPT = -std=c++11 -O3 -pthread
INC = -I src/.
LIB = -lgzstream -lsnappy

examples: split

split: examples/split.cc src/gpe_utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)
