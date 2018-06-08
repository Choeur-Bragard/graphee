CC = g++
OPT = -std=c++11 -O3 -pthread
#OPT = -std=c++11 -g -pthread
INC = -I src/.
LIB = src/snappy/build/libsnappy.a -lgzstream -lz -lm

all: examples doc

doc: .doxyconf
	doxygen .doxyconf

examples: split

split: examples/split.cc src/utils.cc
	$(CC) $(OPT) $(INC) -o examples/$@ $^ $(LIB)
