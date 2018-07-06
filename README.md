# Graphee :chart_with_upwards_trend:
Graphee is a lightweight piece of code written in C++ to compute graph properties, such as:
PageRank, TrustRank, harmonic centrality, etc...

## What's the purpose of Graphee?
In order to improve search engine results, one has compute some values from the graph of the web.
Actually, some famous works are known like the
[PageRank](http://infolab.stanford.edu/pub/papers/google.pdf) [\[1\]](#References) or the
[TrustRank](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf) [\[2\]](#References) algorithms. They were
both indispensable to make their search engines pertinent and spam-proof.

However, with time, the graph of the web has become wider and wider. It clearly contains more than a
dozen billion vertices and at least ten times more edges. Some, has chosen to persist in full-RAM
computations and using distributed framework as *Hadoop* or *Spark*. While, some other approaches as
[GraphChi](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf) [\[3;4\]](#References) or
[m-flash](https://www.cc.gatech.edu/~dchau/papers/16-pkdd-mflash.pdf) [\[5;6\]](#References) are using serialization methods,
but both limited by `uint32_t` for the vertex ids.

### Graphee's innovations
Thus, we propose to develop software able to handle `uint64` and still performing calculations
on a single laptop (sometimes with an external hard-drive for wide graphs).

The graph is converted into an adjacency matrix and then saved with the fully-tested (since 70's)
[Compressed Sparse Row (CSR) matrix format](https://en.wikipedia.org/wiki/Sparse_matrix).
The serialization is done by dividing vectors into slices, and matrices into blocks, with sizes
not exceeding the RAM limit of the computer.

However, nowadays the graph of the web is so wide that one cannot restart calculations from scratch.
The preferred solution would be using delta files, saving only differences between two graph
snapshots. But the CSR format is unfriendly with insertions and deletions. So we propose
to implement the [Dynamic CSR](https://thomas.gilray.org/pdf/dynamic-csr.pdf) [\[7\]](#References) matrix format,
which presents high speedups for the matrix modifications.

**Also a binding of Graphee with `Python` is under development !**

## Efficiency
Technically we use `pthread` and `OpenMP` technologies for the parallel calculations on CPUs.
The full support or `CUDA` standards is one of the major properties to be included soon !

### Benchmarks
We are presenting the results for the Pagerank computaton of the
domain-level graph of the [CommonCrawl](https://commoncrawl.org/2018/05/webgraphs-feb-mar-apr-2018/):
- 98 million vertices,
- 1.5 billion edges.

The graph prepossessing (RAW to CSR) and the 10 iterations of the Pagerank,
runs in **16 minutes and 11 seconds** on a Dell XPS laptop made of an Intel Core i7-7700HQ @ 2.8GHz x 8, 16GB RAM, 512GB SSD.
(We limited the RAM usage to 10GB to avoid freezing because of the RAM usage of other apps)

## Examples & functionalities

### Preliminaries
- Verify that `zlib` is installed on your computer. If not follow the link [ZLib](http://zlib.net).
- Clone `graphee` and synchronize the `submodules` necessary for the code, with the following command:
```
$ git clone --recursive https://github.com/QwantResearch/graphee.git
```
- Compile the customized version of `Snappy` which supports files over `4.2 GB` (again the 2018-problem of `int32`):
```
$ cd src/snappy
$ mkdir -p build
$ cd build && cmake ..
$ make
```

### Pagerank example
Compile the [Pagerank example](examples/pagerank.cpp) and download the data:
```
$ make pagerank
$ cd examples/
$ chmod u+x load_data.sh
$ ./load_data.sh
```
Then launch it ! :tada:
```
$ ./pagerank filelist
```

## Documentation
The documentation is realized with [Doxygen](https://www.stack.nl/~dimitri/doxygen/), create it:
```
$ mkdir -p docs
$ make docs
```

## Contributing
Please first read [`CONTRIBUTING.md`](CONTRIBUTING.md) and propose what you want or you can fix or add functionalities detailed
within [`TODO.md`](TODO.md).

For any questions, comments, or collaborations, please use: **n.martin [at] qwantresearch [dot] com** or also [@stdthread](https://www.twitter.com/stdthread) on Twitter.

## References
1. [The Anatomy of a Large-Scale Hypertextual Web Search Engine, *S. Brin and L. Page*](http://infolab.stanford.edu/pub/papers/google.pdf)
2. [Web Spam Detection with Anti-Trust Rank, *V. Krishnan and R. Raj*](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf)
3. [GraphChi: Large-Scale Graph Computation on Just a PC, *A. Kyrola, G. Blelloch, and C. Guestrin*](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf)
4. [Sources: GraphChi/graphchi-cpp](https://github.com/GraphChi/graphchi-cpp)
5. [M-Flash: Fast Billion-Scale Graph Computation Using a Bimodal Block Processing Model, *H. Gualdron* et al.](https://www.cc.gatech.edu/~dchau/papers/16-pkdd-mflash.pdf)
6. [Sources: M-Flash/m-flash-cpp](https://github.com/M-Flash/m-flash-cpp)
7. [Dynamic Sparse-Matrix Allocation on GPUs, *J. King, T. Gilray, R.M. Kirby and M. Might*](https://thomas.gilray.org/pdf/dynamic-csr.pdf)
