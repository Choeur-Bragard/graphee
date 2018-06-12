# Graphee :chart_with_upwards_trend:
Graphee is a light-weight code written in C++ to compute graph properties, such as:
PageRank, TrustRank, harmonic centrality, etc...

## What's the purpose of Graphee?
In order to improve search engine results, one has compute some values from the graph of the web.
Actually, some famous works are known like the
[*PageRank*](http://infolab.stanford.edu/pub/papers/google.pdf) or the
[*TrustRank*](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf) algorithms. They were
booth indispensable to make their search engines pertinent and spam proof.

However, with time the graph of the web has become wider and wider. It clearly contains more than a
dozen billion vertices and at least ten times more edges. Some, has chosen to persist in full-RAM
computations and using distributed framework as *Hadoop* or *Spark*. While, some other approaches as
[*GraphChi*](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf) or
[*m-flash*](https://www.cc.gatech.edu/~dchau/papers/16-pkdd-mflash.pdf) are using serialization methods,
but both limited by `uint32_t` for the vertex ids.

### Graphee's innovations
Thus, we propose to develop software able to handle `uint64` and still performing calculations
on a single laptop (sometimes with an external hard-drive for wide graphs).

The graph is converted into an adjacency matrix and then saved with the fully-tested (since 70's)
[Compressed Sparse Row (CSR) matrix format](https://en.wikipedia.org/wiki/Sparse_matrix).
The serialization is done by dividing vectors into slices, and matrices into blocks, with sizes
not exceeding the RAM-limit of the computer.

However, nowadays the graph of the web is so wide that one cannot restart calculations from scratch.
The preferred solution would be using delta files, saving only differences between two graph
snapshots. But the CSR format is unfriendly with insertions and deletions. So we propose
to implement the [*Dynamic CSR*](https://thomas.gilray.org/pdf/dynamic-csr.pdf) matrix format,
which presents high speedups for the matrix modifications.

Also a binding of Graphee with `Python` is under development !

## Efficiency
Technically we use `pthread` and `OpenMP` technologies for the parallel calculations on CPUs.
The full support or `CUDA` standards is one of the major properties to be included soon !

## Examples & functionalities
**NOTE: For now on, the code contains all the needed tools for the PageRank.
But I need to write a small script to DL few test data. I'm on my way...**

## Contributing
Please first read `CONRTIBUTING.md` and propose what you want or you can fix or add functionalities detailed
within `TODO.md`.

For any questions, comments, or collaborations, please use: *n.martin [at] qwantresearch [dot] com*.

## References
- [The Anatomy of a Large-Scale Hypertextual Web Search Engine, *S. Brin and L. Page*](http://infolab.stanford.edu/pub/papers/google.pdf)
- [Web Spam Detection with Anti-Trust Rank, *V. Krishnan and R. Raj*](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf)
- [GraphChi: Large-Scale Graph Computation on Just a PC, *A. Kyrola, G. Blelloch, and C. Guestrin*](http://i.stanford.edu/~kvijay/krishnan-raj-airweb06.pdf)
- [M-Flash: Fast Billion-Scale Graph Computation Using a Bimodal Block Processing Model, *H. Gualdron* et al.](https://www.cc.gatech.edu/~dchau/papers/16-pkdd-mflash.pdf)
- [Dynamic Sparse-Matrix Allocation on GPUs, *J. King, T. Gilray, R.M. Kirby and M. Might*](https://thomas.gilray.org/pdf/dynamic-csr.pdf)
