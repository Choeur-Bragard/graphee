# Graphee :chart_with_upwards_trend:
Graphee is a light-weight code written in C++ to compute graphs properties.

## What to do?
The graph study is essential within the web science. Actually, the output of search engine could be clearly improved by a good knowledge of its structure. Two of the flagship studies are the PageRank and the TrustRank. 
Nowadays the extraction of some metrics from the graph of the web is essential. However, the available frameworks are either limited by `int32` and about 4.2B vertices (GraphChi and m-flash), or requires a huge amount of RAM and have to be launched over a Hadoop cluster.

We propose to extend the fields open by GraphChi and m-flash in order to make a more powerful tool for graph studies. We propose to use compressed sparse row (CSR) format to save the graph, and implement the Dynamic CSR (DCSR) which allows the update of the graph. Nowadays, some snapshots of the web contains more than 180B edges. In such cases, one cannot restart a graph calculation from scratch. The solution is to treat delta between two version.

## Technologies
For the moment we restrict ourselves on `POSIX Threads` and `OpenMP` for the parallelization stuff. However, we are preparing a full `CUDA` compatible framework.

## Examples
The basic example concerns the calculations of a PageRank of a 98M vertices graph (from Common Crawl).

## Joining us?
Please before read `CONTRIBUTING.md` and if you don't know what to do the `TODO.md`. Then propose some pull-request !

## References
