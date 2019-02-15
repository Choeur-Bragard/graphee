#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <omp.h>

#include "graphee.hpp"

using namespace graphee;

int main (int argc, char** argv)
{
  /**
   * Necessary for linear algebra operations
   */
  omp_set_nested(1);

  /**
   * Declare the properties of the graph
   */
  graphee::Properties
  props (std::string("cc18q1"), // name of your graph
         98, // number of links in million
         4, // number of slices
         8, // number of threads
         10*graphee::Properties::GB, // max RAM value
         256*graphee::Properties::MB); // max size of sorting vector

  /**
   * Read the list of files, which contains all
   * the files necessary for the graph
   */
  std::ifstream filelist(argv[1]); // open file
  std::vector<std::string> filenames;
  std::string filename;
  while (filelist.good()) // read the file until EOF
  {
    filelist >> filename;
    if(filename.size()>0 && filename.at(0)!='#')
    {
      filenames.push_back(filename);
    }
  }
  filelist.close();

  /**
   * Declare the adjacency sparse matrix
   */
  graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR>
  adjacency_matrix (&props, // graph properties
                    "adj"); // name of the `diskmatrix`

  /**
   * Load from the raw files the structure of
   * the graph
   */
  adjacency_matrix.load_edgelist(filenames);

  /**
   * Declare the Pagerank object
   */
  graphee::Pagerank<graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR>>
      pagerank(&props, //graph properties
               &adjacency_matrix, // give the adress to the adjacency matrix
               0.85); // damping factor of the Pagerank (original value)

  /**
   * Compute the Pagerank with 10 iterations
   */
  pagerank.compute_pagerank(10);

  return 0;
}
