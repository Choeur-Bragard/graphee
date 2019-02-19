#define BOOST_TEST_MODULE example
#include <boost/test/included/unit_test.hpp>
#include "pagerank.hpp"
#include "vector.hpp"
#include <cstdio>

void clean_pagerank_files(const graphee::Properties& props){
    for(int i = 0; i < props.nslices; i++){
        for(int j = 0; j < props.nslices; j++){
            std::remove((props.name+"_adj_dmatblk_"+std::to_string(i)+"_"+std::to_string(j)+".gpe").c_str());
            std::remove((props.name+"_adj_tmpblk_"+std::to_string(i)+"_"+std::to_string(j)+".gpe").c_str());
        }
        std::remove((props.name+"_ob_dvecslc_"+std::to_string(i)+".gpe").c_str());
        std::remove((props.name+"_prp1_dvecslc_"+std::to_string(i)+".gpe").c_str());
        std::remove((props.name+"_pr_dvecslc_"+std::to_string(i)+".gpe").c_str());
    }
}

BOOST_AUTO_TEST_CASE( test_dense_web )
/* Compare with void free_test_function() */
{
  graphee::Properties props(
      std::string("test_pagerank_web"),            // name of your graph
      325729,                              // number of nodes
      5,                         // number of slices
      4,                              // number of threads
      5 * graphee::Properties::GB,    // max RAM value
      32 * graphee::Properties::MB); // max size of sorting vector
  
  graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR> adjacency_matrix(
      &props, // graph properties
      "adj");
  
  std::vector<std::string> filenames;
  filenames.push_back("test/ressources/web-NotreDame.txt.gz");
  adjacency_matrix.load_edgelist(filenames);
  
    graphee::Pagerank<graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR>>
        pagerank(&props,            // graph properties
                 &adjacency_matrix, // give the adress to the adjacency matrix
                 0.85); // damping factor of the Pagerank (original value)

    /**
     * Compute the Pagerank with 100 iterations
     */
    pagerank.compute_pagerank(100);
    long n=0;
    float score_sum=0;
    graphee::Vector<float> vec(&props);
    for(uint64_t slice_i=0; slice_i<props.nslices; slice_i++){
        vec.load("test_pagerank_web_pr_dvecslc_"+std::to_string(slice_i)+".gpe");
        for(float score : vec){
            score_sum+=score;
            n++;
        }
    }
    BOOST_CHECK(abs(score_sum-1.)<0.001);
    std::cout<<"SCORE SUM : "<<score_sum<<std::endl;
    clean_pagerank_files(props);
    std::cout<<"NVERTEX : "<<n<<std::endl;
}

BOOST_AUTO_TEST_CASE( test_smallGraph )
/* Compare with void free_test_function() */
{
  graphee::Properties props(
      std::string("test_smallGraph"),            // name of your graph
      6,                              // number of nodes
      2,                         // number of slices
      1,                              // number of threads
      5 * graphee::Properties::GB,    // max RAM value
      128 * graphee::Properties::MB); // max size of sorting vector
  
  graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR> adjacency_matrix(
      &props, // graph properties
      "adj");
  
  std::vector<std::string> filenames;
  filenames.push_back("test/ressources/test_smallGraph.txt.gz");
  adjacency_matrix.load_edgelist(filenames);
  
    graphee::Pagerank<graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR>>
        pagerank(&props,            // graph properties
                 &adjacency_matrix, // give the adress to the adjacency matrix
                 0.85); // damping factor of the Pagerank (original value)

    /**
     * Compute the Pagerank with 10 iterations
     */
    pagerank.compute_pagerank(10);
    long n=0;
    float score_sum=0;
    graphee::Vector<float> vec(&props);
    float expected[] = {0.21495,0.15189,0.03953,0.26713,0.22387,0.10260};
    for(uint64_t slice_i=0; slice_i<props.nslices; slice_i++){
        vec.load("test_smallGraph_pr_dvecslc_"+std::to_string(slice_i)+".gpe");
        for(float score : vec){
            std::cout<<n<<"\t"<<score<<std::endl;
            score_sum+=score;
            BOOST_CHECK(abs(score-expected[n])<0.00001);
            n++;
        }
    }
    std::cout<<"SCORE SUM : "<<score_sum<<std::endl;
    clean_pagerank_files(props);
}


BOOST_AUTO_TEST_CASE( test_smallGraph_sum_columns )
/* Compare with void free_test_function() */
{
  graphee::Properties props(
      std::string("test_sum_columns"),            // name of your graph
      1000,                              // number of nodes
      20,                         // number of slices
      4,                              // number of threads
      10 * graphee::Properties::GB,    // max RAM value
      1 * graphee::Properties::MB); // max size of sorting vector
  
  graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR> adjacency_matrix(
      &props, // graph properties
      "adj");
  
  std::vector<std::string> filenames;
  filenames.push_back("test/ressources/test_sum_columns.txt.gz");
  adjacency_matrix.load_edgelist(filenames);

    graphee::DiskVector<graphee::Vector<float>> out_bounds = std::move(graphee::DiskVector<graphee::Vector<float>>(&props, "ob", 0.));
out_bounds.dmat_columns_sum(adjacency_matrix);
float n = 0;
for(int slice = 0 ; slice < props.nslices ; slice++){
      graphee::Vector<float>& vec = out_bounds.get_slice(slice);
      for(int i = 0;i<vec.get_lines();i++){
        BOOST_CHECK_EQUAL(vec.at(i),n);
        n+=1;
      }
    }
    clean_pagerank_files(props);
}