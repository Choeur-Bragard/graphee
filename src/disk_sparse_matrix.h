#ifndef GRAPHEE_DISK_SPARSE_MATRIX_H__ 
#define GRAPHEE_DISK_SPARSE_MATRIX_H__ 

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>

#include <gzstream.h>
#include "snappy/snappy.h"

#include "utils.h"
#include "properties.h"
#include "vector.h"

#define SORT_NAME gpe
#define SORT_TYPE uint64_t*
#define SORT_DIM 2
#define SORT_CMP(x, y)  (*(x) < *(y) ? -1 : (*(x) == *(y) ? 0 : 1))
#define SORT_SWAP(x, y) {int i; for (i=0; i<SORT_DIM; i++) {uint64_t __SORT_SWAP = *(x+i); *(x+i) = *(y+i); *(y+i) = __SORT_SWAP;}}
#include "sort/sort.h"

namespace graphee {

template <typename matrixT>
class diskSparseMatrix {
public:
  diskSparseMatrix () {}
  diskSparseMatrix (properties& properties, std::string matrix_name) :
    props(properties), name(matrix_name) 
    {tmpfp = std::vector<std::fstream> (props.nblocks);}

  ~diskSparseMatrix ()
    {close_files ();}

  void load_edgelist (const std::vector<std::string>& filenames, int ftype = utils::GZ, int options = utils::TRANS);

  matrixT& get_block (uint64_t line, uint64_t col);

  using matrixType = matrixT;

  bool empty () const;

private:
  properties& props;

  std::vector<std::fstream> tmpfp;

  std::string name;

  void read_and_split_list (const std::vector<std::string>& filenames, int ftype = utils::GZ);

  /* To be modernize for thread safe new C++11 features */
  static void sort_and_save_list (std::vector<uint64_t>& block, uint64_t nelems,
      std::fstream& fp, std::mutex& mtx);

  /* To be modernize for thread safe new C++11 features */
  static void load_GZ (const std::string filename, std::stringstream& ss,
      std::mutex& read_mtx);

  void diskblock_manager ();
  /* To be modernize for thread safe new C++11 features */
  static void diskblock_builder (diskSparseMatrix<matrixT>* dmat, uint64_t line, uint64_t col,
      std::fstream& tmpfp, size_t& alloc_mem, std::mutex& mtx, std::condition_variable& cond);

  std::string get_block_filename (uint64_t line, uint64_t col);

  void open_files (std::ios_base::openmode mode);
  void close_files ();
};

template <typename matrixT>
void diskSparseMatrix<matrixT>::load_edgelist (const std::vector<std::string>& filenames, int ftype, int options) {
  read_and_split_list (filenames, ftype);
  diskblock_manager ();
}

template <typename matrixT>
matrixT& diskSparseMatrix<matrixT>::get_block (uint64_t line, uint64_t col) {
  std::ostringstream oss;
  oss << "Start to load disk block [" << line << ":" << col << "]";
  print_log (oss.str());

  matrixT mat {};
  mat.load(get_block_filename(line, col));

  return mat;
}

template <typename matrixT>
bool diskSparseMatrix<matrixT>::empty () const {
  return (props.nvertices == 0);
}

/*! Reads the raw edgelist files
 *  Including some in GNU Zip format
 *
 *  It splits the edgelist into blocks in order to make
 *  the D/CSR building faster
 */
template <typename matrixT>
void diskSparseMatrix<matrixT>::read_and_split_list (const std::vector<std::string>& filenames, int ftype) {
  if (2*props.sort_limit*props.nblocks > props.ram_limit) {
    print_error ("To few memory allocated for sorting with respect to the \'ram_limit\' setting");
    exit (-1);
  }

  const uint64_t maxElemsPerSortBlock {props.sort_limit / sizeof(uint64_t)};

  std::vector<std::vector<uint64_t>> edglst_in (props.nblocks, std::vector<uint64_t>(maxElemsPerSortBlock, 0));
  std::vector<std::vector<uint64_t>> edglst_out (props.nblocks, std::vector<uint64_t>(maxElemsPerSortBlock, 0));

  std::vector<uint64_t> edglst_pos (props.nblocks, 0);

  std::vector<std::mutex> write_mtxs (props.nblocks);
  std::mutex read_mtx;

  uint64_t block_id, from_id, to_id;
  open_files (std::ios_base::out | std::ios_base::binary);

  std::stringstream in_stream;
  std::stringstream out_stream;

  load_GZ (filenames[0], in_stream, read_mtx);

  for (uint64_t i = 1; i < filenames.size(); i++) {
    read_mtx.lock();
    std::swap (in_stream, out_stream);
    read_mtx.unlock();

    std::thread gz_reader_thread {load_GZ, filenames[i], std::ref(in_stream), std::ref(read_mtx)};
    gz_reader_thread.detach();

    while (!out_stream.eof()) {
      /* (*ssp_read) >> fromID >> toID; */
      out_stream >> to_id >> from_id;
      if (to_id == from_id) continue; // oriented graph !

      block_id = from_id/props.window + to_id/props.window * props.nslices;

      if (edglst_pos[block_id] < maxElemsPerSortBlock) {
        edglst_in[block_id].at(edglst_pos[block_id]  ) = from_id; 
        edglst_in[block_id].at(edglst_pos[block_id]+1) = to_id;
        edglst_pos[block_id] += 2;

      } else {
        write_mtxs[block_id].lock ();
        std::swap (edglst_in[block_id], edglst_out[block_id]);

        std::thread sort_write_thread {sort_and_save_list, std::ref(edglst_out[block_id]), 
          edglst_pos[block_id], std::ref(tmpfp[block_id]), std::ref(write_mtxs[block_id])};
        sort_write_thread.detach ();

        edglst_pos[block_id] = 0;

        edglst_in[block_id].at(edglst_pos[block_id]  ) = from_id; 
        edglst_in[block_id].at(edglst_pos[block_id]+1) = to_id;
        edglst_pos[block_id] += 2;
      }
    }

    std::ostringstream oss;
    oss << "Sorted file " << filenames[i];
    print_strong_log (oss.str());
  }

  read_mtx.lock();
  std::swap (in_stream, out_stream);
  read_mtx.unlock();

  while (!out_stream.eof()) {
    /* (*ssp_read) >> fromID >> toID; */
    out_stream >> to_id >> from_id;
    if (to_id == from_id) continue; // oriented graph !

    block_id = from_id/props.window + to_id/props.window * props.nslices;

    if (edglst_pos[block_id] < maxElemsPerSortBlock) {
      edglst_in[block_id].at(edglst_pos[block_id]  ) = from_id; 
      edglst_in[block_id].at(edglst_pos[block_id]+1) = to_id;
      edglst_pos[block_id] += 2;

    } else {
      write_mtxs[block_id].lock ();
      std::swap (edglst_in[block_id], edglst_out[block_id]);

      std::thread sort_write_thread {sort_and_save_list, std::ref(edglst_out[block_id]), 
        edglst_pos[block_id], std::ref(tmpfp[block_id]), std::ref(write_mtxs[block_id])};
      sort_write_thread.detach ();

      edglst_pos[block_id] = 0;

      edglst_in[block_id].at(edglst_pos[block_id]  ) = from_id; 
      edglst_in[block_id].at(edglst_pos[block_id]+1) = to_id;
      edglst_pos[block_id] += 2;
    }
  }
  
  std::vector<std::thread> sort_write_threads;
  for (uint64_t i = 0; i < props.nblocks; i++) {
    if (edglst_pos[i] > 0) {
      write_mtxs[i].lock();
      sort_write_threads.push_back(std::thread {sort_and_save_list, std::ref(edglst_in[block_id]), 
        std::ref(edglst_pos[block_id]), std::ref(tmpfp[block_id]), std::ref(write_mtxs[block_id])});
    }
  }

  for (auto& thd : sort_write_threads) {
    if (thd.joinable()) {
      thd.join();
    }
  }

  close_files ();
}

template <typename matrixT>
void diskSparseMatrix<matrixT>::load_GZ (const std::string filename, std::stringstream& ss, std::mutex& read_mtx) {
  igzstream gz_ifp (filename.data());

  if (gz_ifp) {
    read_mtx.lock();
    ss.str(""); // empty string stream
    ss.clear(); // reset the position of EOS

    std::ostringstream oss;
    oss << "Extracting file " << filename;
    print_log (oss.str());

    ss << gz_ifp.rdbuf();
    gz_ifp.close();
    read_mtx.unlock();
  } else {
    std::ostringstream oss;
    oss << "Cannot open file " << filename << ", exiting...";
    print_error (oss.str());
  }
}

template <typename matrixT>
void diskSparseMatrix<matrixT>::sort_and_save_list (std::vector<uint64_t>& block, uint64_t nelems,
    std::fstream& ofp, std::mutex& mtx) {

  if (nelems%2 != 0) {
    print_error ("Wrong number of edges");
    return;
  }

  size_t nedges = nelems/2;
  uint64_t** edgeptrs = new uint64_t* [nedges];
  for (uint64_t i = 0; i < nedges; i++) {
    edgeptrs[i] = &(block[2*i]);
  }

  gpe_heap_sort(edgeptrs, nedges);

  ofp.write(reinterpret_cast<const char*> (block.data()), sizeof(uint64_t)*nelems);

  mtx.unlock();

  delete[] edgeptrs;
}

template <typename matrixT>
void diskSparseMatrix<matrixT>::diskblock_manager () {
  std::vector<std::thread> diskblock_threads;

  open_files (std::ios_base::in | std::ios_base::binary);

  std::mutex mtx;
  std::condition_variable cond;

  size_t alloc_mem {0};
  size_t ram_limit {props.ram_limit};
  uint64_t block_id;

  for (uint64_t line = 0; line < props.nslices; line++) {
    for (uint64_t col = 0; col < props.nslices; col++) {
      block_id = line + col*props.nslices;

      diskblock_threads.push_back(std::thread {diskblock_builder, this, line, col, std::ref(tmpfp[block_id]), 
          std::ref(alloc_mem), std::ref(mtx), std::ref(cond)});
    }
  }

  for (auto& thd : diskblock_threads) {
    if (thd.joinable ()) {
      thd.join ();
    }
  }

  close_files ();
}

template <typename matrixT>
void diskSparseMatrix<matrixT>::diskblock_builder (diskSparseMatrix<matrixT>* dmat, uint64_t line, uint64_t col,
      std::fstream& tmpfp, size_t& alloc_mem, std::mutex& mtx, std::condition_variable& cond) {
  properties& props {dmat->props};

  tmpfp.seekg (0, tmpfp.end);
  size_t filelen = tmpfp.tellg ();
  tmpfp.seekg (0, tmpfp.beg);

  uint64_t nnz = filelen/(2*sizeof(uint64_t));
  uint64_t nsections = filelen/props.sort_limit + (filelen%props.sort_limit == 0 ? 0 : 1);

  size_t alloc_needs {0};
  if (typeid (typename matrixT::valueType) == typeid (bool)) {
    alloc_needs = (props.window+1+nnz)*sizeof(uint64_t);
  } else {
    alloc_needs = (props.window+1+nnz)*sizeof(uint64_t) + nnz*sizeof(typename matrixT::valueType);
  }

  if (alloc_needs > props.ram_limit) {
    mtx.lock();
    std::ostringstream err;
    err << "Disk block [" << line << ";" << col << "] needs " << alloc_needs/(1UL << 30) << "GB";
    print_error (err.str());
    err.str("");
    err << "which is more memory than \'ram_limit\' " << props.ram_limit/(1UL << 30) << "GB";
    print_error (err.str());
    mtx.unlock();
    return; // not exit(-1); because we are within a thread
  } else {
    mtx.lock();
    std::ostringstream log;
    log << "Starting disk block conversion [" << line << ";" << col << "]";
    print_log (log.str());
    mtx.unlock();
  }

  std::unique_lock<std::mutex> mlock(mtx);
  while (alloc_mem + alloc_needs > props.ram_limit) {
    cond.wait(mlock);
  }
  alloc_mem += alloc_needs;
  mlock.unlock();

  matrixT mat {props, props.window, props.window, nnz};

  std::vector<uint64_t> offsets {nsections};
  for (uint64_t i = 0; i < nsections; i++) {
    offsets[i] = i*props.sort_limit;
  }
  std::vector<uint64_t> ids {nsections, line*props.window};

  uint64_t currentid {line*props.window};
  uint64_t minid {(line+1)*props.window};

  uint64_t edge[2];
  uint64_t offl = line*props.window;
  uint64_t offc = col*props.window;

  while (currentid < (line+1)*props.window) {
    minid = (line+1)*props.window;

    for (uint64_t sec = 0; sec < nsections; sec++) {
      if (ids[sec] == currentid && std::min(filelen, (sec+1)*props.sort_limit)) {
        if (tmpfp.tellg() != offsets[sec]) {
          tmpfp.seekg(offsets[sec], tmpfp.beg);
        }

        while (offsets[sec] < std::min(filelen, (sec+1)*props.sort_limit)) {
          tmpfp.read((char*)edge, 2*sizeof(uint64_t));

          if (edge[0] == currentid) {
            mat.fill (edge[0], edge[1], true);
            offsets[sec] += 2*sizeof(uint64_t);
          } else {
            ids[sec] = edge[0];
            break;
          }
        }
      }

      if (ids[sec] > currentid) {
        minid = std::min(minid, ids[sec]);
      }
    }
    currentid = minid;
  }

  if (!mat.verify()) {
    mtx.lock();
    std::ostringstream err;
    err << "Block [" << line << ";" << col << "] conversion to \'" << mat.matrixType << "\' failed !";
    print_error (err.str());
    mtx.unlock();
    return;
  } else {
    mtx.lock();
    std::ostringstream log;
    log << "Block [" << line << ";" << col << "] conversion to \'" << mat.matrixType << "\' succeed !";
    print_log (log.str());
    mtx.unlock();
  }

  mat.save(dmat->get_block_filename(line, col), utils::SNAPPY);

  mtx.lock ();
  alloc_mem -= alloc_needs;
  mtx.unlock ();

  cond.notify_one ();
}

template <typename matrixT>
std::string diskSparseMatrix<matrixT>::get_block_filename (uint64_t line, uint64_t col) {
  std::ostringstream matrixname;
  matrixname << props.name << "_" << name << "_dmatblk_" << line << "_" << col << ".gpe";
  return matrixname.str();
}

template <typename matrixT>
void diskSparseMatrix<matrixT>::open_files (std::ios_base::openmode mode) {
  std::ostringstream blockname;

  for (uint64_t line = 0; line < props.nslices; line++) {
    for (uint64_t col = 0; col < props.nslices; col++) {
      uint64_t bid = line + col*props.nslices;
      blockname.str("");
      blockname << props.name << "_" << name << "_tmpblk_" << line << "_" << col << ".gpe";
      tmpfp[bid].open(blockname.str(), mode);

      if (!tmpfp[bid].is_open()) {
        std::ostringstream err;
        err << "Could not open file: " << blockname.str();
        print_error (err.str());
        exit(-1);
      }
    }
  }
}

template <typename matrixT>
void diskSparseMatrix<matrixT>::close_files () {
  for (auto& fp : tmpfp) {
    if (fp.is_open()) {
      fp.close();
    }
  }
}

} // namespace graphee

#endif // GRAPHEE_DISK_SPARSE_MATRIX_H__
