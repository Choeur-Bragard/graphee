#ifndef GPE_DISKMAT_H
#define GPE_DISKMAT_H

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
#include <chrono>

#include <gzstream.h>
#include <snappy.h> // SNAPPY Compression Google Inc.

#include "gpe_bsmat_csr.h"
#include "gpe_props.h"

#define SORT_NAME gpe
#define SORT_TYPE uint64_t*
#define SORT_DIM 2
#define SORT_CMP(x, y)  (*(x) < *(y) ? -1 : (*(x) == *(y) ? 0 : 1))
#define SORT_SWAP(x, y) {int i; for (i=0; i<SORT_DIM; i++) {uint64_t __SORT_SWAP = *(x+i); *(x+i) = *(y+i); *(y+i) = __SORT_SWAP;}}
#include "sort/sort.h"

namespace graphee {

template <class gpe_mat_t, class idx_t>
class gpe_diskmat {
public:
  gpe_diskmat (gpe_props in_props);
  ~gpe_diskmat ();

  void load_edgelist (const std::vector<std::string>& filenames, int ftype, int options);

  void preLoadLine (const uint64_t lineID);
  void preLoadCol (const uint64_t colID);
  void getCSRBlock (const uint64_t blockID, uint64_t& m, uint64_t** ia, uint64_t& nnz, uint64_t** ja);

  enum {PLAIN, GZ};
  static const int DIRECT = 0x00000001; 
  static const int TRANS  = 0x00000010; 
  static const int UO     = 0x00000100; 

private:
  std::fstream *tmpfp;

  gpe_props props;

  void read_and_split_list (const std::vector<std::string>& filenames, int ftype);

  static void sort_and_save_list (uint64_t* block, uint64_t nelems,
      std::fstream* fp, std::mutex* mtx);

  static void load_GZ (const std::string filename, std::stringstream* ss,
      std::mutex* read_mtx);

  void csr_manager ();
  static void csr_builder (uint64_t m, uint64_t line, uint64_t col, std::fstream* tmpfp,
      gpe_props* props, size_t* alloc_mem, std::mutex* mtx, std::condition_variable* cond);

  void open_tmp_blocks (std::ios_base::openmode mode);
  void close_tmp_blocks ();
};

template <class gpe_mat_t, class idx_t>
gpe_diskmat<gpe_mat_t, idx_t>::gpe_diskmat (gpe_props in_props) {
  props = in_props;
  tmpfp = new std::fstream [props.nblocks];
}

template <class gpe_mat_t, class idx_t>
gpe_diskmat<gpe_mat_t, idx_t>::~gpe_diskmat () {
  close_tmp_blocks ();
  delete [] tmpfp;
}

template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::load_edgelist (const std::vector<std::string>& filenames, int ftype, int options) {
  read_and_split_list (filenames, ftype);
  csr_manager ();
}

/*! Reads the raw edgelist files
 *  Including some in GNU Zip format
 *
 *  It splits the edgelist into blocks in order to make
 *  the D/CSR building faster
 */
template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::read_and_split_list (const std::vector<std::string>& filenames, int ftype) {
  uint64_t** edlI = new uint64_t* [props.nblocks];
  uint64_t** edlO = new uint64_t* [props.nblocks];

  uint64_t* edlI_pos = new uint64_t [props.nblocks];

  if (props.sort_limit*2*props.nblocks > props.ram_limit) {
    std::cerr << "[GRAPHEE] To few memory allocated for sorting" << std::endl;
    std::cerr << "          with respect to the \'ram_limit\' setting" << std::endl;
    exit (-1);
  }

  const uint64_t maxElemsPerRamBlock = 
    props.sort_limit / sizeof(uint64_t);

  for (uint64_t i = 0; i < props.nblocks; i++) {
    edlI[i] = new uint64_t [maxElemsPerRamBlock];
    edlI_pos[i] = 0;
    edlO[i] = new uint64_t [maxElemsPerRamBlock];
  }

  std::vector<std::mutex> mutexes(props.nblocks);
  std::mutex read_mtx;

  size_t bytes;
  uint64_t blockID, fromID, toID;
  open_tmp_blocks (std::ios_base::out | std::ios_base::binary);

  std::stringstream* ssp_read;
  std::stringstream* ssp_load;
  std::stringstream* ssp_temp;
  std::stringstream  ss1;
  std::stringstream  ss2;

  ssp_load = &ss1;
  ssp_read = &ss2;

  load_GZ (filenames[0], ssp_load, &read_mtx);

  for (uint64_t i = 1; i < filenames.size(); i++) {
    read_mtx.lock();
    ssp_temp = ssp_read; ssp_read = ssp_load; ssp_load = ssp_temp;
    read_mtx.unlock();

    std::thread gzReadThd(load_GZ, filenames[i], ssp_load, &read_mtx);
    gzReadThd.detach();
    bytes = 0;

    while (!ssp_read->eof()) {
      //(*ssp_read) >> fromID >> toID;
      (*ssp_read) >> toID >> fromID;
      if (toID == fromID) continue; // oriented graph !
      blockID = fromID/props.window + toID/props.window*props.nslices;

      if (edlI_pos[blockID] < maxElemsPerRamBlock-1) {
        edlI[blockID][edlI_pos[blockID]  ] = fromID;
        edlI[blockID][edlI_pos[blockID]+1] = toID;
        edlI_pos[blockID] += 2;
        bytes += 2*sizeof(uint64_t);
      } else {
        mutexes[blockID].lock();
        std::memcpy(edlO[blockID], edlI[blockID],
            sizeof(uint64_t)*edlI_pos[blockID]);
        std::thread t(sort_and_save_list, edlO[blockID], edlI_pos[blockID],
            &(tmpfp[blockID]), &(mutexes[blockID]));
        t.detach();

        edlI_pos[blockID] = 0;
        edlI[blockID][edlI_pos[blockID]  ] = fromID;
        edlI[blockID][edlI_pos[blockID]+1] = toID;
        edlI_pos[blockID] += 2;
        bytes += 2*sizeof(uint64_t);
      }
    }
    std::cout << "[GRAPHEE] Sorting file " << filenames[i-1] << " with "
      << bytes/1024/1024 << " MB" << std::endl;
  }

  read_mtx.lock();
  ssp_temp = ssp_read; ssp_read = ssp_load; ssp_load = ssp_temp;
  read_mtx.unlock();

  while (!ssp_read->eof()) {
    //(*ssp_read) >> fromID >> toID;
    (*ssp_read) >> toID >> fromID;
    if (toID == fromID) continue; // oriented graph !
    blockID = fromID/props.window + toID/props.window*props.nslices;

    if (edlI_pos[blockID] < maxElemsPerRamBlock-1) {
      edlI[blockID][edlI_pos[blockID]  ] = fromID;
      edlI[blockID][edlI_pos[blockID]+1] = toID;
      edlI_pos[blockID] += 2;
      bytes += 2*sizeof(uint64_t);
    } else {
      mutexes[blockID].lock();
      std::memcpy(edlO[blockID], edlI[blockID],
          sizeof(uint64_t)*maxElemsPerRamBlock);
      std::thread t(sort_and_save_list, edlO[blockID], edlI_pos[blockID],
          &(tmpfp[blockID]), &(mutexes[blockID]));
      t.detach();

      edlI_pos[blockID] = 0;
      edlI[blockID][edlI_pos[blockID]  ] = fromID;
      edlI[blockID][edlI_pos[blockID]+1] = toID;
      edlI_pos[blockID] += 2;
      bytes += 2*sizeof(uint64_t);
    }
  }

  for (uint64_t i = 0; i < props.nblocks; i++) {
    if (edlI_pos[i] > 0) {
      mutexes[i].lock();
      std::memcpy(edlO[i], edlI[i], sizeof(uint64_t)*edlI_pos[i]);
      std::thread t(sort_and_save_list, edlO[i], edlI_pos[i],
          &(tmpfp[i]), &(mutexes[i]));
      t.detach();
      edlI_pos[i] = 0;
    }
  }

  for (uint64_t i = 0; i < props.nblocks; i++) {
    mutexes[i].lock();
    delete[] edlI[i];
    delete[] edlO[i];
    mutexes[i].unlock();
  }

  delete[] edlI;
  delete[] edlI_pos;
  delete[] edlO;

  close_tmp_blocks ();
}

template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::load_GZ (const std::string filename, std::stringstream* ss, std::mutex* read_mtx) {
  igzstream gz_ifp (filename.data());

  if (gz_ifp) {
    read_mtx->lock();
    ss->str(""); // empty string stream
    ss->clear(); // reset the position of EOS

    std::cout << "[GRAPHEE] [LOAD_GZ] Extracting file " << filename << "... ";
    (*ss) << gz_ifp.rdbuf();
    gz_ifp.close();

    std::cout << "extracted !" << std::endl;
    read_mtx->unlock();
  } else {
    std::cerr << "[GRAPHEE] [LOAD_GZ] Cannot open file " << filename << ", exiting...";
    exit (-1);
  }
}

template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::sort_and_save_list (uint64_t* block, uint64_t nelems,
    std::fstream* ofp, std::mutex* mtx) {

  if (nelems%2 != 0) {
    std::cerr << "[GRAPHEE] Wrong number of edges" << std::endl;
    exit (-1);
  }

  size_t nedges = nelems/2;
  uint64_t** edgeptrs = new uint64_t* [nedges];
  for (uint64_t i = 0; i < nedges; i++) {
    edgeptrs[i] = &(block[2*i]);
  }

  gpe_heap_sort(edgeptrs, nedges);

  ofp->write((const char*)block, sizeof(uint64_t)*nelems);

  mtx->unlock();

  delete[] edgeptrs;
}

template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::csr_manager () {
  std::vector<std::thread> csr_threads;

  open_tmp_blocks (std::ios_base::in | std::ios_base::binary);

  std::mutex mtx;
  std::condition_variable cond;

  size_t alloc_mem {0};
  size_t ram_limit {props.ram_limit};

  for (uint64_t line = 0; props.nslices; line++) {
    for (uint64_t col = 0; props.nslices; col++) {
      uint64_t bid = col + line*props.nslices;
      csr_threads.push_back(std::thread(csr_builder, props.window, line, col, &(tmpfp[bid]),
            &props, &alloc_mem, &mtx, &cond));
    }
  }

  for (uint64_t i = 0; i < props.nblocks; i++) {
    if (csr_threads[i].joinable()) {
      csr_threads[i].join();
    }
  }

  close_tmp_blocks ();
}

template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::csr_builder (uint64_t m, uint64_t line, uint64_t col, std::fstream* tmpfp,
    gpe_props* props, size_t* alloc_mem, std::mutex* mtx, std::condition_variable* cond) {
  tmpfp->seekg (0, tmpfp->end);
  size_t filelen = tmpfp->tellg ();
  tmpfp->seekg (0, tmpfp->beg);

  uint64_t nnz = filelen/(2*sizeof(uint64_t));
  uint64_t nsections = filelen/props->sort_limit + (filelen%props->sort_limit == 0 ? 0 : 1);

  size_t alloc_needs = (nnz + m+1)*sizeof(idx_t);

  if (alloc_needs > props->ram_limit) {
    mtx->lock();
    std::cerr << "[GRAPHEE] CSR block [" << line << ";" << col << "] needs " << alloc_needs/(1UL << 30) << "GB" << std::endl;
    std::cerr << "          which is more memory than \'ram_limit\' " << props->ram_limit/(1UL << 30) << "GB" << std::endl;
    mtx->unlock();
    return; // not exit(-1); because we are within a thread
  }

  std::unique_lock<std::mutex> mlock(*mtx);
  while ((*alloc_mem) + alloc_needs > props->ram_limit) {
    cond->wait(mlock);
  }
  (*alloc_mem) += alloc_needs;
  mlock.unlock();

  gpe_mat_t mat (*props, props->window, nnz);

  uint64_t* offsets = new uint64_t [nsections];
  uint64_t* ids = new uint64_t [nsections];

  for (uint64_t i = 0; i < nsections; i++) {
    offsets[i] = i*props->sort_limit;
    ids[i] = line*m;
  }

  uint64_t currentid = line*m;
  uint64_t minid;

  uint64_t edge[2];
  uint64_t offl = line*m;
  uint64_t offc = col*m;

  while (currentid < (line+1)*m) {
    minid = (line+1)*m;

    for (uint64_t sec = 0; sec < nsections; sec++) {
      if (ids[sec] == currentid && std::min(filelen, (sec+1)*props->sort_limit)) {
        if (tmpfp->tellg() != offsets[sec]) {
          tmpfp->seekg(offsets[sec], tmpfp->beg);
        }

        while (offsets[sec] < std::min(filelen, (sec+1)*props->sort_limit)) {
          tmpfp->read((char*)edge, 2*sizeof(uint64_t));

          /* Substract the offset */
          edge[0] -= offl;
          edge[1] -= offc;

          if (edge[0] == currentid) {
            if (edge[0] < std::numeric_limits<idx_t>::max() && edge[1] < std::numeric_limits<idx_t>::max()) {
              mat.sorted_fill ((idx_t) edge[0], (idx_t) edge[1]);
            } else {
              std::cerr << "[GRAPHEE] [CSR_BUILDER] Overflow of \'idx_t\', exiting..." << std::endl;
              return;
            }
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

  delete[] ids;
  delete[] offsets;

  if (!mat.verify()) {
    mtx->lock();
    std::cerr << "[GRAPHEE] Block [" << line << ";" << col << "] wrong conversion to CSR" << std::endl;
    mtx->unlock();
    return;
  }

  std::ostringstream filename;
  filename << props->name << "_csrblk_" << line << "_" << col << ".gpe";
  mat.save(filename.str(), gpe_mat_t::SNAPPY, offl, offc);

  mtx->lock();
  (*alloc_mem) -= alloc_needs;
  mtx->unlock();

  cond->notify_one();
}

template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::open_tmp_blocks (std::ios_base::openmode mode) {
  std::ostringstream oss;
  for (uint64_t i = 0; i < props.nblocks; i++) {
    oss.str("");
    oss << props.name << "_tmpblk" << i << ".gphee";
    tmpfp[i].open(oss.str(), mode);
    if (!tmpfp[i].is_open()) {
      std::cerr << "[GRAPHEE] Could not open file" << oss.str() << std::endl;
      exit(-1);
    }
  }
}

template <class gpe_mat_t, class idx_t>
void gpe_diskmat<gpe_mat_t, idx_t>::close_tmp_blocks () {
  for (uint64_t i = 0; i < props.nblocks; i++) {
    if (!tmpfp[i].is_open()) {
      std::cerr << "[GRAPHEE] Could not close file" << std::endl;
      exit(-1);
    }
    tmpfp[i].close();
  }
}

} // namespace graphee

#endif // GPE_DISKMAT_H
