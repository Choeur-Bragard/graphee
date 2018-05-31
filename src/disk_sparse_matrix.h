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

#include <gzstream.h>
#include <snappy.h>

#include "graphee.h"

#define SORT_NAME gpe
#define SORT_TYPE uint64_t*
#define SORT_DIM 2
#define SORT_CMP(x, y)  (*(x) < *(y) ? -1 : (*(x) == *(y) ? 0 : 1))
#define SORT_SWAP(x, y) {int i; for (i=0; i<SORT_DIM; i++) {uint64_t __SORT_SWAP = *(x+i); *(x+i) = *(y+i); *(y+i) = __SORT_SWAP;}}
#include "sort/sort.h"

namespace graphee {

template <typename gpe_mat_t>
class gpe_diskmat {
public:
  gpe_diskmat ();
  gpe_diskmat (gpe_props props, std::string matrix_name);

  ~gpe_diskmat ();

  void init_diskmat (gpe_props props, std::string matrix_name);

  void load_edgelist (const std::vector<std::string>& filenames, int ftype, int options);

  void get_matrix_block (uint64_t line, uint64_t col, gpe_mat_t& mat);

  enum {PLAIN, GZ};
  static const int DIRECT = 0x00000001; 
  static const int TRANS  = 0x00000010; 
  static const int UO     = 0x00000100; 

  typedef gpe_mat_t matrix_type;

  bool empty ();

private:
  std::ostringstream log;
  std::ostringstream wrn;
  std::ostringstream err;

  gpe_props props;

  std::fstream *tmpfp;

  std::string matrix_name;

  void read_and_split_list (const std::vector<std::string>& filenames, int ftype);

  static void sort_and_save_list (uint64_t* block, uint64_t nelems,
      std::fstream* fp, std::mutex* mtx);

  static void load_GZ (const std::string filename, std::stringstream* ss,
      std::mutex* read_mtx);

  void diskblock_manager ();
  static void diskblock_builder (gpe_diskmat<gpe_mat_t>* dmat, uint64_t line, uint64_t col,
      std::fstream* tmpfp, size_t* alloc_mem, std::mutex* mtx, std::condition_variable* cond);

  std::string get_block_filename (uint64_t line, uint64_t col);

  void open_tmp_blocks (std::ios_base::openmode mode);
  void close_tmp_blocks ();
};

template <typename gpe_mat_t>
gpe_diskmat<gpe_mat_t>::gpe_diskmat () {
}

template <typename gpe_mat_t>
gpe_diskmat<gpe_mat_t>::gpe_diskmat (gpe_props arg_props, std::string arg_matrix_name) {
  props = arg_props;
  matrix_name = arg_matrix_name;
  tmpfp = new std::fstream [props.nblocks];
}

template <typename gpe_mat_t>
gpe_diskmat<gpe_mat_t>::~gpe_diskmat () {
  close_tmp_blocks ();
  delete [] tmpfp;
}

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::init_diskmat (gpe_props arg_props, std::string arg_matrix_name) {
  props = arg_props;
  matrix_name = arg_matrix_name;
  tmpfp = new std::fstream [props.nblocks];
}

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::load_edgelist (const std::vector<std::string>& filenames, int ftype, int options) {
  read_and_split_list (filenames, ftype);
  diskblock_manager ();
}

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::get_matrix_block (uint64_t line, uint64_t col, gpe_mat_t& mat) {
  log.str("");
  log << "Start to load disk block [" << line << ":" << col << "]";
  gpe_log (log.str());

  mat.load(get_block_filename(line, col));

  log.str("");
  log << "Loaded disk block [" << line << ":" << col << "]";
  gpe_log (log.str());
}

template <typename gpe_mat_t>
bool gpe_diskmat<gpe_mat_t>::empty () {
  return (props.nvertices == 0);
}

/*! Reads the raw edgelist files
 *  Including some in GNU Zip format
 *
 *  It splits the edgelist into blocks in order to make
 *  the D/CSR building faster
 */
template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::read_and_split_list (const std::vector<std::string>& filenames, int ftype) {
  uint64_t** edlI = new uint64_t* [props.nblocks];
  uint64_t** edlO = new uint64_t* [props.nblocks];

  uint64_t* edlI_pos = new uint64_t [props.nblocks];

  if (props.sort_limit*2*props.nblocks > props.ram_limit) {
    err.str("");
    err << "To few memory allocated for sorting with respect to the \'ram_limit\' setting";
    gpe_error (err.str());
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

    log.str("");
    log << "Sorting file " << filenames[i-1] << " with "
      << bytes/1024/1024 << " MB";
    gpe_log (log.str());
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

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::load_GZ (const std::string filename, std::stringstream* ss, std::mutex* read_mtx) {
  std::ostringstream log;
  std::ostringstream err;

  igzstream gz_ifp (filename.data());

  if (gz_ifp) {
    read_mtx->lock();
    ss->str(""); // empty string stream
    ss->clear(); // reset the position of EOS

    log.str("");
    log << "Extracting file " << filename;
    gpe_log (log.str());
    (*ss) << gz_ifp.rdbuf();
    gz_ifp.close();

    read_mtx->unlock();
  } else {
    err.str("");
    err << "Cannot open file " << filename << ", exiting...";
    gpe_error (err.str());
    exit (-1);
  }
}

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::sort_and_save_list (uint64_t* block, uint64_t nelems,
    std::fstream* ofp, std::mutex* mtx) {

  if (nelems%2 != 0) {
    gpe_error ("Wrong number of edges");
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

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::diskblock_manager () {
  std::vector<std::thread> diskblock_threads;

  open_tmp_blocks (std::ios_base::in | std::ios_base::binary);

  std::mutex mtx;
  std::condition_variable cond;

  size_t alloc_mem {0};
  size_t ram_limit {props.ram_limit};

  for (uint64_t line = 0; line < props.nslices; line++) {
    for (uint64_t col = 0; col < props.nslices; col++) {
      uint64_t bid = line + col*props.nslices;
      diskblock_threads.push_back(std::thread(diskblock_builder, this, line, col, (&tmpfp[bid]), &alloc_mem, &mtx, &cond));
    }
  }

  for (uint64_t i = 0; i < props.nblocks; i++) {
    if (diskblock_threads[i].joinable()) {
      diskblock_threads[i].join();
    }
  }

  close_tmp_blocks ();
}

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::diskblock_builder (gpe_diskmat<gpe_mat_t>* dmat, uint64_t line, uint64_t col,
      std::fstream* tmpfp, size_t* alloc_mem, std::mutex* mtx, std::condition_variable* cond) {
  std::ostringstream log;
  std::ostringstream err;

  gpe_props props = dmat->props;

  tmpfp->seekg (0, tmpfp->end);
  size_t filelen = tmpfp->tellg ();
  tmpfp->seekg (0, tmpfp->beg);

  uint64_t nnz = filelen/(2*sizeof(uint64_t));
  uint64_t nsections = filelen/props.sort_limit + (filelen%props.sort_limit == 0 ? 0 : 1);

  size_t alloc_needs = (props.window+1)*sizeof(uint64_t) + (nnz)*sizeof(typename gpe_mat_t::index_type);

  if (alloc_needs > props.ram_limit) {
    mtx->lock();
    err.str("");
    err << "Disk block [" << line << ";" << col << "] needs " << alloc_needs/(1UL << 30) << "GB";
    gpe_error (err.str());
    err.str("");
    err << "which is more memory than \'ram_limit\' " << props.ram_limit/(1UL << 30) << "GB";
    gpe_error (err.str());
    mtx->unlock();
    return; // not exit(-1); because we are within a thread
  } else {
    mtx->lock();
    log.str("");
    log << "Starting disk block conversion [" << line << ";" << col << "]";
    gpe_log (log.str());
    mtx->unlock();
  }

  std::unique_lock<std::mutex> mlock(*mtx);
  while ((*alloc_mem) + alloc_needs > props.ram_limit) {
    cond->wait(mlock);
  }
  (*alloc_mem) += alloc_needs;
  mlock.unlock();

  gpe_mat_t mat (props, props.window, nnz);

  uint64_t* offsets = new uint64_t [nsections];
  uint64_t* ids = new uint64_t [nsections];

  for (uint64_t i = 0; i < nsections; i++) {
    offsets[i] = i*props.sort_limit;
    ids[i] = line*props.window;
  }

  uint64_t currentid = line*props.window;
  uint64_t minid;

  uint64_t edge[2];
  uint64_t offl = line*props.window;
  uint64_t offc = col*props.window;

  while (currentid < (line+1)*props.window) {
    minid = (line+1)*props.window;

    for (uint64_t sec = 0; sec < nsections; sec++) {
      if (ids[sec] == currentid && std::min(filelen, (sec+1)*props.sort_limit)) {
        if (tmpfp->tellg() != offsets[sec]) {
          tmpfp->seekg(offsets[sec], tmpfp->beg);
        }

        while (offsets[sec] < std::min(filelen, (sec+1)*props.sort_limit)) {
          tmpfp->read((char*)edge, 2*sizeof(uint64_t));

          if (edge[0] == currentid) {
            /* Substract the offset */
            edge[0] -= offl;
            edge[1] -= offc;
            if (edge[0] < std::numeric_limits<typename gpe_mat_t::index_type>::max() 
                && edge[1] < std::numeric_limits<typename gpe_mat_t::index_type>::max()) {

              mat.sorted_fill (static_cast<typename gpe_mat_t::index_type>(edge[0]), 
                  static_cast<typename gpe_mat_t::index_type>(edge[1]));
            } else {
              gpe_error ("Overflow of \'gpe_mat_t::index_type\', exiting...");
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
    err.str("");
    err << "Block [" << line << ";" << col << "] conversion to \'" << mat.matrix_type << "\' failed !";
    gpe_error (err.str());
    mtx->unlock();
    return;
  } else {
    mtx->lock();
    log.str("");
    log << "Block [" << line << ";" << col << "] conversion to \'" << mat.matrix_type << "\' succeed !";
    gpe_log (log.str());
    mtx->unlock();
  }

  mat.save(dmat->get_block_filename(line, col), gpe_mat_t::SNAPPY, offl, offc);

  mtx->lock();
  (*alloc_mem) -= alloc_needs;
  mtx->unlock();

  cond->notify_one();
}

template <typename gpe_mat_t>
std::string gpe_diskmat<gpe_mat_t>::get_block_filename (uint64_t line, uint64_t col) {
  std::ostringstream matrixname;
  matrixname << props.name << "_" << matrix_name << "_dmatblk_" << line << "_" << col << ".gpe";
  return matrixname.str();
}

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::open_tmp_blocks (std::ios_base::openmode mode) {
  std::ostringstream blockname;

  for (uint64_t line = 0; line < props.nslices; line++) {
    for (uint64_t col = 0; col < props.nslices; col++) {
      uint64_t bid = line + col*props.nslices;
      blockname.str("");
      blockname << props.name << "_" << matrix_name << "_tmpblk_" << line << "_" << col << ".gpe";
      tmpfp[bid].open(blockname.str(), mode);

      if (!tmpfp[bid].is_open()) {
        err.str("");
        err << "Could not open file: " << blockname.str();
        gpe_error (err.str());
        exit(-1);
      }
    }
  }
}

template <typename gpe_mat_t>
void gpe_diskmat<gpe_mat_t>::close_tmp_blocks () {
  for (uint64_t bid = 0; bid < props.nblocks; bid++) {
    if (tmpfp[bid].is_open()) {
      tmpfp[bid].close();
    }
  }
}

} // namespace graphee

#endif // GPE_DISKMAT_H
