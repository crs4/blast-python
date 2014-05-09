// BEGIN_COPYRIGHT
// 
// Copyright (C) 2014 CRS4.
// 
// This file is part of blast-python.
// 
// blast-python is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
// 
// blast-python is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
// 
// You should have received a copy of the GNU General Public License along
// with blast-python.  If not, see <http://www.gnu.org/licenses/>.
// 
// END_COPYRIGHT

#include <algo/blast/api/bl2seq.hpp>
#include <algo/blast/core/blast_options.h> // for the default values defines
#include <corelib/ncbitime.hpp>
#include <objmgr/util/sequence.hpp>
#include <objmgr/util/seq_loc_util.hpp>
#include <objtools/readers/fasta.hpp>

#include <serial/iterator.hpp>
#include <objects/seq/Bioseq.hpp>

#include <boost/python.hpp>
//#include <boost/numeric_cast.hpp>

#include <iostream>
#include <assert.h>
#include <map>


using namespace boost::python;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Helpers
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "cseq_sequence_extractor.hpp"

static inline int convert_index(int size, int idx) {
  idx = ((idx<0) ? size : 0) + idx; 
  return std::min(size, std::max(0, idx));
}

// from boost.python

static inline void extract_slice_data(PySliceObject* slice, std::size_t size, int& beg_, int& end_) {
  if (Py_None != slice->step) { 
    PyErr_SetString( PyExc_IndexError, "slice step size not supported.");
    throw_error_already_set();
  }
  std::size_t min_index = 0;
  std::size_t max_index = size;
  
  if (Py_None == slice->start) {
    beg_ = min_index;
  }  else {
    long from = extract<long>( slice->start);
    if (from < 0) // Negative slice index
      from += max_index;
    if (from < 0) // Clip lower bounds to zero
      from = 0;
    //    beg_ = boost::numeric_cast<std::size_t>(from);
    beg_ = static_cast<int>(from);
    if (beg_ > max_index) // Clip upper bounds to max_index.
      beg_ = max_index;
  }
  if (Py_None == slice->stop) {
    end_ = max_index;
  } else {
    long to = extract<long>( slice->stop);
    if (to < 0)
      to += max_index;
    if (to < 0)
      to = 0;
    //end_ = boost::numeric_cast<std::size_t>(to);
    end_ = static_cast<int>(to);
    if (end_ > max_index)
      end_ = max_index;
  }
}        


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ncbi_blast_SSeqLoc_wrapper
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//namespace { // put here a namespace if needed. This is only to avoid indenting if we change idea...

  struct ncbi_blast_SSeqLoc_wrapper : public ncbi::blast::SSeqLoc {

    ncbi_blast_SSeqLoc_wrapper(PyObject* py_self) :
      ncbi::blast::SSeqLoc(), _py_self(py_self) { }

    ncbi_blast_SSeqLoc_wrapper(PyObject* py_self, const ncbi::blast::SSeqLoc& ssl) :
      ncbi::blast::SSeqLoc(ssl), _py_self(py_self) { }

    int get_len() const {
      return ncbi::objects::sequence::GetLength(*seqloc, scope);
    }

    ncbi::objects::ENa_strand get_strand() const {
      return ncbi::objects::sequence::GetStrand(*seqloc, scope);
    }

    std::string get_id() const {
      ncbi::objects::CBioseq_Handle  bsh = ncbi::objects::sequence::GetBioseqFromSeqLoc(*seqloc, *scope);
      return  ncbi::objects::CSeq_id::GetStringDescr(*bsh.GetBioseqCore(), ncbi::objects::CSeq_id::eFormat_FastA);
    }

    std::string get_title() const {
      ncbi::objects::CBioseq_Handle  bsh = ncbi::objects::sequence::GetBioseqFromSeqLoc(*seqloc, *scope);
      return ncbi::objects::sequence::GetTitle(bsh, 0);
    }

    std::string get_sequence() const {
      ncbi::objects::CBioseq_Handle  bsh = ncbi::objects::sequence::GetBioseqFromSeqLoc(*seqloc, *scope);
      return ncbi::objects::sequence::extract_sequence_from_bioseq(bsh, true, 0, 0);
    }

    std::string get_item(PyObject* i) const {
      ncbi::objects::CBioseq_Handle  bsh = ncbi::objects::sequence::GetBioseqFromSeqLoc(*seqloc, *scope);
      int beg, end;
      if (PySlice_Check(i)) {
	PySliceObject* slice = reinterpret_cast<PySliceObject*>(i);
	extract_slice_data(slice, get_len(), beg, end);
      } else {
	int idx = convert_index(get_len(), extract<int>(i));
	beg = idx;
	end = idx+1;
      }
      if (beg == end) {
        return std::string("");
      } else {
        return ncbi::objects::sequence::extract_sequence_from_bioseq(bsh, 
                                                                     true, 
                                                                     beg, end);
      }
    }
    
    bool is_masked() const  {
      return ! mask.Empty();
    }

    PyObject* _py_self;
  };

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// export_blast_sseq
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void export_blast_sseq() 
{
  class_< ncbi::blast::SSeqLoc, ncbi_blast_SSeqLoc_wrapper>("SSeqLoc", 
      "C++ glue to ncbi::blast::SSeqLoc",
      init<>("Initializes as an empty sequence. FIXME docs...")
      )
    .add_property("title",  &ncbi_blast_SSeqLoc_wrapper::get_title)
    .add_property("id",     &ncbi_blast_SSeqLoc_wrapper::get_id)
    .add_property("strand", &ncbi_blast_SSeqLoc_wrapper::get_strand)
    .add_property("length", &ncbi_blast_SSeqLoc_wrapper::get_len)
    .add_property("is_masked", &ncbi_blast_SSeqLoc_wrapper::is_masked)
    .def("__len__", &ncbi_blast_SSeqLoc_wrapper::get_len)
    .def("__getitem__", &ncbi_blast_SSeqLoc_wrapper::get_item)
    .def("get_sequence", &ncbi_blast_SSeqLoc_wrapper::get_sequence)
    ;
}


//}
