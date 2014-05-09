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

#include "blast_sseq_loc_from_fasta.hpp"
#include "blast_sseq_loc_from_str.hpp"

#include <boost/python.hpp>

#include <iostream>
#include <assert.h>
#include <map>

using namespace boost::python;


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// blast_sseq_loc_from_string_wrapper
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct blast_sseq_loc_from_str_wrapper : public blast_sseq_loc_from_str, wrapper<blast_sseq_loc_from_str> {
  blast_sseq_loc_from_str_wrapper(PyObject* py_self) : blast_sseq_loc_from_str(), _py_self(py_self) {
  }

  ncbi::blast::SSeqLoc make(const std::string seq_data, bool is_prot, int seq_gi, const std::string title,
			    ncbi::objects::ENa_strand strand, int from, int to) {
    return blast_sseq_loc_from_str::make(seq_data, is_prot, seq_gi, title, strand, from, to);
  }
  
  PyObject* _py_self;
};


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// blast_sseq_loc_from_fasta_wrapper 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct blast_sseq_loc_from_fasta_wrapper : public blast_sseq_loc_from_fasta, wrapper<blast_sseq_loc_from_fasta> {

  blast_sseq_loc_from_fasta_wrapper(PyObject* py_self) : blast_sseq_loc_from_fasta(), _py_self(py_self) {
  }

  ncbi::blast::SSeqLoc make(str fasta, int s, 
			    ncbi::objects::ENa_strand strand, int from, int to, bool get_lcase_mask) {
    const char* fstr = extract<const char*>(fasta);
    return blast_sseq_loc_from_fasta::make(fstr, s, strand, from, to, get_lcase_mask);
  }

  std::string make_dummy(str fasta, ncbi::objects::ENa_strand strand, 
			 int from, int to, bool get_lcase_mask) {
    const char* fstr = extract<const char*>(fasta);

    std::string s(fstr);
    return s;
  }

  ncbi::blast::SSeqLoc make_zero_seq(str fasta, ncbi::objects::ENa_strand strand, 
				     int from, int to, bool get_lcase_mask) {
    return make(fasta, 0, strand, from, to, get_lcase_mask) ;
  }
  
  PyObject* _py_self;
};


static char make_from_fasta_doc[] = "make(FASTA, S=0, STRAND, FROM, TO, LCASE_MASK)\n\
 Returns a blast::SSeq_loc object from chars[FROM:TO[ of the S sequence contained in string FASTA\n\
 assign STRAND to it and do a lower case mask if LCASE_MASK is True. To select the whole sequence put FROM=0 and TO=0.";

static char make_from_str_doc[] = "make(SEQ_DATA, IS_PROT, GID, TITLE, STRAND, FROM, TO)\n\
 Returns a blast::SSeq_loc objectfrom chars[FROM:TO[ of the S sequence contained\
 in string SEQ_DATA\n if IS_PROT == False, assign STRAND to it. The underlying Bioseq will\
 have gid GID and title TITLE. To select the whole sequence put FROM=0\
 and TO=0.";

void export_blast_sseq_factories() 
{
  class_< blast_sseq_loc_from_fasta,
    boost::noncopyable, 
    blast_sseq_loc_from_fasta_wrapper
    >("blast_sseq_loc_from_fasta", 
      "A SSeq_Loc factory that generates ncbi::blast::SSeqLoc from fasta strings",
      init<>("Ready to generate!.")
      )
    .def("make", &blast_sseq_loc_from_fasta_wrapper::make, make_from_fasta_doc)
    .def("make", &blast_sseq_loc_from_fasta_wrapper::make_zero_seq)
    .def("make_dummy", &blast_sseq_loc_from_fasta_wrapper::make_dummy)
    ;

  class_< blast_sseq_loc_from_str,
    boost::noncopyable, 
    blast_sseq_loc_from_str_wrapper
    >("blast_sseq_loc_from_str", 
      "A SSeq_Loc factory that generates ncbi::blast::SSeqLoc from seq strings",
      init<>("Ready to generate!.")
      )
    .def("make", &blast_sseq_loc_from_str_wrapper::make, make_from_str_doc)
    ;
}


