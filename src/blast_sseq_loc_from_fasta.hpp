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
#ifndef _BLAST_SSEQ_LOC_FROM_FASTA_HPP_
#define _BLAST_SSEQ_LOC_FROM_FASTA_HPP_

#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/sseqloc.hpp>
#include <objmgr/object_manager.hpp>

class blast_sseq_loc_from_fasta {
 public:
  blast_sseq_loc_from_fasta() {
    _objmngr= ncbi::objects::CObjectManager::GetInstance(); 
    if (!_objmngr) {
      throw std::runtime_error("Could not initialize object manager");
    }
    _counter= 0;
  }

  blast_sseq_loc_from_fasta(const blast_sseq_loc_from_fasta& o)  {
    _objmngr = o._objmngr;
    _counter = o._counter;
  }
  
  ncbi::blast::SSeqLoc make(const char* fasta_str, ncbi::objects::ENa_strand strand, 
			    std::size_t from, std::size_t to, bool get_lcase_mask) {
    return make(fasta_str, 0, strand, from, to, get_lcase_mask);
  }

  ncbi::blast::SSeqLoc make(const char* fasta_str, std::size_t s, ncbi::objects::ENa_strand strand, 
			    std::size_t from, std::size_t to, bool get_lcase_mask);

 protected:
  ncbi::CRef<ncbi::objects::CObjectManager> _objmngr;
  int                                       _counter;
};

#endif // _BLAST_SSEQ_LOC_FROM_FASTA_HPP_
