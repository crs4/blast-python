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

// Based on ncbi_cxx--Aug_27_2007/src/algo/blast/api/demo/blast_input.cpp by
// Christiam Camacho

#include "blast_sseq_loc_from_fasta.hpp"

#include <corelib/ncbitime.hpp>
#include <objmgr/util/sequence.hpp>
#include <objtools/readers/fasta.hpp>

#include <serial/iterator.hpp>
#include <objects/seq/Bioseq.hpp>

#include <sstream>

using namespace ncbi;
using namespace ncbi::objects;
using namespace ncbi::blast;

/** Reads FASTA file and creates a TSeqLocVector type for the read sequence(s). 
 * Restricts sequences to an interval, if requested. 
 * @param in Input file stream [in]
 * @param objmgr Object manager reference [in]
 * @param strand What strand to use if it is a nucleotide sequence [in]
 * @param from Starting offset of an interval [in]
 * @param to Ending offset of an interval (end of sequence if 0) [in]
 * @param counter What index to start assigning local ids from? First unused 
 *                index on exit. [in] [out]
 * @param get_lcase_mask Should lower case be masked? [in]
 * @return Vector of sequence location structures.
 */
TSeqLocVector
BLASTGetSeqLocFromStream(std::istream& in, CObjectManager& objmgr, 
                         ENa_strand strand, int from, int to, 
                         int *counter, bool get_lcase_mask)
{
  TSeqLocVector retval;
  CRef<CSeq_entry> seq_entry;

  vector<CConstRef<CSeq_loc> > lcase_mask;

  CRef<CScope> scope(new CScope(objmgr));
  scope->AddDefaults();

  if (get_lcase_mask) {
    if ( !(seq_entry = ReadFasta(in, fReadFasta_AllSeqIds, counter, 
				 &lcase_mask)))
      throw std::runtime_error("Could not retrieve seq entry");
  } else {
    if ( !(seq_entry = ReadFasta(in, fReadFasta_AllSeqIds, counter)))
      throw std::runtime_error("Could not retrieve seq entry");
  }

  int index = 0;
  scope->AddTopLevelSeqEntry(*seq_entry);

  from = std::max(from - 1, 0);
  to = std::max(to - 1, 0);

  for (CTypeConstIterator<CBioseq> itr(ConstBegin(*seq_entry)); itr; ++itr) {

    CRef<CSeq_loc> seqloc(new CSeq_loc());
    TSeqPos seq_length = ncbi::objects::sequence::GetLength(*itr->GetId().front(), 
							    scope) - 1;

    if (to > 0 && to < seq_length)
      seqloc->SetInt().SetTo(to);
    else
      seqloc->SetInt().SetTo(seq_length);

    if (from > 0 && from < seq_length && from < to)
      seqloc->SetInt().SetFrom(from);
    else
      seqloc->SetInt().SetFrom(0);

    seqloc->SetInt().SetStrand(strand);
    seqloc->SetInt().SetId().Assign(*itr->GetId().front());

    //CRef<CScope> s(scope);
    SSeqLoc sl(seqloc, scope);

    if (get_lcase_mask) {
#if 0
      sl.mask.Reset(lcase_mask[index++]);
#else
      CSeq_loc* cs = const_cast<CSeq_loc*>(lcase_mask[index++].GetPointer());
      sl.mask.Reset(cs);
#endif
    }
    retval.push_back(sl);
  }

  return retval;
}


SSeqLoc blast_sseq_loc_from_fasta::make(const char* fasta_str, std::size_t s, ENa_strand strand, 
					std::size_t from, std::size_t to, bool get_lcase_mask) {
  std::istringstream is(fasta_str);
  TSeqLocVector vec = BLASTGetSeqLocFromStream(is, *_objmngr, strand, from, to, &_counter, get_lcase_mask);
  if (s < 0 ) {
    throw std::runtime_error("Only non negative sequence selections are allowed.");
  }
  if (vec.size() <= s ) {
    throw std::runtime_error("Could not read requested sequence from fasta file.");
  }
  return vec[s];
}
