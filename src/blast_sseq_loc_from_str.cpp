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

// based on ncbi_cxx--Aug_27_2007/src/objtools/readers/seqdb/seqdbvol.cpp by
// Kevin Bealer

#include "blast_sseq_loc_from_str.hpp"

#include <algorithm>

#include <corelib/ncbitime.hpp>

#include <util/sequtil/sequtil.hpp>
#include <util/sequtil/sequtil_convert.hpp>

#include <objmgr/util/sequence.hpp>

#include <serial/iterator.hpp>
#include <objects/seq/Bioseq.hpp>
#include <objects/seq/Seqdesc.hpp>
#include <objects/seq/Seq_descr.hpp>
#include <objects/seq/Seq_data.hpp>

#include <sstream>

// TODO: generalize to alphanumeric IDs (full NCBI or even non-NCBI)

using namespace ncbi;
using namespace ncbi::objects;
using namespace ncbi::blast;


static void load_protein_sequence(CSeq_inst& seqinst, const std::string& seq_data){
  std::vector<char> aa_data(seq_data.begin(), seq_data.end());
  seqinst.SetSeq_data().SetNcbistdaa().Set().swap(aa_data);
  seqinst.SetMol(CSeq_inst::eMol_aa);
}

static void load_nucleotide_sequence(CSeq_inst& seqinst, const std::string& seq_data){
  std::vector<char> packed;
  CSeqUtil::TCoding coding = CSeqUtil::e_not_set;
  CSeqConvert::Pack(seq_data, CSeqUtil::e_Iupacna, packed, coding, seq_data.size());
  if (coding == CSeqUtil::e_Ncbi2na){
    seqinst.SetSeq_data().SetNcbi2na().Set().swap(packed);
  } else if (coding == CSeqUtil::e_Ncbi2na){
    seqinst.SetSeq_data().SetNcbi4na().Set().swap(packed);  
  } else {
    throw std::runtime_error("Cannot pack seq_data.");
  }
  seqinst.SetMol(CSeq_inst::eMol_na);
}

CRef<CBioseq>
build_bioseq(const std::string& seq_data,
	     bool is_prot,             
	     int  seq_gi,
	     const std::string& title) {
  
  CRef<CBioseq> bioseq(new CBioseq);

  CSeq_inst& seqinst = bioseq->SetInst();
  if (is_prot) {
    load_protein_sequence(seqinst, seq_data);
  } else {
    load_nucleotide_sequence(seqinst, seq_data);
  }
  seqinst.SetLength(seq_data.size());
  seqinst.SetRepr(CSeq_inst::eRepr_raw);
  
  std::list< CRef<CSeq_id> > seqids;
  CRef<CSeq_id> seqid(new CSeq_id(CSeq_id::e_Gi, seq_gi));
  seqids.push_back(seqid);
  bioseq->SetId().swap(seqids);
  
  CRef<CSeqdesc> desc(new CSeqdesc);
  desc->SetTitle(title);
  bioseq->SetDescr().Set().push_back(desc);

  return bioseq;
}

SSeqLoc blast_sseq_loc_from_str::make(const std::string& seq_data, bool is_prot, int  seq_gi, const std::string& title,
				      ENa_strand strand, int from, int to) {
  
  if (seq_gi < 0) {
    throw std::runtime_error("seq_id cannot be negative.");
  }
  if (seq_data.empty()) {
    throw std::runtime_error("seq_data cannot be empty.");
  }
  if (title.empty()) {
    throw std::runtime_error("title cannot be empty.");
  }
  
  CRef<CScope> scope(new CScope(*_objmngr));
  scope->AddDefaults();

  from = std::max(from - 1, 0);
  to   = std::max(to - 1,   0);
  if (to <= 0 || to > seq_data.size()){
    to = seq_data.size() -1; // this is because is [from:to]  and NOT [from:to[
  }
  if (from > to) {
    from = to;
  }

  CRef<CSeq_entry> entry(new CSeq_entry);
  CRef<CBioseq> bioseq = build_bioseq(seq_data, is_prot, seq_gi, title);
  entry->SetSeq(*bioseq);
  scope->AddTopLevelSeqEntry(*entry);
  
  CRef<CSeq_loc> seqloc(new CSeq_loc());

  seqloc->SetInt().SetTo(to);
  seqloc->SetInt().SetFrom(from);
  seqloc->SetStrand(strand);
  seqloc->SetInt().SetId().Assign(*bioseq->GetId().front());
  
  SSeqLoc sl(seqloc, scope);
  return sl;
}

