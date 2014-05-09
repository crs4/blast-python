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

// Based on ncbi_cxx--Apr_22_2005/src/objmgr/util/sequence.cpp by
// Clifford Clausen

#include <ncbi_pch.hpp>
#include <serial/iterator.hpp>
#include <util/static_map.hpp>

#include <objmgr/object_manager.hpp>
#include <objmgr/scope.hpp>
#include <objmgr/seq_vector.hpp>
#include <objmgr/seq_vector_ci.hpp>
#include <objmgr/seqdesc_ci.hpp>
#include <objmgr/feat_ci.hpp>
#include <objmgr/bioseq_ci.hpp>
#include <objmgr/seq_entry_handle.hpp>
#include <objmgr/impl/handle_range_map.hpp>
#include <objmgr/impl/synonyms.hpp>

#include <objects/general/Int_fuzz.hpp>
#include <objects/general/Dbtag.hpp>
#include <objects/general/Object_id.hpp>
#include <objects/general/User_object.hpp>
#include <objects/general/User_field.hpp>

#include <objects/seq/Bioseq.hpp>
#include <objects/seq/Delta_ext.hpp>
#include <objects/seq/Delta_seq.hpp>
#include <objects/seq/MolInfo.hpp>
#include <objects/seq/Seg_ext.hpp>
#include <objects/seq/Seq_ext.hpp>
#include <objects/seq/Seq_inst.hpp>
#include <objects/seq/Seq_literal.hpp>

#include <objects/seqloc/Packed_seqpnt.hpp>
#include <objects/seqloc/Seq_bond.hpp>
#include <objects/seqloc/Seq_id.hpp>
#include <objects/seqloc/Seq_interval.hpp>
#include <objects/seqloc/Seq_loc.hpp>
#include <objects/seqloc/Seq_loc_equiv.hpp>
#include <objects/seqloc/Seq_loc_mix.hpp>
#include <objects/seqloc/Seq_point.hpp>

#include <objects/seqset/Seq_entry.hpp>

#include <objects/seqfeat/Org_ref.hpp>
#include <objects/seqfeat/BioSource.hpp>
#include <objects/seqfeat/Cdregion.hpp>
#include <objects/seqfeat/Code_break.hpp>
#include <objects/seqfeat/Genetic_code.hpp>
#include <objects/seqfeat/Genetic_code_table.hpp>
#include <objects/seqfeat/Seq_feat.hpp>

#include <objmgr/seq_loc_mapper.hpp>
#include <objmgr/seq_entry_ci.hpp>
#include <objmgr/util/sequence.hpp>
#include <util/strsearch.hpp>

#include <list>
#include <algorithm>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(objects)
BEGIN_SCOPE(sequence)

struct SGap {
    SGap(TSeqPos start, TSeqPos length) : m_Start(start), m_Length(length) { }
    TSeqPos GetEnd(void) const { return m_Start + m_Length - 1; }

    TSeqPos m_Start, m_Length;
};
typedef list<SGap> TGaps;


static bool s_IsGap(const CSeq_loc& loc, CScope& scope)
{
    if (loc.IsNull()) {
        return true;
    }

    CTypeConstIterator<CSeq_id> id(loc);
    CBioseq_Handle handle = scope.GetBioseqHandle(*id);
    if (handle  &&  handle.GetInst_Repr() == CSeq_inst::eRepr_virtual) {
        return true;
    }

    return false; // default
}

static TGaps s_FindGaps(const CSeq_ext& ext, CScope& scope)
{
    TSeqPos pos = 0;
    TGaps   gaps;

    switch (ext.Which()) {
    case CSeq_ext::e_Seg:
        ITERATE (CSeg_ext::Tdata, it, ext.GetSeg().Get()) {
            TSeqPos length = sequence::GetLength(**it, &scope);
            if (s_IsGap(**it, scope)) {
                gaps.push_back(SGap(pos, length));
            }
            pos += length;
        }
        break;

    case CSeq_ext::e_Delta:
        ITERATE (CDelta_ext::Tdata, it, ext.GetDelta().Get()) {
            switch ((*it)->Which()) {
            case CDelta_seq::e_Loc:
            {
                const CSeq_loc& loc = (*it)->GetLoc();
                TSeqPos length = sequence::GetLength(loc, &scope);
                if (s_IsGap(loc, scope)) {
                    gaps.push_back(SGap(pos, length));
                }
                pos += length;
                break;
            }

            case CDelta_seq::e_Literal:
            {
                const CSeq_literal& lit    = (*it)->GetLiteral();
                TSeqPos             length = lit.GetLength();
                if ( !lit.IsSetSeq_data() ) {
                    gaps.push_back(SGap(pos, length));
                }
                pos += length;
                break;
            }

            default:
                ERR_POST(Warning << "CFastaOstream::WriteSequence: "
                         "unsupported Delta-seq selection "
                         << CDelta_seq::SelectionName((*it)->Which()));
                break;
            }
        }

    default:
        break;
    }

    return gaps;
}


static TGaps s_AdjustGaps(const TGaps& gaps, const CSeq_loc& location)
{
    // assume location matches handle
    const TSeqPos         kMaxPos = numeric_limits<TSeqPos>::max();
    TSeqPos               pos     = 0;
    TGaps::const_iterator gap_it  = gaps.begin();
    TGaps                 adjusted_gaps;
    SGap                  new_gap(kMaxPos, 0);

    for (CSeq_loc_CI loc_it(location);  loc_it  &&  gap_it != gaps.end();
         pos += loc_it.GetRange().GetLength(), ++loc_it) {
        CSeq_loc_CI::TRange range = loc_it.GetRange();

        if (new_gap.m_Start != kMaxPos) {
            // in progress
            if (gap_it->GetEnd() < range.GetFrom()) {
                adjusted_gaps.push_back(new_gap);
                new_gap.m_Start = kMaxPos;
                ++gap_it;
            } else if (gap_it->GetEnd() <= range.GetTo()) {
                new_gap.m_Length += gap_it->GetEnd() - range.GetFrom() + 1;
                adjusted_gaps.push_back(new_gap);
                new_gap.m_Start = kMaxPos;
                ++gap_it;
            } else {
                new_gap.m_Length += range.GetLength();
                continue;
            }
        }

        while (gap_it != gaps.end()  &&  gap_it->GetEnd() < range.GetFrom()) {
            ++gap_it; // skip
        }
        // we may be out of gaps now...
        if (gap_it == gaps.end()) {
            break;
        }

        if (gap_it->m_Start <= range.GetFrom()) {
            if (gap_it->GetEnd() <= range.GetTo()) {
                adjusted_gaps.push_back
                    (SGap(pos, gap_it->GetEnd() - range.GetFrom() + 1));
                ++gap_it;
            } else {
                new_gap.m_Start  = pos;
                new_gap.m_Length = range.GetLength();
                continue;
            }
        }

        while (gap_it->m_Start <= range.GetTo()) {
            TSeqPos pos2 = pos + gap_it->m_Start - range.GetFrom();
            if (gap_it->GetEnd() <= range.GetTo()) {
                adjusted_gaps.push_back(SGap(pos2, gap_it->m_Length));
                ++gap_it;
            } else {
                new_gap.m_Start  = pos2;
                new_gap.m_Length = range.GetTo() - gap_it->m_Start + 1;
            }
        }
    }

    if (new_gap.m_Start != kMaxPos) {
        adjusted_gaps.push_back(new_gap);
    }

    return adjusted_gaps;
}

/*
enum EFlags {
  eAssembleParts   = 0x1,
  eInstantiateGaps = 0x2
};
*/
// 
/** Returns a Iupac coded string containing the subsequence [beg, end[.
    FIXME: it is probable that it will not work with AA sequences (not tested!)
 */
std::string extract_sequence_from_bioseq(const CBioseq_Handle& handle,
					 bool instantiate_gaps,
					 int beg, int end)
{
  TSeqPos start = static_cast<TSeqPos>(beg);
  TSeqPos stop  = static_cast<TSeqPos>(end);


  std::string result;
  CConstRef<CBioseq> seq  = handle.GetCompleteBioseq();
  const CSeq_inst& inst = seq->GetInst();
  if (!inst.IsSetSeq_data()) {
    return result;
  }
  CSeqVector v;
  v = handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac);
  bool is_na = inst.GetMol() != CSeq_inst::eMol_aa;
  // autodetection is sometimes broken (!)
  v.SetCoding(is_na ? CSeq_data::e_Iupacna : CSeq_data::e_Iupacaa);

  TSeqPos              pos  = start;
  CSeqVector::TResidue gap  = v.GetGapChar();
  TGaps                gaps;
  CScope&              scope   = handle.GetScope();
  const TSeqPos        kMaxPos = numeric_limits<TSeqPos>::max();


  if ( !inst.IsSetSeq_data()  &&  inst.IsSetExt() ) {
    gaps = s_FindGaps(inst.GetExt(), scope);
  }
  gaps.push_back(SGap(kMaxPos, 0));

  stop = (stop == 0)? v.size() : stop;

  if (start >= stop) {
    return result;
  }

  std::size_t buffer_size = stop - start;

  TGaps filtered_gaps;
  for(TGaps::iterator it = gaps.begin(); it != gaps.end(); ++it){
    if (it->GetEnd() >= pos) {
      TSeqPos gap_start = std::max(pos, it->m_Start);
      filtered_gaps.push_back(SGap(gap_start, it->m_Length - (gap_start - it->m_Start)));
    }
  }
  gaps = filtered_gaps;
#if 0
  //FIXME: is this really needed?
  for(TGaps::iterator it = gaps.begin(); it != gaps.end(); ++it){
    buffer_size += it->m_Length;
  }
#endif
  result.resize(buffer_size);
  while (pos < stop) {
    unsigned int limit = std::min(stop, gaps.front().m_Start) - pos;
    v.GetSeqData(pos, pos + limit, result);
    pos += limit;
    if (pos == gaps.front().m_Start) {
      if (instantiate_gaps){
	result += std::string(gaps.front().m_Length, gap);
      }
      pos += gaps.front().m_Length;
      gaps.pop_front();
    }
  }
  return result;
}

END_SCOPE(sequence)
END_SCOPE(objects)
END_NCBI_SCOPE
