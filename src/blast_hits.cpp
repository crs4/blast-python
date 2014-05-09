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
#include <algo/blast/core/blast_hits.h>
#include <algo/blast/core/gapinfo.h>

#include <iostream>
#include <boost/python.hpp>

#include <assert.h>

using namespace boost::python;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Helpers
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BlastHitList
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


struct ncbi_blast_BlastHitList_wrapper : public BlastHitList {

  ncbi_blast_BlastHitList_wrapper(PyObject* py_self) : _py_self(py_self) {}

  ncbi_blast_BlastHitList_wrapper(PyObject* py_self, 
				  const BlastHitList& hl) :
    BlastHitList(hl), _py_self(py_self) {}
  
  static int    get_len(const BlastHitList* hl) { return hl->hsplist_count; }
  static int    get_allocated_size(const BlastHitList* hl) { return hl->hsplist_max; }
  static double get_worst_evalue(const BlastHitList* hl) { return hl->worst_evalue; }
  static int    get_low_score(const BlastHitList* hl) { return hl->low_score; }
  
  static const BlastHSPList* get_item(const BlastHitList* hl, int i) { 
    if (i < 0) { throw std::runtime_error("Negative index"); }
    if (i >= hl->hsplist_count) {
      boost::python::objects::stop_iteration_error();
    }
    return hl->hsplist_array[i];
  }
  PyObject* _py_self;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BlastHSPList
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct ncbi_blast_BlastHSPList_wrapper : BlastHSPList {

  ncbi_blast_BlastHSPList_wrapper(PyObject* py_self) : _py_self(py_self) {}

  ncbi_blast_BlastHSPList_wrapper(PyObject* py_self, 
				  const BlastHSPList& hl) :
    BlastHSPList(hl), _py_self(py_self) {}
  
  static int    get_len(const BlastHSPList* hl) { return hl->hspcnt;}
  static int    get_ordinal_id_subject_sequence(const BlastHSPList* hl) { return hl->oid;}
  static int    get_query_index(const BlastHSPList* hl) { return hl->query_index;}
  static double get_best_evalue(const BlastHSPList* hl) { return hl->best_evalue;}
  
  static const BlastHSP* get_item(const BlastHSPList* hl, int i) { 
    if (i < 0) {
      throw std::runtime_error("Negative index");
    }
    if (i == hl->hspcnt) {
      boost::python::objects::stop_iteration_error();
    }
    return hl->hsp_array[i];
  }
  PyObject* _py_self;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BlastHSP
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct ncbi_blast_BlastHSP_wrapper : BlastHSP {

  ncbi_blast_BlastHSP_wrapper(PyObject* py_self) : _py_self(py_self) {}

  ncbi_blast_BlastHSP_wrapper(PyObject* py_self, 
			      const BlastHSP& hl) :
    BlastHSP(hl), _py_self(py_self) {}
  
  static int      get_score(const BlastHSP* hl)          { return hl->score;}
  static int      get_num_ident(const BlastHSP* hl)      { return hl->num_ident;}
  static double   get_bit_score(const BlastHSP* hl)      { return hl->bit_score;}
  static double   get_evalue(const BlastHSP* hl)         { return hl->evalue;}

  static tuple    map_to_tuple(const BlastSeg& bs) {
    return make_tuple(bs.frame, bs.offset, bs.end, bs.gapped_start);
  }

  static tuple    get_query(const BlastHSP* hl)          { return map_to_tuple(hl->query);}
  static tuple    get_subject(const BlastHSP* hl)        { return map_to_tuple(hl->subject);}


  static int      get_context(const BlastHSP* hl)        { return hl->context;}
  static int      get_num_linked_hsp(const BlastHSP* hl) { return hl->num;}
  //  static Uint4    get_pattern_length(const BlastHSP* hl) { return hl->pattern_length;}

  static GapEditScript* get_gap_info(const BlastHSP* hl) { return hl->gap_info;}

  PyObject* _py_self;
};

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// GapEditScript
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

typedef struct gap_edit_script_iterator_s {
  gap_edit_script_iterator_s(const GapEditScript* p) : _p(p), _index(0) {}

  struct gap_edit_script_iterator_s get_iter() const { return *(this); }
  
  tuple get_next() {
    
    if (_index == _p->size) {
      boost::python::objects::stop_iteration_error();
    }
    tuple t = make_tuple(_p->op_type[_index], _p->num[_index]);
    _index++;
    return t;
  }
  const GapEditScript* _p;
  Uint4                _index;
} gap_edit_script_iterator;


struct ncbi_blast_GapEditScript_wrapper : GapEditScript {

  ncbi_blast_GapEditScript_wrapper(PyObject* py_self) : 
    _py_self(py_self) {}
  
  ncbi_blast_GapEditScript_wrapper(PyObject* py_self, 
                                   const GapEditScript& hl) :
    GapEditScript(hl), _py_self(py_self) {}

  static int get_len(const GapEditScript* g) { return g->size;}
  static gap_edit_script_iterator get_iterator(const GapEditScript* g) {
    return gap_edit_script_iterator(g);
  }

  PyObject* _py_self;
};



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Exporting class definitions.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void export_blast_hits() 
{

  enum_<EGapAlignOpType> ("EGapAlignOpType")
    .value("eGapAlignDel", eGapAlignDel) // "Deletion: a gap in query ")
    .value("eGapAlignDel2", eGapAlignDel2) // "Frame shift deletion of two nucleotides ")
    .value("eGapAlignDel1", eGapAlignDel1) // "Frame shift deletion of one nucleotide ")
    .value("eGapAlignSub", eGapAlignSub) // "Substitution ")
    .value("eGapAlignIns1", eGapAlignIns1) // "Frame shift insertion of one nucleotide ")
    .value("eGapAlignIns2", eGapAlignIns2) // "Frame shift insertion of two nucleotides ")
    .value("eGapAlignIns", eGapAlignIns) // "Insertion: a gap in subject ")
    .value("eGapAlignDecline", eGapAlignDecline) // "Non-aligned region ")
    .value("eGapAlignInvalid", eGapAlignInvalid) // "Invalid operation ")
    ;

  class_< BlastHitList,
    ncbi_blast_BlastHitList_wrapper>("BlastHitList",
				     "C++ glue to access Blast raw results",
				     init<>("Initializes with junk, most likely..")
				     )
      .def("__len__",      &ncbi_blast_BlastHitList_wrapper::get_len)
      .def("__getitem__",  &ncbi_blast_BlastHitList_wrapper::get_item,  return_internal_reference<>())
      .add_property("allocated_size", &ncbi_blast_BlastHitList_wrapper::get_allocated_size, 
		    "Maximal allowed size of the HSP lists array")
      .add_property("worst_evalue", &ncbi_blast_BlastHitList_wrapper::get_worst_evalue, 
		    "Highest of the best e-values among the HSP lists")
      .add_property("low_score", &ncbi_blast_BlastHitList_wrapper::get_low_score, 
		    "Highest of the best e-values among the HSP lists")
      ;

  class_< BlastHSPList,
    ncbi_blast_BlastHSPList_wrapper>("BlastHSPList",
				     "C++ glue to access Blast HSP lists for a given sequence after the gapped alignment",
				     init<>("Initializes with junk, most likely..")
				     )
      .def("__len__",      &ncbi_blast_BlastHSPList_wrapper::get_len)
      .def("__getitem__",  &ncbi_blast_BlastHSPList_wrapper::get_item,  return_internal_reference<>())
      .add_property("ordinal_id_subject_sequence", &ncbi_blast_BlastHSPList_wrapper::get_ordinal_id_subject_sequence,
		    "The ordinal id of the subject sequence this HSP list is for.")
      .add_property("query_index", &ncbi_blast_BlastHSPList_wrapper::get_query_index,
		    "Index of the query which this HSPList corresponds to. Set to 0 if not applicable")
      .add_property("best_evalue", &ncbi_blast_BlastHSPList_wrapper::get_best_evalue,
		    "Smallest e-value for HSPs in this list. Filled after e-values are calculated.\
 Necessary because HSPs are sorted by score, but highest scoring HSP may not have the lowest \
 e-value if sum statistics is used")
      ;

  class_< BlastHSP,
    ncbi_blast_BlastHSP_wrapper>("BlastHSP",
				 "C++ glue to access Blast HSP info (the actal hit!)",
				 init<>("Initializes with junk, most likely..")
				 )
    .def("get_gap_info", &ncbi_blast_BlastHSP_wrapper::get_gap_info, return_internal_reference<>(),
	 "Get the gap editing script containing all information returned from a  single gapped extension")

    .add_property("score", &ncbi_blast_BlastHSP_wrapper::get_score,
		  "This HSP's raw score ")
    .add_property("num_ident", &ncbi_blast_BlastHSP_wrapper::get_num_ident,
		  "Number of identical base pairs in this HSP ")
    .add_property("bit_score", &ncbi_blast_BlastHSP_wrapper::get_bit_score,
		  "Bit score, calculated from score ")
    .add_property("evalue", &ncbi_blast_BlastHSP_wrapper::get_evalue,
		  "This HSP's e-value ")
    .add_property("query", &ncbi_blast_BlastHSP_wrapper::get_query,
		  "Query sequence info. A tuple: \n\
 (<Translation frame>, <Start of hsp>, <End of HSP>, <Where the gapped extension started>)")
    .add_property("subject", &ncbi_blast_BlastHSP_wrapper::get_subject,
		  "Subject sequence info. A tuple: \n\
 (<Translation frame>, <Start of hsp>, <End of HSP>, <Where the gapped extension started>)")
    .add_property("context", &ncbi_blast_BlastHSP_wrapper::get_context,
		  "Context number of query")
    .add_property("num_linked_hsp", &ncbi_blast_BlastHSP_wrapper::get_num_linked_hsp,
		  "How many HSP's are linked together for sum statistics evaluation? If unset (0), this HSP is\n\
 not part of a linked set, i.e. value 0 is treated the same way as 1.")
    ;
  
  class_< gap_edit_script_iterator >("gap_edit_script_iterator",
                                     "Iterates along the elements of a GapEditScript",
                                     init<const GapEditScript*>("Initializes root pointer..")
                                     )
    .def("__iter__", &gap_edit_script_iterator::get_iter)
    .def("next",     &gap_edit_script_iterator::get_next)
    ;


  class_< GapEditScript, 
    ncbi_blast_GapEditScript_wrapper>("GapEditScript",
                                      "C++ glue to access Blast GapEditScript info.",
                                      init<>("Initializes with junk, most likely..")
                                      )
    .def("__len__",  &ncbi_blast_GapEditScript_wrapper::get_len)
    .def("__iter__", &ncbi_blast_GapEditScript_wrapper::get_iterator)
    ;
}



