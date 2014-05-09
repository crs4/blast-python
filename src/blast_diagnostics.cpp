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
#include <algo/blast/core/blast_diagnostics.h>
#include <boost/python.hpp>

#include <assert.h>

using namespace boost::python;

struct ncbi_blast_BlastDiagnostic_wrapper : BlastDiagnostics {
  typedef BlastDiagnostics super_t;

  ncbi_blast_BlastDiagnostic_wrapper(PyObject*& py_self, const BlastDiagnostics& bd) : 
    BlastDiagnostics(bd), _py_self(py_self){}

  ncbi_blast_BlastDiagnostic_wrapper(PyObject* py_self) : _py_self(py_self) {}

  static tuple get_ungapped_stats(const super_t& s) {
    return make_tuple(s.ungapped_stat->lookup_hits, s.ungapped_stat->num_seqs_lookup_hits,
		      s.ungapped_stat->init_extends, s.ungapped_stat->good_init_extends,
		      s.ungapped_stat->num_seqs_passed);
  }

  static tuple get_gapped_stats(const super_t& s) {
    return make_tuple(s.gapped_stat->seqs_ungapped_passed, /* s.gapped_stat->extra_extensions, */
		      s.gapped_stat->extensions, s.gapped_stat->good_extensions,
		      s.gapped_stat->num_seqs_passed);
  }
  
  static tuple get_raw_cutoffs(const super_t& s) {
    return make_tuple(s.cutoffs->x_drop_ungapped, s.cutoffs->x_drop_gap,
		      s.cutoffs->x_drop_gap_final, s.cutoffs->ungapped_cutoff, 
		      s.cutoffs->cutoff_score);
  }

  PyObject* _py_self;
};

static const char get_ungapped_stats_docs[] = "\
Statistics on the hit counts from the ungapped stage of a BLAST  search. Returns the tuple\n\
(<Number of successful lookup table hits>, \n\
 <Number of sequences which had at least one lookup table hit.>,\n\
 <Number of initial words found and extended>,\n\
 <Number of successful initial extensions, i.e. number of HSPs saved after ungapped stage.>,\n\
 <Number of sequences with at least one HSP saved after ungapped stage>)";

static const char get_gapped_stats_docs[] = "\
Statistics on the hit counts from the gapped stage of a BLAST  search. Returns the tuple\n\
(<Number of sequences with top HSP after ungapped extension passing the e-value threshold.>,\n\
 <Total number of gapped extensions performed>,\n\
 <Number of HSPs below the e-value threshold after gapped extension>,\n\
 <Number of sequences with top HSP passing the e-value threshold.>)";

static const char get_raw_cutoffs_docs[] = "\
Raw cutoff and gap-x-drop values. It returns the tuple\n\
(<Raw value of the x-dropoff for ungapped extensions>,\n\
 <Raw value of the x-dropoff for preliminary gapped extensions>,\n\
 <Raw value of the x-dropoff for gapped extensions with traceback>,\n\
 <Minimal raw score for starting gapped extension>,\n\
 <Cutoff score corresponding to given evalue.>)";

void export_blast_diagnostics()
{
  class_< BlastDiagnostics, 
    ncbi_blast_BlastDiagnostic_wrapper>("BlastDiagnostics", 
					"C++ glue to BlastDiagnostics",
					init<>("Initializes as an empty container.")
					)
    .def("get_ungapped_stats", &ncbi_blast_BlastDiagnostic_wrapper::get_ungapped_stats, get_ungapped_stats_docs)
    .def("get_gapped_stats",   &ncbi_blast_BlastDiagnostic_wrapper::get_gapped_stats, get_gapped_stats_docs)
    .def("get_raw_cutoffs",    &ncbi_blast_BlastDiagnostic_wrapper::get_raw_cutoffs, get_raw_cutoffs_docs)
    ;
}


