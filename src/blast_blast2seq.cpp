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

// 
#define private protected
#include <algo/blast/api/bl2seq.hpp>
#undef private


#include <boost/python.hpp>

#include <assert.h>

using namespace boost::python;

// for compatibility with old boost versions
#ifndef len
#define len(l) PyList_Size(l.ptr())
#endif

struct ncbi_blast_CBl2Seq_wrapper : public ncbi::blast::CBl2Seq {
  typedef ncbi::blast::CBl2Seq super_t;

  ncbi_blast_CBl2Seq_wrapper(PyObject* py_self, 
			     const ncbi::blast::SSeqLoc& query, const ncbi::blast::SSeqLoc& subject,
			     const ncbi::blast::EProgram prog) :
    super_t(query, subject, prog), _py_self(py_self) {}


  list get_results() const {
    list l;
    BlastHSPResults* r = GetResults();
    for(int i = 0; i < r->num_queries; ++i){
      l.append(r->hitlist_array[i]);
    }
    return l;
  }
  list get_results_static(const ncbi::blast::CBl2Seq* p) const {
    const ncbi_blast_CBl2Seq_wrapper* cp = static_cast<const ncbi_blast_CBl2Seq_wrapper*>(p);
    return cp->get_results();
  }

  void SetupSearchProxy(bool query_already_setup = false) {   
    this->mi_bQuerySetUpDone = query_already_setup;
    this->SetupSearch();
  }

  void SetQueryProxy(const ncbi::blast::SSeqLoc& query){
    SetQuery(query);
  }
  void SetQueriesProxy(list queries){
    size_t L = len(queries);
    ncbi::blast::TSeqLocVector vec;
    vec.reserve(L);
    for(std::size_t i = 0; i < L; ++i){
      extract<const ncbi::blast::SSeqLoc&> x(queries[i]);
      if (x.check()){
	const ncbi::blast::SSeqLoc& s = x();
	vec.push_back(s);
      } else {
	throw std::runtime_error("Illegal list element.");
      }
    }
#if 1
    this->x_ResetQueryDs();
//     std::cerr << "Done with reset" << std::endl;
    this->m_tQueries.clear();
    this->m_tQueries.reserve(vec.size());
    for(std::size_t i = 0; i < vec.size(); ++i) {
      this->m_tQueries.push_back(vec[i]);
    }
#else
    SetQueries(vec);
#endif
  }

  void SetSubjectProxy(const ncbi::blast::SSeqLoc& subject){
    SetSubject(subject);
  }
  void SetSubjectsProxy(list subjects){
    size_t L = len(subjects);
    ncbi::blast::TSeqLocVector vec;
    vec.reserve(L);
    for(std::size_t i = 0; i < L; ++i){
      extract<const ncbi::blast::SSeqLoc&> x(subjects[i]);
      if (x.check()){
	const ncbi::blast::SSeqLoc& s = x();
	vec.push_back(s);
      } else {
	throw std::runtime_error("Illegal list element.");
      }
    }
#if 1
    this->x_ResetSubjectDs();
//     std::cerr << "Done with reset" << std::endl;
    this->m_tSubjects.clear();
    this->m_tSubjects.reserve(vec.size());
    for(std::size_t i = 0; i < vec.size(); ++i) {
      this->m_tSubjects.push_back(vec[i]);
    }
#else
    SetSubjects(vec);
#endif
  }

  void ScanDBProxy() {   
    this->RunFullSearch();  // TODO: check why this is not called directly.
  }

  const ncbi::blast::CBlastOptions& GetOptionsProxy() const {   
    return (this->GetOptionsHandle()).GetOptions();
  }

  ncbi::blast::CBlastOptions& SetOptionsProxy() {   
    return (this->SetOptionsHandle()).SetOptions();
  }

  
  PyObject* _py_self;
};

static char set_query_doc[] = "Set the query to a given SSeqLoc";
static char get_query_doc[] = "Get the query";
static char set_subject_doc[] = "Set the subject to a given SSeqLoc";
static char get_subject_doc[] = "Get the subject";
static char set_options_doc[] = "Returns a modifiable CBlastOptions object.";
static char get_options_doc[] = "Returns a NON modifiable CBlastOptions object.";

static char setup_search_doc[] = "Process the queries, do setup, and build the lookup table. If it passed True, assumes that query SSeqLoc has alreaby been setup. WARNING pass False the first time around :-).";
static char RunWithoutSeqalignGeneration_doc[] = "Runs the search but does not produce seqalign output (useful if the raw search results are needed, rather than a set of complete Seq-aligns) (it amounts to  SetupQuery(False) + ScanDB())";
static char get_diagnostics_doc[] = "Retrieves the diagnostics information returned from the engine.";
static char get_results_doc[] = "Retrieves the list of HSP results from the engine (to be used after RunWithoutSeqalignGeneration method)";
static char get_error_message[] = "Returns error messages/warnings.";

void export_blast_blast2seq() 
{
  class_< ncbi::blast::CBl2Seq, boost::noncopyable, 
    ncbi_blast_CBl2Seq_wrapper>("CBl2Seq", 
				"C++ glue to ncbi::blast::CBl2Seq",
				init<const ncbi::blast::SSeqLoc&, 
				const ncbi::blast::SSeqLoc&,
				const ncbi::blast::EProgram>("Initializes as an empty sequence.")
				)
    .def("SetQuery", &ncbi_blast_CBl2Seq_wrapper::SetQueryProxy)
    .def("SetQuery", &ncbi_blast_CBl2Seq_wrapper::SetQueriesProxy)
    .def("GetQuery", &ncbi::blast::CBl2Seq::GetQuery, return_internal_reference<>())
    .def("SetSubject", &ncbi_blast_CBl2Seq_wrapper::SetSubjectProxy)
    .def("SetSubject", &ncbi_blast_CBl2Seq_wrapper::SetSubjectsProxy)
    .def("GetSubject", &ncbi::blast::CBl2Seq::GetSubject, return_internal_reference<>())
    
    .def("SetOptions", &ncbi_blast_CBl2Seq_wrapper::SetOptionsProxy, 
         return_internal_reference<>())
    .def("GetOptions", &ncbi_blast_CBl2Seq_wrapper::GetOptionsProxy, 
         return_internal_reference<>())
    .def("SetupSearch", &ncbi_blast_CBl2Seq_wrapper::SetupSearchProxy)
    .def("ScanDB", &ncbi_blast_CBl2Seq_wrapper::ScanDBProxy)
    .def("RunWithoutSeqalignGeneration", &ncbi::blast::CBl2Seq::RunWithoutSeqalignGeneration)
    .def("GetDiagnostics", &ncbi::blast::CBl2Seq::GetDiagnostics, return_internal_reference<>())
    .def("GetResults", &ncbi_blast_CBl2Seq_wrapper::get_results)
    .def("GetMessages",  &ncbi::blast::CBl2Seq::GetMessages, return_value_policy<manage_new_object>())
    ;
}
