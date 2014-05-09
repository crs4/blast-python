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
#include <algo/blast/api/blast_options.hpp>
#include <algo/blast/core/blast_options.h> // for the default values defines
#include <ncbi_source_ver.h> // get NCBI_PRODUCTION_VER

#include <boost/python.hpp>

#include <assert.h>
#include <map>

using namespace boost::python;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Helpers
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

struct get_holder_base {
  typedef ncbi::blast::CBlastOptions target_class_type;
  virtual object operator()(const target_class_type* ci) = 0;
};
struct set_holder_base {
  typedef ncbi::blast::CBlastOptions target_class_type;
  virtual void operator()(target_class_type* ci, const object& v) = 0;
  virtual void operator()(target_class_type* ci) = 0;
};

template<typename R>
struct get_holder : get_holder_base {
  typedef R (target_class_type::*method_type)();
  
  get_holder(method_type m) : _m(m) {}
  object operator()(const target_class_type* ci) {
    target_class_type* n_ci = const_cast<target_class_type*>(ci);
    R v = (n_ci->*_m)();
    return object(v);
  }
  method_type _m;
};

template<typename R>
struct set_holder : set_holder_base {
  typedef R (target_class_type::*method_type)(R);
  
  set_holder(method_type m, R v) : _m(m), _v(v) {}

  void operator()(target_class_type* ci, const object& ov) {
    R v = extract<R>(ov);
    (ci->*_m)(v);
  }
  void operator()(target_class_type* ci) {
    (ci->*_m)(_v);
  }

  method_type _m;
  R           _v;
};


typedef struct command_s {
  std::string name;
  std::string doc;
  set_holder_base* setter;
  get_holder_base* getter;
  command_s (const char* n, const char* d, set_holder_base* s, get_holder_base* g) : 
    name(n), doc(d), setter(s), getter(g) {}
  ~command_s () {
    delete setter;
    delete getter;
  }
} command_t;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Option list definitions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define get_method_sig(r,o,m) (r (o::*)())(&o::m)
#define set_method_sig(r,o,m) (r (o::*)(r))(&o::m)
#define set_method(r,m,d) new set_holder<r>(set_method_sig(r,ncbi::blast::CBlastOptions,m), d)
#define get_method(r,m)   new get_holder<r>(get_method_sig(r,ncbi::blast::CBlastOptions,m))
#define define_command(n,d,t,v) command_s(#n,d,set_method(t, Set ## n, v), get_method(t, Get ## n))


#include "blast_options_list.incl"

//namespace { // put here a namespace if needed. This is only to avoid indenting if we change idea...

  typedef std::map<std::string, const command_t*> command_map_t;

  struct ncbi_blast_CBlastOptions_wrapper : public ncbi::blast::CBlastOptions, wrapper<ncbi::blast::CBlastOptions> {

    ncbi_blast_CBlastOptions_wrapper(PyObject* py_self, 
				     ncbi::blast::CBlastOptions::EAPILocality p0 = ncbi::blast::CBlastOptions::eLocal,
				     const char* progname = "blastn" ):
      ncbi::blast::CBlastOptions(p0), _py_self(py_self) {
	build_g_command_map_if_needed();
	for(command_map_t::iterator it = g_command_map.begin(); it != g_command_map.end(); ++it){
	  (*(it->second->setter))(this);
	}
      }

    int len() const {
      return(g_command_map.size());
    }

    bool has_key(str key) const { 
      std::string k = extract<std::string>(key);
      return g_command_map.find(k) != g_command_map.end();
    }
    
    object get_item(str key) const {
      std::string k = extract<std::string>(key);
      if (g_command_map.find(k) != g_command_map.end()) {
	const command_t* com = g_command_map[k];
	assert( com->name == k);
	return (*(com->getter))(this);
      } else {  
	PyErr_SetString(PyExc_KeyError, "Invalid key");
	throw_error_already_set();
      }
    }
    void set_item(str key, object v) {
      std::string k = extract<std::string>(key);
      if (g_command_map.find(k) != g_command_map.end()) {
	const command_t* com = g_command_map[k];
	assert( com->name == k);
	return (*(com->setter))(this, v);
      } else {  
	PyErr_SetString(PyExc_KeyError, "Invalid key");
	throw_error_already_set();
      }
    }
    
    list keys() const {
      list l;
      for(command_map_t::iterator it = g_command_map.begin(); it != g_command_map.end(); ++it){
	l.append(str(it->first));
      }
      return l;
    }


    static bool has_key_static(const ncbi::blast::CBlastOptions* p, str k) {
      build_g_command_map_if_needed();
      const ncbi_blast_CBlastOptions_wrapper* cp = static_cast<const ncbi_blast_CBlastOptions_wrapper*>(p);
      return cp->has_key(k);
    }

    static object get_item_static(const ncbi::blast::CBlastOptions* p, str k) {
      build_g_command_map_if_needed();
      const ncbi_blast_CBlastOptions_wrapper* cp = static_cast<const ncbi_blast_CBlastOptions_wrapper*>(p);
      return cp->get_item(k);
    }
    static void set_item_static(ncbi::blast::CBlastOptions* p, str k, object v) {
      build_g_command_map_if_needed();
      ncbi_blast_CBlastOptions_wrapper* cp = static_cast<ncbi_blast_CBlastOptions_wrapper*>(p);
      cp->set_item(k, v);
    }
    static list keys_static(const ncbi::blast::CBlastOptions* p) {
      build_g_command_map_if_needed();
      const ncbi_blast_CBlastOptions_wrapper* cp = static_cast<const ncbi_blast_CBlastOptions_wrapper*>(p);
      return cp->keys();
    }

    static void build_g_command_map_if_needed() {
      if (g_command_map.size() == 0 && sizeof(commands) > 0) {
	for(std::size_t i = 0; i < sizeof(commands)/sizeof(command_t); ++i){
	  g_command_map.insert(std::make_pair(commands[i].name, &(commands[i])));
	}
      }
    }
    static command_map_t  g_command_map;

    PyObject* _py_self;
  };

  command_map_t  ncbi_blast_CBlastOptions_wrapper::g_command_map;



void export_blast_options() 
{
    enum_< ncbi::blast::CBlastOptions::EAPILocality>("EAPILocality")
      .value("eLocal", ncbi::blast::CBlastOptions::eLocal)
      .value("eRemote", ncbi::blast::CBlastOptions::eRemote)
      .value("eBoth", ncbi::blast::CBlastOptions::eBoth)
      ;

    enum_ < ncbi::blast::EProgram>("EProgram")
      .value("eBlastn", ncbi::blast::eBlastn)
      .value("eBlastp", ncbi::blast::eBlastp)
      .value("eBlastx", ncbi::blast::eBlastx)
      .value("eTblastn", ncbi::blast::eTblastn)
      .value("eTblastx", ncbi::blast::eTblastx)
      .value("eRPSBlast", ncbi::blast::eRPSBlast)
      .value("eRPSTblastn", ncbi::blast::eRPSTblastn)
      .value("eMegablast", ncbi::blast::eMegablast)
      .value("eDiscMegablast", ncbi::blast::eDiscMegablast)
      .value("ePSIBlast", ncbi::blast::ePSIBlast)
      .value("eBlastProgramMax", ncbi::blast::eBlastProgramMax)
      ;

    enum_ < EBlastTbackExt>("EBlastTbackExt")
      .value("eDynProgTbck", eDynProgTbck)
      .value("eGreedyTbck",  eGreedyTbck)
      .value("eSmithWatermanTbck", eSmithWatermanTbck)
      ;

    /*
    enum_< EBlastPrelimGapExt > ("EBlastPrelimGapExt")
      .value("eDynProgExt", eDynProgExt)
      .value("eGreedyExt", eGreedyExt)
      .value("eGreedyWithTracebackExt", eGreedyWithTracebackExt)
      ;
    */

    enum_< ncbi::objects::ENa_strand >("strand")
      .value("both_rev", ncbi::objects::eNa_strand_both_rev)
      .value("plus", ncbi::objects::eNa_strand_plus)
      .value("minus", ncbi::objects::eNa_strand_minus)
      .value("unknown", ncbi::objects::eNa_strand_unknown)
      .value("other", ncbi::objects::eNa_strand_other)
      .value("both", ncbi::objects::eNa_strand_both)
      ;    

    class_< ncbi::blast::CBlastOptions,
    boost::noncopyable, 
    ncbi_blast_CBlastOptions_wrapper
    >("CBlastOptions", 
	"C++ glue to ncbi::blast::CBlastOptions",
	init<>("Initializes with blastn default values")
	)
      .def("__len__", &ncbi_blast_CBlastOptions_wrapper::len)
      .def("__getitem__", &ncbi_blast_CBlastOptions_wrapper::get_item)
      .def("__getitem__", &ncbi_blast_CBlastOptions_wrapper::get_item_static)
      .def("__setitem__", &ncbi_blast_CBlastOptions_wrapper::set_item)
      .def("__setitem__", &ncbi_blast_CBlastOptions_wrapper::set_item_static)
      .def("keys",   &ncbi_blast_CBlastOptions_wrapper::keys)
      .def("keys",   &ncbi_blast_CBlastOptions_wrapper::keys_static)
      .def("has_key",  &ncbi_blast_CBlastOptions_wrapper::has_key)
      .def("has_key",  &ncbi_blast_CBlastOptions_wrapper::has_key_static)
      .def("contains", &ncbi_blast_CBlastOptions_wrapper::has_key)
      .def("contains", &ncbi_blast_CBlastOptions_wrapper::has_key_static)
      .def("validate", &ncbi::blast::CBlastOptions::Validate)
      ;
}


//}
