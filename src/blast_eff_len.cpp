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
#include <algo/blast/core/blast_stat.h>
#include <boost/python.hpp>

//+++ Helper Class +++
class blast_eff_len {
public:
  blast_eff_len() {}

  Int8 calc_eff_len(const char* matrix_name, Int4 query_length, Int8 db_length, Int4 db_num_seqs) {

    return query_length + db_num_seqs + db_length;  // DUMMY - real code TBD
    
//     Int4 length_adjustment = 0;

//     // code for k, logK, alpha_d_lambda, beta

//     BLAST_ComputeLengthAdjustment(K,
// 				  logK,
// 				  alpha_d_lambda,
// 				  beta,
// 				  query_length,
// 				  db_length,
// 				  db_num_seqs,
// 				  &length_adjustment);

//     return db_length - (Int8)db_num_seqs * length_adjustment;
    
  }
};
//++++++++++++++++++++

using namespace boost::python;

struct blast_eff_len_wrapper : public blast_eff_len, wrapper<blast_eff_len> {

  blast_eff_len_wrapper(PyObject* py_self): _py_self(py_self) {}

  long_ calc(str matrix_name, int query_length, long_ db_length, int db_num_seqs) {
    std::string mn = extract<std::string>(matrix_name);
    Int8 dbl = extract<Int8>(db_length);

    Int8 res = calc_eff_len(mn.c_str(), query_length, dbl, db_num_seqs);
    return long_(res);
  }

  PyObject* _py_self;
};

void export_blast_eff_len() {
  class_< blast_eff_len, boost::noncopyable, blast_eff_len_wrapper
    >("blast_eff_len", "Helper class for computing db effective length")
    .def("calc", &blast_eff_len_wrapper::calc)
    ;
}
