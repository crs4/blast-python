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
#ifndef _CSEQ_SEQUENCE_EXTRACTOR_HPP_
#define _CSEQ_SEQUENCE_EXTRACTOR_HPP_

#include <string>
#include <objmgr/bioseq_handle.hpp>

#include <objects/seq/Bioseq.hpp>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(objects)
BEGIN_SCOPE(sequence)

std::string extract_sequence_from_bioseq(const ncbi::objects::CBioseq_Handle& handle,
					 bool instantiate_gaps,
					 int beg, int end);

END_SCOPE(sequence)
END_SCOPE(objects)
END_NCBI_SCOPE


#endif // _CSEQ_SEQUENCE_EXTRACTOR_HPP_
