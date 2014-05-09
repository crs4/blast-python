# BEGIN_COPYRIGHT
# 
# Copyright (C) 2014 CRS4.
# 
# This file is part of blast-python.
# 
# blast-python is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# blast-python is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
# 
# You should have received a copy of the GNU General Public License along with
# blast-python.  If not, see <http://www.gnu.org/licenses/>.
# 
# END_COPYRIGHT
import ncbi_toolkit


class seq_factory_from_fasta(ncbi_toolkit.blast_sseq_loc_from_fasta):

  def __init__(self, strand) :
    ncbi_toolkit.blast_sseq_loc_from_fasta.__init__(self)
    self.strand   = strand

  def make(self, fasta):
    return super(seq_factory_from_fasta, self).make(
      fasta, self.strand, 0, 0, False
      )


class seq_factory_from_str(ncbi_toolkit.blast_sseq_loc_from_str) :

  def __init__(self, strand) :
    ncbi_toolkit.blast_sseq_loc_from_str.__init__(self)
    self.strand   = strand

  def make(self, s):
    i, seq_data = s.split()[:2]
    return super(seq_factory_from_str, self).make(
      seq_data, False, int(i), 'fake title', self.strand, 0, 0
      )
