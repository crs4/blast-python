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
import time

class blast_result_stream(object):

  def __init__(self, blaster, seq_stream):
    self.blaster = blaster
    self.in_stream  = seq_stream
    self.total_time = 0

  def __iter__(self) :
    return self

  def next(self):
    seq = self.in_stream.next()
    start = time.time()
    r = self.blaster.blast(seq)
    self.total_time += time.time() - start
    return r
