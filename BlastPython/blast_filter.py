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
import logging


class base_blast_filter(object):

  def __init__(self, in_stream):
    self.in_stream = in_stream

  def __iter__(self):
    return self

  def is_acceptable(self, r):
    return True

  def next(self):
    while True:
      r=self.in_stream.next()
      if self.is_acceptable(r):
        return r


class blast_filter(base_blast_filter):
  """
  Pass all filter implementation. Saves in .total_time the total time
  spent, and in .counter the number of objects counted. A value of 0 for
  max_count means infinite.
  """
  def __init__(self, in_stream, max_count = 0):
    base_blast_filter.__init__(self, in_stream)
    self.max_count  = max_count
    self.total_time = 0
    self.counter = 0
    self.start_time = time.time()
    self.logger = logging.getLogger('blast_filter')

  def is_acceptable(self, r):
    self.counter += 1
    self.logger.debug("doing %d: %s " % (self.counter, r[0].id))
    if self.counter > self.max_count > 0:
      raise StopIteration
    delta_t = time.time() - self.start_time
    self.total_time += delta_t
    self.logger.debug("delta t: %g" % delta_t)
    self.start_time = time.time()
    return True
