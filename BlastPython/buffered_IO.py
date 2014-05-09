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

"""
buffered_IO adds readline to IO objects that only support read and write (e.g.
srb streams). If your IO objects already have a readline method (e.g. ordinary
python file objects) you can fall back to python's readline by setting bufsize
to 0.

*IMPORTANT*: if you do need to use this wrapper with bufsize != 0, make sure
you set it to a value at least equal to maximum input line size! This is the
normal use case which the wrapper is designed to optimize (buffered data is
kept in a string at all times, therefore you could get a performance-killing
iterative string join with very long lines).
"""


class buffered_IO:

  def __init__(self, IOint, bufsize):
    self.fp = IOint
    self.bufsize = bufsize
    self.data = ''
    self.offset = 0

  def readline(self):
    if self.bufsize <= 0:
      return self.fp.readline()
    eol = self.data.find('\n', self.offset)
    if eol == -1:
      new_chunk = self.fp.read(2*self.bufsize)
      self.data = self.data[self.offset:] + new_chunk
      self.offset = 0
      if not new_chunk:  # last line
        line, self.data = self.data, ''
        return line
      if self.data:
        return self.readline()
      else:
        return ''
    else:
      st, self.offset = self.offset, eol+1
      return self.data[st:eol+1]
