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
Drivers to access data files -- standard file system access.
"""

from buffered_IO import buffered_IO


class FS_IO(file):

    def __init__(self, filename, bufsize):
        file.__init__(self,filename)
        self.bufsize = bufsize
        self.buffer_obj = buffered_IO(self, bufsize)

    def next(self):
        if not self.bufsize:
            return file.next(self)
        else:
            line = self.buffer_obj.readline()
            if line:
                return line
            else:
                raise StopIteration
