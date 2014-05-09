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
Drivers to access data files -- SRB access.
"""

import srb
from buffered_IO import buffered_IO


class SRB_IO:
    def __init__(self, path, filename, host, port, domain, user, passwd,
                 bufsize):
        self.conn = srb.connect(
            host, port, domain, "ENCRYPT1", user, passwd, ""
            )
        self.fd = srb.obj_open(self.conn, path, filename, 0)
        self.buffer_obj = buffered_IO(self, bufsize)

    def __iter__(self):
        return self

    def read(self, size):
        return srb.obj_read(self.conn, self.fd, size)

    def write(self):
        pass

    def close(self):
        srb.disconnect(self.conn)

    def next(self):
        line = self.buffer_obj.readline()
        if line:
            return line
        else:
            raise StopIteration
