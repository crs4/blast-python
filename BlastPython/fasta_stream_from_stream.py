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

from fasta_stream import fasta_stream


class fasta_stream_from_stream(fasta_stream):
    
    def __init__(self, io_obj):
        """
        Generate a stream of 'FASTA strings' from an io stream.
        """
        self.IO = io_obj
        self.buffer = []

    def __del__(self) :
        self.IO.close()

    def next(self):
        while 1:
            try:
                l = self.IO.next()
                if l.startswith('>') and self.buffer:
                    fasta = "".join(self.buffer)
                    self.buffer = [l]
                    return fasta
                else:
                    self.buffer.append(l)
            except StopIteration:
                if self.buffer:
                    fasta = "".join(self.buffer)
                    self.buffer = []
                    return fasta
                else:
                    raise StopIteration
