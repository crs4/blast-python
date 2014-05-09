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

import sys
from timeit import Timer


BUFSIZE_EXP_MIN = 1  # min. bufsize = 10**BUFSIZE_EXP_MIN
BUFSIZE_EXP_MAX = 9  # max. bufsize = 10**BUFSIZE_EXP_MAX
NRUNS = 1000


def line_length(filename):
    f = open(filename)
    maxlen = 0
    lensum = 0
    for i, line in enumerate(f):
        L = len(line)
        if L > maxlen:
            maxlen = L
        lensum += L
    f.close()
    return maxlen, lensum/float(i)


def main(argv):
    try:
        filename = sys.argv[1]
    except IndexError:
        print "Usage: %s FILENAME" % sys.argv[0]
        sys.exit(2)

    IMPORT = "from BlastPython.buffered_IO import buffered_IO"

    CODE = """f = open("%s")
s = buffered_IO(f, %d)
l = s.readline()
while l:
    l = s.readline()
f.close()
"""
    
    print "max/average line length: %d/%d" % line_length(filename)
    print "-" * 79
    bufsizes = [10**i for i in xrange(BUFSIZE_EXP_MIN, BUFSIZE_EXP_MAX)] + [0]
    for bs in bufsizes:
        print "%9d:" % bs, Timer(CODE % (filename, bs), IMPORT).timeit(NRUNS)


if __name__ == "__main__":
    main(sys.argv)
