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

import random, time
import ncbi_toolkit


NITER = 1000
LEN = 10000


r = random.Random()

def make_seq(n):
    return ''.join([r.choice(['A', 'C', 'G', 'T']) for i in xrange(n)])

s = make_seq(LEN)
factory_fasta = ncbi_toolkit.blast_sseq_loc_from_fasta()
factory_str = ncbi_toolkit.blast_sseq_loc_from_str()
start = time.time()
for i in xrange(NITER):
    sseq = factory_fasta.make(
        '>xxxx\n%s' % s, ncbi_toolkit.strand.plus, 0, 0, False
        )
print 'sseq construction (Fasta) ', (time.time() - start)/NITER

start = time.time()
for i in xrange(NITER):
    sseq = factory_str.make(
        s, False, 10022, 'title xxx', ncbi_toolkit.strand.plus, 0, 0
        )
print 'sseq construction (str)', (time.time() - start)/NITER

start = time.time()
for i in xrange(NITER):
    sseq = factory_fasta.make_dummy(
        '>xxxx\n%s' % s, ncbi_toolkit.strand.plus, 0, 0, False
        )
print 'dummy construction ', (time.time() - start)/NITER

start = time.time()
for i in xrange(NITER):
    sseq = 'x' * LEN
print 'string construction ', (time.time() - start)/NITER
