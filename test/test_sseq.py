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

import unittest, random
import ncbi_toolkit

N = 10
M = 10
GID_OFFSET = 1000


def make_seq(n):
    return ''.join([random.choice(['A', 'C', 'G', 'T']) for i in xrange(n)])


class sseq_tc(unittest.TestCase):
    def test_id(self):
        for t, s in self.sseqs:
            self.assertEqual('gi|%d' % t[2], s.id)

    def test_title(self):
        for t, s  in self.sseqs:
            self.assertEqual(t[3], s.title)

    def test_len_by_len(self):
        for t, s  in self.sseqs:
            self.assertEqual(len(t[0]), len(s))

    def test_len_by_length(self):
        for t, s  in self.sseqs:
            self.assertEqual(len(t[0]), s.length)

    def test_seq(self):
        for t, s  in self.sseqs:
            self.assertEqual(t[0], s.get_sequence())

    def test_strand(self):
        for t, s  in self.sseqs:
            self.assertEqual(t[4], s.strand)

    def test_slice(self):
        for t, s  in self.sseqs:
            x = t[0]
            self.assertEqual(x, s[:])
            for i in range(len(x)):
                self.assertEqual(x[0:i], s[0:i])
                self.assertEqual(x[i:], s[i:])
                self.assertEqual(x[:-i], s[:-i])


class sseq_from_str_tc(sseq_tc):
    def setUp(self):
        factory = ncbi_toolkit.blast_sseq_loc_from_str()
        self.sseqs = []
        for i in xrange(M):
            t = (make_seq(N),  False, GID_OFFSET + i, 'title %s' % i,
                 random.choice(
                     [ncbi_toolkit.strand.plus, ncbi_toolkit.strand.minus]
                     ), 0, 0)
            self.sseqs.append(
                (t, factory.make(t[0], t[1], t[2], t[3], t[4], t[5], t[6]))
                )


class sseq_from_fasta_tc(sseq_tc):
    def setUp(self):
        factory = ncbi_toolkit.blast_sseq_loc_from_fasta()
        self.sseqs = []
        for i in xrange(M):
            t = (make_seq(N), False, GID_OFFSET + i, 'title %s' % i,
                 random.choice(
                     [ncbi_toolkit.strand.plus, ncbi_toolkit.strand.minus]
                     ), 0, 0)
            self.sseqs.append((t, factory.make(
                '>gi|%s %s \n%s' % (t[2], t[3], t[0]), t[4], t[5], t[6], False
                )))


def suite():
    suite = unittest.TestSuite()
    suite.addTest(sseq_from_str_tc('test_id'))
    suite.addTest(sseq_from_str_tc('test_title'))
    suite.addTest(sseq_from_str_tc('test_len_by_len'))
    suite.addTest(sseq_from_str_tc('test_len_by_length'))
    suite.addTest(sseq_from_str_tc('test_seq'))
    suite.addTest(sseq_from_str_tc('test_strand'))
    suite.addTest(sseq_from_str_tc('test_slice'))
    #--
    suite.addTest(sseq_from_fasta_tc('test_id'))
    suite.addTest(sseq_from_fasta_tc('test_title'))
    suite.addTest(sseq_from_fasta_tc('test_len_by_len'))
    suite.addTest(sseq_from_fasta_tc('test_len_by_length'))
    suite.addTest(sseq_from_fasta_tc('test_seq'))
    suite.addTest(sseq_from_fasta_tc('test_strand'))
    suite.addTest(sseq_from_fasta_tc('test_slice'))
    #---
    return suite


if __name__ == '__main__':
  runner = unittest.TextTestRunner(verbosity=2)
  runner.run((suite()))
