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

import itertools, unittest
from BlastPython import *
from ncbi_toolkit import *

seq_tuples = [
    ('gi|67678568|gb|DR109384.1|DR109384',
     "RTS1_1_D03.g1_A029 Roots minus sulfur Pinus taeda cDNA clone RTS1_1_D03_A029 5', mRNA sequence",
     strand.plus,
     "ACCGAAGAGGAGACTGAAATGAATAAAACTTCTCCCGTGTCAGAGAATGTAACTTTAGAGAAGAATGTGTCAAATTATTTGGAAGTAGTCCCACAGCCATCGTTATGTACTGGAACACAAAATGAACAGCACAAGCAGATGGAAGACTCCCCTGCAAGGGCAAACTGTGACAGCTCGGCAGGCAATATACAGGTACGAGATAAATATGAGGAACTGCAAACCGAAGAGGAGATTGAAATGAAAAAAGCTTCTCCTGTTTCAGATAATGTAACTCTAGAGATTAATGTCCCGAAGTGTTTGGAAGTCGTCCCACTGCCATCGATATGCACTGAAATACGAGATGAACGGCACAAGCAGCAGATGGAAGACTCACCTGTAAGGGAAATCTGTGACAGCTCAGCAGGCATTAAACATGTACAAGATAAATGTGAGGAACTGCGTAATGAAACTCGGGCTATATGTGTGGAAATGTGGCAGGACTGTGAGGAAAAGTTTCAAGACATGCAGTCTAGGATTAATGAAATGCTGGATAACAAGATGAAGGAGTTAAAGGAAACTCACGCTATCCAGGACAGCTGTTTTAATGAAAGGAGACTGTTCTTTGGAGAGGAGTCTTCTGCTAGCAATGGCAGTCGTGAAATCCATCTTTTAATGGACAAACATCT"),
    ('lcl|BM970451',
     "UI-CF-EC1-abs-k-20-0-UI.s1 UI-CF-EC1 Homo sapiens cDNA clone UI-CF-EC1-abs-k-20-0-UI 3', mRNA sequence. #CRS4 tissue_type-Lung- organ-Lung-#",
     strand.plus,
     "TTTTTTTTTTTTCTTTTTCACGCATTTGCTTTATTCGAAAAGAGGCTTTTAAAATGTGCATGTTTAGAAACAAAATTTCTTCATGGAAATCATATACATTAGAAAATCACAGTCAGATGTTTAATCAATCCAAAATGTCCACTATTTCTTATGTCATTCGTTAGTCTACATGTTTCTAAACATATAAATGTGAATTTAATCAATTCCTTTCATAGTTTTATAATTCTCTGGCAGTTCCTTATGATAGAGTTTATAAAACAGTCCTGTGTAAACTGCTGGAAGTTCTTCCACAGTCAGGTCAATTTTGTCAAACCCTTCTCTGTACCCATACAGCAGCAGCCTAGCAACTCTGCTGGTGATGGGAGTTGTATTTTCAGTCTTCGCCAGGTCATTGAGATCCATCCACTCACATCTTAAGCATTCTTCCTGGCAAAAATTTATGGTGAATGAATATGGCTTTAGGCGGCAGATGATATACATATCTGACTTCCCAAAAGCTCCAGGATTTGTGTGCTGTTGCCGAATACTCAGGACGGACCTGAATTCTGATTTTATACCAGTCTCTTCAAAAACTTCTCGAACCGCTGTGTCTCCAATATCTTCTTCAGGCTCTGACAGGCCTCCTGGAAACTTCCACATATTTTTCAATTTATTTCGATCTTGTACAACCAGTATTTTTCTAGTAC"),
    ('gi|67678568|gb|DR109384.1|DR109384',
     "RTS1_1_D03.g1_A029 Roots minus sulfur Pinus taeda cDNA clone RTS1_1_D03_A029 5', mRNA sequence",
     strand.plus,
     "ACCGAAGAGGAGACT"),
    ('gi|1111111',
     "foobar 3'",
     strand.minus,
     "ACCGAAGAGGAGACT"),
    ('lcl|AI370586',
     "ta40b12.x1 Soares_total_fetus_Nb2HF8_9w Homo sapiens cDNA clone IMAGE:2046527 3' similar to gb:J04513 HEPARIN-BINDING GROWTH FACTOR PRECURSOR 2 (HUMAN);, mRNA sequence. #CRS4 tissue_type-not reported- organ-not reported-#",
     strand.minus,
     "TTTTTTTCTTTTTCACGCATTTGCTTTATTCGAAAAGAGGCTTTTAAAATGTGCATGTTTAGAAACAAAATTTCTTCATGGAAATCATATACATTAGAAAATCACAGTCAGATGTTTAATCAATCCAAAATGTCCACTATTTCTTATGTCATTCGTTAGTCTACATGTTTCTAAACATATAAATGTGAATTTAATCAATTCCTTTCATAGTTTTATAATTCTCTGGCAGTTCCTTATGATAGAGTTTATAAAACAGTCCTGTGTAAACTGCTGGAAGTTCTTCCACAGTCAGGTCAATTTTGTCAAACCTTTTTTGGACCCATACAGCAGCAGCCTAACACTTTCTGGGATGGGGAGTTGTATTTTAAGTCTCGCCAGGGCATTGAGATCATCCACTCACATCTAAAGCATTCTTCCTGCCAAATATATGGGGAATGAATATGGCTTTAGCCGGAGATGAATTACATTTTTGACTTCCCAAA"),
    ('lcl|BM970451',
     "UI-CF-EC1-abs-k-20-0-UI.s1 UI-CF-EC1 Homo sapiens cDNA clone UI-CF-EC1-abs-k-20-0-UI 3', mRNA sequence. #CRS4 tissue_type-Lung- organ-Lung-#",
     strand.minus,
     "TTTTTTTTTTTTCTTTTTCACGCATTTGCTTTATTCGAAAAGAGGCTTTTAAAATGTGCATGTTTAGAAACAAAATTTCTTCATGGAAATCATATACATTAGAAAATCACAGTCAGATGTTTAATCAATCCAAAATGTCCACTATTTCTTATGTCATTCGTTAGTCTACATGTTTCTAAACATATAAATGTGAATTTAATCAATTCCTTTCATAGTTTTATAATTCTCTGGCAGTTCCTTATGATAGAGTTTATAAAACAGTCCTGTGTAAACTGCTGGAAGTTCTTCCACAGTCAGGTCAATTTTGTCAAACCCTTCTCTGTACCCATACAGCAGCAGCCTAGCAACTCTGCTGGTGATGGGAGTTGTATTTTCAGTCTTCGCCAGGTCATTGAGATCCATCCACTCACATCTTAAGCATTCTTCCTGGCAAAAATTTATGGTGAATGAATATGGCTTTAGGCGGCAGATGATATACATATCTGACTTCCCAAAAGCTCCAGGATTTGTGTGCTGTTGCCGAATACTCAGGACGGACCTGAATTCTGATTTTATACCAGTCTCTTCAAAAACTTCTCGAACCGCTGTGTCTCCAATATCTTCTTCAGGCTCTGACAGGCCTCCTGGAAACTTCCACATATTTTTCAATTTATTTCGATCTTGTACAACCAGTATTTTTCTAGTAC"),
    ('lcl|AA193558',
     "zr41f02.r1 Soares_NhHMPu_S1 Homo sapiens cDNA clone IMAGE:665979 5' similar to gb:J04513 HEPARIN-BINDING GROWTH FACTOR PRECURSOR 2 (HUMAN);, mRNA sequence. #CRS4 tissue_type-Pooled human melanocyte, fetal heart, and pregnant uterus- organ-mixed (see below)-#",
     strand.plus,
     'GCAGGAGCTGTATTTGATGAAAGTACTAGAAAAATACTGGTTGTACAAGATCGAAATAAATTGAAAAATATGTGGAAGTTTCCAGGAGGCCTGTCAGAGCCTGAAGAAGATATTGGAGACACAGCGGTTCGAGAAGTTTTTGAAGAGACTGGTATAAAATCAGAATTCAGGTCCGTCCTGAGTATTCAGCAACAGCACACAAATCCTGGAGCTTTTGGGAAGTCAGATATGTATATCATCTGCCGCCTAAAGCCATATTCATTCACCATAAATTTTTGCCAGGAAGAATGCTTAAGATGTGAGTGGATGGATCTCAATGACCTGGCGAAGACTGAANATTACAACTCCCATCACCAGCAGAGTTGCTAAGGCTGATGCNGGTATGGGTACAGAGAAAGGNTTTGACCAAAATTGACCTGGCCTGGGGGAAGAACTTCCCAGCAGNTTT')
]


class streams_tc(unittest.TestCase):

    def setUp(self):
        self.sseqs = map (lambda s : (">%s %s\n%s" % (s[0], s[1], s[3]),
                                      s[0], s[1], s[2], s[3]), seq_tuples)
        self.fname = 'ciccio.fa'
        fout = open(self.fname, 'wb')
        for s in self.sseqs:
            fout.write(s[0] + '\n')
        fout.close()

    def test_FS_IO(self):
        for bufsize in 0, 50, 500:
            FS_OBJ = FS_IO(self.fname, bufsize)
            lines = [l for l in FS_OBJ]
            for l1, l2 in itertools.izip(
                lines,
                itertools.chain(*[s[0].splitlines() for s in self.sseqs])
                ):
                self.assertEqual(l1, l2+"\n")

    def test_FS_IO_no_newline_at_EOF(self):
        fname = 'ciccio.no_newline.fa'
        fout = open(fname, 'wb')
        for s in self.sseqs[:-1]:
            fout.write(s[0] + '\n')
        fout.write(self.sseqs[-1][0])
        fout.close()
        for bufsize in 0, 50, 500:
            FS_OBJ = FS_IO(fname, bufsize)
            lines = [l for l in FS_OBJ]
            ok_lines = list(itertools.chain(*[s[0].splitlines()
                                              for s in self.sseqs]))
            for l1, l2 in itertools.izip(lines[:-1], ok_lines[:-1]):
                self.assertEqual(l1, l2+"\n")
            self.assertEqual(lines[-1], ok_lines[-1])

    def test_fasta_stream_from_file_unbuffered(self):
        FS_OBJ = FS_IO(self.fname, 0)
        fsff = fasta_stream_from_stream(FS_OBJ)
        counter = 0
        for s in fsff:
            self.assertEqual(s[-1], '\n')
            self.assertEqual(s[:-1], self.sseqs[counter][0])
            counter += 1
    def test_fasta_stream_from_file_buffered(self):
        FS_OBJ = FS_IO(self.fname, 1000)
        fsff = fasta_stream_from_stream(FS_OBJ)
        counter = 0
        for s in fsff:
            self.assertEqual(s[-1], '\n')
            self.assertEqual(s[:-1], self.sseqs[counter][0])
            counter += 1
    def test_blast_seq_stream(self):
        FS_OBJ = FS_IO(self.fname, 1000)
        sf     = seq_factory_from_fasta(strand.unknown)
        dbEST  = fasta_stream_from_stream(FS_OBJ)
        seqEST = blast_seq_stream(sf, dbEST)
        counter = 0
        for s in seqEST:
            self.assertEqual(s.get_sequence(), self.sseqs[counter][4])
            self.assertEqual(s.id, self.sseqs[counter][1])
            self.assertEqual(s.title, self.sseqs[counter][2])
            counter += 1


def suite():
    suite = unittest.TestSuite()
    suite.addTest(streams_tc('test_fasta_stream_from_file_unbuffered'))
    suite.addTest(streams_tc('test_fasta_stream_from_file_buffered'))
    suite.addTest(streams_tc('test_blast_seq_stream'))
    suite.addTest(streams_tc('test_FS_IO'))
    suite.addTest(streams_tc('test_FS_IO_no_newline_at_EOF'))
    return suite


if __name__ == '__main__':
  runner = unittest.TextTestRunner(verbosity=2)
  runner.run((suite()))
