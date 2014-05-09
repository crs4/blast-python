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

### blastn tests fail when query == subject (condition named
### 'diagonal' because corresponding results lie in the diagonal of
### the N by N matrix of all sequences vs themselves). This is
### probably due to some weird bug in blast / bl2seq: instead of
### reporting a single match corresponding to the whole sequence with
### evalue = 0, they output several partial matches with lower
### evalues.


import unittest, os, copy
from itertools import izip
from operator import attrgetter, itemgetter

import ncbi_toolkit
from BlastPython import blaster


# possible values: both_rev, plus, minus, unknown, other, both
# both_rev and other crash the program (exits with "Aborted.")
STRAND = ncbi_toolkit.strand.both  # NCBI bl2seq default
SEQDIR="sequences"


try:
    from blastn_results import BLASTN_RESULTS
    from tblastx_results import TBLASTX_RESULTS
    from blastn_results_multi_query import BLASTN_RESULTS_MULTI_QUERY
    from tblastx_results_multi_query import TBLASTX_RESULTS_MULTI_QUERY    
    from blastn_results_multi_subject import BLASTN_RESULTS_MULTI_SUBJECT
    from tblastx_results_multi_subject import TBLASTX_RESULTS_MULTI_SUBJECT
    from blastn_results_multi_both import BLASTN_RESULTS_MULTI_BOTH
    from tblastx_results_multi_both import TBLASTX_RESULTS_MULTI_BOTH
except ImportError:
    # try to generate expected results - requires NCBI toolkit
    import make_expected_results
    for prog in "blastn", "tblastx":
        make_expected_results.main(["make_expected_results.py",
                                    "-d", SEQDIR, "-p", prog])
        make_expected_results.main(["make_expected_results.py",
                                    "-d", SEQDIR, "-p", prog, "-m", "query"])
        make_expected_results.main(["make_expected_results.py",
                                    "-d", SEQDIR, "-p", prog, "-m", "subject"])
        make_expected_results.main(["make_expected_results.py",
                                    "-d", SEQDIR, "-p", prog, "-m", "both"])
        
    from blastn_results import BLASTN_RESULTS
    from tblastx_results import TBLASTX_RESULTS
    from blastn_results_multi_query import BLASTN_RESULTS_MULTI_QUERY
    from tblastx_results_multi_query import TBLASTX_RESULTS_MULTI_QUERY
    from blastn_results_multi_subject import BLASTN_RESULTS_MULTI_SUBJECT
    from tblastx_results_multi_subject import TBLASTX_RESULTS_MULTI_SUBJECT
    from blastn_results_multi_both import BLASTN_RESULTS_MULTI_BOTH
    from tblastx_results_multi_both import TBLASTX_RESULTS_MULTI_BOTH
    
for exp_res_dict in (BLASTN_RESULTS, TBLASTX_RESULTS,
                     BLASTN_RESULTS_MULTI_QUERY, TBLASTX_RESULTS_MULTI_QUERY,
                     BLASTN_RESULTS_MULTI_SUBJECT,
                     TBLASTX_RESULTS_MULTI_SUBJECT,
                     BLASTN_RESULTS_MULTI_BOTH, TBLASTX_RESULTS_MULTI_BOTH):
    for v in exp_res_dict.itervalues():
        v.sort(key=itemgetter("subject_end"))
        v.sort(key=itemgetter("subject_start"))
        v.sort(key=itemgetter("query_end"))
        v.sort(key=itemgetter("query_start"))


def get_sequences(filenames):
    sequences = []
    for fn in filenames:
        fp = open(fn)
        sequences.append(fp.read())
        fp.close()
    return sequences


def get_modified_value(v):
    t = type(v)
    if t == type(1):
        if v > 0:
            nv = v - 1
        else:
            nv = 1
        return nv
    elif t == type(0.5):
        return (v + 1.33) * 0.5
    else:
        return v


def get_program(progname):
    if progname == "blastn":
        return ncbi_toolkit.EProgram.eBlastn
    elif progname == "tblastx":
        return ncbi_toolkit.EProgram.eTblastx
    else:
        return ncbi_toolkit.EProgram.eBlastn


# hsp.query/subject tuples: (<Translation frame>, <Start of hsp>, <End
# of HSP>, <Where the gapped extension started>)
def get_match_endpoints(q_tuple, s_tuple, q_len, s_len, progname="blastn"):
    q_frame, q_start, q_end, q_gap_start = q_tuple
    s_frame, s_start, s_end, s_gap_start = s_tuple
    if progname == "blastn":
        if q_frame == -1:
            q_start, q_end = q_len - q_end + 1, q_len - q_start
            s_start, s_end = s_end, s_start + 1
        else:
            q_start += 1
            s_start += 1
    elif progname == "tblastx":
        if q_frame >= 0:
            q_start, q_end = 3*q_start + q_frame, 3*q_end + q_frame - 1
        else:
            q_start, q_end = q_len - 3*q_start + q_frame + 1, \
                             q_len - 3*q_end + q_frame + 2
        if s_frame >= 0:
            s_start, s_end = 3*s_start + s_frame, 3*s_end + s_frame - 1
        else:
            s_start, s_end = s_len - 3*s_start + s_frame + 1, \
                             s_len - 3*s_end + s_frame + 2
    else:
        raise ValueError("%s program not supported" % progname)

    return q_start, q_end, s_start, s_end


factory = ncbi_toolkit.blast_sseq_loc_from_fasta()
seq_files = [os.path.join(SEQDIR,f) for f in os.listdir(SEQDIR)
             if f.endswith(".fa")]
sequences = [factory.make(s, STRAND, 0, 0, False)
             for s in get_sequences(seq_files)]


class set_options_tc(unittest.TestCase):
    
    def setUp(self):
        self.sequences = sequences

    def test_set_options(self):
        for progname in "blastn", "tblastx":
            for seq in self.sequences:
                prog = get_program(progname)
                self.blaster = blaster(seq, Program=prog)
                og = self.blaster.get_options()
                os = self.blaster.set_options()
                for k in og.keys():
                    v = og[k]
                    nv = get_modified_value(v)
                    os[k] = nv
                    self.assertEqual(nv, og[k])


class base_blast_tc(object):
    """
    Base class for blast test cases. Derived test classes must define
    the following attributes in setUp: self.sequences,
    self.exp_res_dict, self.blast_options. setUp must also call
    compute_results AFTER defining these attributes.
    """
    def compute_results(self):
        """
        Compute results dict.
        """
        self.res_dict = {}
        for s1 in self.sequences:
            b = blaster(s1, **self.blast_options)
            for s2 in self.sequences:
                subject, r = b.blast(s2)
                self.assertEqual(subject, s2)
                hit_list = r[0]
                if hit_list is None:
                    continue
                hsp_list = hit_list[0]
                results = []
                for hsp in hsp_list:
                    res = {}
                    res['bit_score'] = hsp.bit_score
                    res['e_value'] = hsp.evalue
                    (res['query_start'],
                     res['query_end'],
                     res['subject_start'],
                     res['subject_end']) = get_match_endpoints(
                        hsp.query, hsp.subject, s1.length, s2.length,
                        self.progname
                        )
                    results.append(res)
                results.sort(key=itemgetter("subject_end"))
                results.sort(key=itemgetter("subject_start"))
                results.sort(key=itemgetter("query_end"))
                results.sort(key=itemgetter("query_start"))
            self.res_dict[(s1.id, s2.id)] = results

    def check_result(self, r, exp_r):
        self.assertAlmostEqual(r['bit_score']/exp_r['bit_score'], 1, 1)
        self.assertAlmostEqual(r['e_value']/exp_r['e_value'], 1, 0)
        for k in 'query_start', 'query_end', 'subject_start', 'subject_end':
            self.assertEqual(int(r[k]), int(exp_r[k]))

    def test_n_hits(self, diagonal=False):
        for s1 in self.sequences:
            for s2 in self.sequences:
                if (diagonal and s1 != s2) or (not diagonal and s1 == s2):
                    continue
                try:
                    results = self.res_dict[(s1.id, s2.id)]
                    exp_results = self.exp_res_dict[(s1.id, s2.id)]
                except KeyError:
                    continue  # no match for this combination
                self.assertEqual(len(results), len(exp_results))

    def test_all_hits(self, diagonal=False):
        for s1 in self.sequences:
            for s2 in self.sequences:
                if (diagonal and s1 != s2) or (not diagonal and s1 == s2):
                    continue
                try:
                    results = self.res_dict[(s1.id, s2.id)]
                    exp_results = self.exp_res_dict[(s1.id, s2.id)]
                except KeyError:
                    continue  # no match for this combination
                for r, exp_r in izip(results, exp_results):
                    self.check_result(r, exp_r)


class blastn_tc(base_blast_tc, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "blastn"
        self.exp_res_dict = BLASTN_RESULTS
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eBlastn,
                              'MatchReward': 1}
        self.compute_results()

    def test_blastn_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_blastn_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_blastn_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_blastn_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


class tblastx_tc(base_blast_tc, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "tblastx"
        self.exp_res_dict = TBLASTX_RESULTS
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eTblastx}
        self.compute_results()

    def test_tblastx_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_tblastx_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_tblastx_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_tblastx_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


class blast_tc_multi_query(base_blast_tc):
    """
    Base class for multi-query blast test cases. Derived test classes
    must define the following attributes in setUp: self.sequences,
    self.exp_res_dict, self.blast_options. setUp must also call
    compute_results AFTER defining these attributes.
    """
    def compute_results(self):
        """
        Compute results dict.
        """
        self.res_dict = {}
        b = blaster(self.sequences, **self.blast_options)
        for s2 in self.sequences:
            subject, r = b.blast(s2)
            self.assertEqual(len(r), len(self.sequences))            
            for hit_list in r:
                if hit_list is None:
                    continue
                hsp_list = hit_list[0]
                query = self.sequences[hsp_list.query_index]
                #print "%s vs %s" % (query.id, s2.id)
                results = []
                for hsp in hsp_list:
                    res = {}
                    res['bit_score'] = hsp.bit_score
                    res['e_value'] = hsp.evalue
                    (res['query_start'],
                     res['query_end'],
                     res['subject_start'],
                     res['subject_end']) = get_match_endpoints(
                        hsp.query, hsp.subject, query.length, s2.length,
                        self.progname
                        )
                    results.append(res)
                results.sort(key=itemgetter("subject_end"))
                results.sort(key=itemgetter("subject_start"))
                results.sort(key=itemgetter("query_end"))
                results.sort(key=itemgetter("query_start"))
            self.res_dict[(query.id, s2.id)] = results


class blast_tc_multi_subject(base_blast_tc):
    """
    Base class for multi-subject blast test cases. Derived test classes
    must define the following attributes in setUp: self.sequences,
    self.exp_res_dict, self.blast_options. setUp must also call
    compute_results AFTER defining these attributes.
    """
    def compute_results(self):
        """
        Compute results dict.
        """
        self.res_dict = {}
        for s1 in self.sequences:
            b = blaster(s1, **self.blast_options)
            subjects, r = b.blast(self.sequences)
            hit_list = r[0]
            if hit_list is None:
                continue            
            self.assertEqual(subjects, self.sequences)
            for hsp_list in hit_list:
                subject = self.sequences[hsp_list.ordinal_id_subject_sequence]
                #print "%s vs %s" % (s1.id, subject.id)
                results = []
                for hsp in hsp_list:
                    res = {}
                    res['bit_score'] = hsp.bit_score
                    res['e_value'] = hsp.evalue
                    (res['query_start'],
                     res['query_end'],
                     res['subject_start'],
                     res['subject_end']) = get_match_endpoints(
                        hsp.query, hsp.subject, s1.length, subject.length,
                        self.progname
                        )
                    results.append(res)
                results.sort(key=itemgetter("subject_end"))
                results.sort(key=itemgetter("subject_start"))
                results.sort(key=itemgetter("query_end"))
                results.sort(key=itemgetter("query_start"))
            self.res_dict[(s1.id, subject.id)] = results


class blast_tc_multi_both(base_blast_tc):
    """
    Base class for multi-query, multi-subject blast test
    cases. Derived test classes must define the following attributes
    in setUp: self.sequences, self.exp_res_dict,
    self.blast_options. setUp must also call compute_results AFTER
    defining these attributes.
    """
    def compute_results(self):
        """
        Compute results dict.
        """
        self.res_dict = {}
        b = blaster(self.sequences, **self.blast_options)
        subjects, r = b.blast(self.sequences)
        self.assertEqual(len(r), len(self.sequences))
        self.assertEqual(subjects, self.sequences)
        for hit_list in r:
            if hit_list is None:
                continue
            for hsp_list in hit_list:
                query = self.sequences[hsp_list.query_index]
                subject = self.sequences[hsp_list.ordinal_id_subject_sequence]
                #print "%s vs %s" % (query.id, subject.id) 
                results = []
                for hsp in hsp_list:
                    res = {}
                    res['bit_score'] = hsp.bit_score
                    res['e_value'] = hsp.evalue
                    (res['query_start'],
                     res['query_end'],
                     res['subject_start'],
                     res['subject_end']) = get_match_endpoints(
                        hsp.query, hsp.subject, query.length, subject.length,
                        self.progname
                        )
                    results.append(res)
                results.sort(key=itemgetter("subject_end"))
                results.sort(key=itemgetter("subject_start"))
                results.sort(key=itemgetter("query_end"))
                results.sort(key=itemgetter("query_start"))
            self.res_dict[(query.id, subject.id)] = results


class blastn_tc_multi_query(blast_tc_multi_query, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "blastn"
        self.exp_res_dict = BLASTN_RESULTS_MULTI_QUERY
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eBlastn,
                              'MatchReward': 1}
        self.compute_results()

    def test_blastn_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_blastn_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_blastn_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_blastn_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


class tblastx_tc_multi_query(blast_tc_multi_query, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "tblastx"
        self.exp_res_dict = TBLASTX_RESULTS_MULTI_QUERY
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eTblastx,
                              'MatchReward': 1}
        self.compute_results()

    def test_tblastx_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_tblastx_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_tblastx_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_tblastx_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


class blastn_tc_multi_subject(blast_tc_multi_subject, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "blastn"
        self.exp_res_dict = BLASTN_RESULTS_MULTI_SUBJECT
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eBlastn,
                              'MatchReward': 1}
        self.compute_results()

    def test_blastn_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_blastn_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_blastn_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_blastn_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


class tblastx_tc_multi_subject(blast_tc_multi_subject, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "tblastx"
        self.exp_res_dict = TBLASTX_RESULTS_MULTI_SUBJECT
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eTblastx,
                              'MatchReward': 1}
        self.compute_results()

    def test_tblastx_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_tblastx_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_tblastx_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_tblastx_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


class blastn_tc_multi_both(blast_tc_multi_both, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "blastn"
        self.exp_res_dict = BLASTN_RESULTS_MULTI_BOTH
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eBlastn,
                              'MatchReward': 1}
        self.compute_results()

    def test_blastn_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_blastn_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_blastn_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_blastn_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


class tblastx_tc_multi_both(blast_tc_multi_both, unittest.TestCase):

    def setUp(self):
        self.sequences = sequences
        self.progname = "tblastx"
        self.exp_res_dict = TBLASTX_RESULTS_MULTI_BOTH
        self.blast_options = {'Program': ncbi_toolkit.EProgram.eTblastx,
                              'MatchReward': 1}
        self.compute_results()

    def test_tblastx_no_diagonal_n_hits(self):
        self.test_n_hits(diagonal=False)

    def test_tblastx_no_diagonal_all_hits(self):
        self.test_all_hits(diagonal=False)

    def test_tblastx_diagonal_n_hits(self):
        self.test_n_hits(diagonal=True)

    def test_tblastx_diagonal_all_hits(self):
        self.test_all_hits(diagonal=True)


def suite():
    suite = unittest.TestSuite()
    
    suite.addTest(set_options_tc('test_set_options'))
    
    suite.addTest(blastn_tc('test_blastn_no_diagonal_n_hits'))
    suite.addTest(blastn_tc('test_blastn_no_diagonal_all_hits'))
##     suite.addTest(blastn_tc('test_blastn_diagonal_n_hits'))
##     suite.addTest(blastn_tc('test_blastn_diagonal_all_hits'))
    suite.addTest(tblastx_tc('test_tblastx_no_diagonal_n_hits'))
    suite.addTest(tblastx_tc('test_tblastx_no_diagonal_all_hits'))
    suite.addTest(tblastx_tc('test_tblastx_diagonal_n_hits'))
    suite.addTest(tblastx_tc('test_tblastx_diagonal_all_hits'))

    suite.addTest(blastn_tc_multi_query('test_blastn_no_diagonal_n_hits'))
    suite.addTest(blastn_tc_multi_query('test_blastn_no_diagonal_all_hits'))
##     suite.addTest(blastn_tc_multi_query('test_blastn_diagonal_n_hits'))
##     suite.addTest(blastn_tc_multi_query('test_blastn_diagonal_all_hits'))
    suite.addTest(tblastx_tc_multi_query('test_tblastx_no_diagonal_n_hits'))
    suite.addTest(tblastx_tc_multi_query('test_tblastx_no_diagonal_all_hits'))
    suite.addTest(tblastx_tc_multi_query('test_tblastx_diagonal_n_hits'))
    suite.addTest(tblastx_tc_multi_query('test_tblastx_diagonal_all_hits'))

    suite.addTest(blastn_tc_multi_subject('test_blastn_no_diagonal_n_hits'))
    suite.addTest(blastn_tc_multi_subject('test_blastn_no_diagonal_all_hits'))
##     suite.addTest(blastn_tc_multi_subject('test_blastn_diagonal_n_hits'))
##     suite.addTest(blastn_tc_multi_subject('test_blastn_diagonal_all_hits'))
    suite.addTest(tblastx_tc_multi_subject('test_tblastx_no_diagonal_n_hits'))
    suite.addTest(tblastx_tc_multi_subject
                  ('test_tblastx_no_diagonal_all_hits'))
    suite.addTest(tblastx_tc_multi_subject('test_tblastx_diagonal_n_hits'))
    suite.addTest(tblastx_tc_multi_subject('test_tblastx_diagonal_all_hits'))

    suite.addTest(blastn_tc_multi_both('test_blastn_no_diagonal_n_hits'))
    suite.addTest(blastn_tc_multi_both('test_blastn_no_diagonal_all_hits'))
##     suite.addTest(blastn_tc_multi_both('test_blastn_diagonal_n_hits'))
##     suite.addTest(blastn_tc_multi_both('test_blastn_diagonal_all_hits'))
    suite.addTest(tblastx_tc_multi_both('test_tblastx_no_diagonal_n_hits'))
    suite.addTest(tblastx_tc_multi_both('test_tblastx_no_diagonal_all_hits'))
    suite.addTest(tblastx_tc_multi_both('test_tblastx_diagonal_n_hits'))
    suite.addTest(tblastx_tc_multi_both('test_tblastx_diagonal_all_hits'))
        
    return suite


if __name__ == '__main__':
  runner = unittest.TextTestRunner(verbosity=2)
  runner.run((suite()))
