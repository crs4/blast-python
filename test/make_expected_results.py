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
Generate expected blast results for all query/subject combinations of
sequences found in the given subdirectory. Outputs a .py file from
which a results dict can be imported.
"""

# reference results for multi-query/subject blasts are computed with
# blast instead of bl2seq, because the latter is not able to perform
# multi-query/subject blasts and results are slightly different from
# those obtained blasting each sequence pair separately (namely we
# observed that, in some cases, when doing a multi-query blast
# additional results appear which are identical to other results but
# with query and subject roles reversed).

# both standard and tabular blast output files are parsed to get
# reference results, because the former has a bug that makes it miss
# the coefficient in the scientific notation of evalues and the latter
# does a weird clustering of bit scores which divides them into
# macro-groups. Therefore we take everything from standard output and
# then replace evalues with those obtained from tabular output.

import sys, os, shutil, optparse, pprint, re
from itertools import izip
from operator import itemgetter


class BlastResultsParser(object):
    """
    Parser for bl2seq output.
    """

##     FIELDS = ["percent_identity",
##               "alignment_length",
##               "mismatches",
##               "gap_openings",
##               "query_start",
##               "query_end",
##               "subject_start",
##               "subject_end",
##               "e_value",
##               "bit_score"]

    FIELDS = ["query_start",
              "query_end",
              "subject_start",
              "subject_end",
              "e_value",
              "bit_score"]

    RE_FLAGS = [re.IGNORECASE]

    def __init__(self):
        self.reset()

    def reset(self):
        self.current_query_id = None
        self.current_subject_id = None
        self.current_match = {}
        self.std_results = {}
        self.tab_results = {}

    def parse_dir(self, resdir, outtype="std"):
        parse = getattr(self, "parse_%s_file" % outtype, self.parse_std_file)
        resfiles = [os.path.join(resdir,f) for f in os.listdir(resdir)]
        for resfile in resfiles:
            resfp = open(resfile)
            parse(resfp)
            resfp.close()

    def store_match(self):
        for field in self.FIELDS:
            assert self.current_match[field] is not None
        res = self.std_results.setdefault(
            (self.current_query_id, self.current_subject_id), []
            )
        res.append(self.current_match)

    def parse_std_file(self, resfp):
        """
        Parse blast results from standard (verbose) output file.
        """
        for line in resfp:
            line = line.strip()
            if not line:
                continue
            m = re.search(r"query\s*=\s*(\S+)", line, *self.RE_FLAGS)
            if m is not None:
                if self.current_match:
                    self.store_match()
                    self.current_match = {}
                self.current_query_id = m.groups()[0]
            m = re.search(r">(\S+)", line, *self.RE_FLAGS)
            if m is not None:
                if self.current_match:
                    self.store_match()
                    self.current_match = {}
                self.current_subject_id = m.groups()[0]
            m = re.search(r"score\s*=\s*(\S+).+expect[^=]*=\s*(\S+)",
                          line, *self.RE_FLAGS)
            if m is not None:
                if self.current_match:
                    self.store_match()
                    self.current_match = {}
                self.current_match["bit_score"] = float(m.groups()[0])
                try:  # quickfix... Take evalues from tabular output instead
                    self.current_match["e_value"] = float(m.groups()[1])
                except ValueError:
                    self.current_match["e_value"] = float("1"+m.groups()[1])
            m = re.search(r"query:\s*(\d+)\D+(\d+)", line, *self.RE_FLAGS)
            if m is not None:
                self.current_match.setdefault("query_start",
                                              float(m.groups()[0]))
                self.current_match["query_end"] = float(m.groups()[1])
            m = re.search(r"sbjct:\s*(\d+)\D+(\d+)", line, *self.RE_FLAGS)
            if m is not None:
                self.current_match.setdefault("subject_start",
                                              float(m.groups()[0]))
                self.current_match["subject_end"] = float(m.groups()[1])
        if self.current_match:
            self.store_match()
            self.current_match = {}

    def parse_tab_file(self, resfp):
        """
        Parse blast results from tabular output file.
        """
        for line in resfp:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            line = line.split("\t")
            query_id, subject_id = line[:2]
            res = self.tab_results.setdefault((query_id, subject_id), [])
            res.append(dict(zip(self.FIELDS, [float(v) for v in line[-6:]])))


def get_results(seq_files, progname, multi="none"):
    p = BlastResultsParser()
    blast_args = "-i %s -j %s -p %s -o %s"
    tmpfile = None
    if multi == "none":
        tool = "bl2seq"
    else:
        tool = "blast"
        all_sequences = []
        seq_ids = []
        for i, fn in enumerate(seq_files):
            fp = open(fn)
            seq = fp.read().strip()
            seq_ids.append(re.search(r">(\S+)", seq).groups()[0])
            all_sequences.append(seq)
            fp.close()
        all_sequences = "\n".join(all_sequences) + "\n"
        tmpfile = "all_sequences.fa"
        tmpfp = open(tmpfile, "w")
        tmpfp.write(all_sequences)
        tmpfp.close()

    blast_cl = tool + " " +  blast_args
    
    for outtype in "std", "tab":
        tmpdir = "%s_output_%s" % (tool, outtype)
        if outtype == "tab":
            if multi == "none":
                blast_cl += " -D 1"
            else:
                blast_cl += " -m 9"
        try:
            shutil.rmtree(tmpdir)
        except OSError:
            pass
        os.mkdir(tmpdir)
        if multi == "none":
            for fn1 in seq_files:
                for fn2 in seq_files:
                    outfn = os.path.join(tmpdir, "%s_vs_%s.out" %
                                         (os.path.basename(fn1),
                                          os.path.basename(fn2)))
                    retval = os.system(blast_cl % (fn1, fn2, progname, outfn))
                    if retval != 0:
                        sys.exit("This script requires NCBI blast")
        elif multi == "query":
            for i, fn in enumerate(seq_files):
                outfn = os.path.join(tmpdir, "ALL_vs_%s.out" %
                                     os.path.basename(fn))
                retval = os.system(blast_cl % (tmpfile, fn, progname, outfn))
                if retval != 0:
                    sys.exit("This script requires NCBI blast")
                outfp = open(outfn)
                content = outfp.read()
                outfp.close()
                content = content.replace("1_0", seq_ids[i])
                outfp = open(outfn, "w")
                outfp.write(content)
                outfp.close()
        elif multi == "subject":
            for fn in seq_files:
                outfn = os.path.join(tmpdir, "%s_vs_ALL.out" %
                                     os.path.basename(fn))
                retval = os.system(blast_cl % (fn, tmpfile, progname, outfn))
                if retval != 0:
                    sys.exit("This script requires NCBI blast")
                outfp = open(outfn)
                content = outfp.read()
                outfp.close()
                for i in xrange(len(seq_ids)):
                    ordinal = "%d_0" % (i+1)
                    content = content.replace(ordinal, seq_ids[i])
                outfp = open(outfn, "w")
                outfp.write(content)
                outfp.close()
        elif multi == "both":
            outfn = os.path.join(tmpdir, "ALL_vs_ALL.out") 
            retval = os.system(blast_cl % (tmpfile, tmpfile, progname, outfn))
            if retval != 0:
                sys.exit("This script requires NCBI blast")
            outfp = open(outfn)
            content = outfp.read()
            outfp.close()
            for i in xrange(len(seq_ids)):
                ordinal = "%d_0" % (i+1)
                content = content.replace(ordinal, seq_ids[i]) 
            outfp = open(outfn, "w")
            outfp.write(content)
            outfp.close()
        p.parse_dir(tmpdir, outtype)
        shutil.rmtree(tmpdir)
        
    if tmpfile is not None:
        os.remove(tmpfile)

    keys = p.std_results.keys()
    assert keys == p.tab_results.keys()
    for k in keys:
        try:
            assert len(p.std_results[k]) == len(p.tab_results[k])
        except AssertionError:
            print k, len(p.std_results[k]), len(p.tab_results[k])
            raise
    for res_dict in p.std_results, p.tab_results:
        for v in res_dict.itervalues():
            v.sort(key=itemgetter("subject_end"))
            v.sort(key=itemgetter("subject_start"))
            v.sort(key=itemgetter("query_end"))
            v.sort(key=itemgetter("query_start"))
    for k in keys:
        for std_match, tab_match in izip(p.std_results[k], p.tab_results[k]):
            std_match["e_value"] = tab_match["e_value"]
    return p.std_results


def make_parser():
    """
    Build command line parser.
    """
    usage = "%prog [OPTIONS]"
    parser = optparse.OptionParser(usage)
    parser.set_description(__doc__.lstrip())
    parser.add_option("-d", type="string", dest="seqdir", metavar="STRING",
                      help="directory with sequences in fasta format")
    parser.add_option("-p", type="choice", dest="program", metavar="STRING",
                      help="blast program to use",
                      choices = ["blastn", "tblastx"])
    parser.add_option("-m", type="choice", dest="multi", metavar="STRING",
                      help="multi sequence blast",
                      choices = ["none", "query", "subject", "both"])

    parser.set_defaults(
        seqdir="sequences",
        program="blastn",
        multi="none"
        )

    for s in "-p", "-m":
        opt = parser.get_option(s)
        opt.help += " %r" % opt.choices
    
    return parser


def main(argv):

    parser = make_parser()
    opt, args = parser.parse_args(argv)

    if opt.multi == "none":
        tool = "bl2seq"
    else:
        tool = "blast"
    
    seq_files = [os.path.join(opt.seqdir,f) for f in os.listdir(opt.seqdir)
                 if f.endswith(".fa")]

    results = get_results(seq_files, opt.program, opt.multi)

    header = """# Expected results for test_blaster.py. Computed with NCBI
# %s, using the %s program, for all query/subject combinations of the
# following queries:\n%s
""" % (tool, opt.program,
       "\n".join(["# "+os.path.basename(fn) for fn in seq_files]))

    if opt.multi == "none":
        tag = ""
    else:
        tag = "_multi_%s" % opt.multi
    res_fn = "%s_results%s.py" % (opt.program, tag)
    dict_n = "%s_RESULTS%s = " % (opt.program.upper(), tag.upper())
    res_fp = open(res_fn, "w")
    pp = pprint.PrettyPrinter(stream=res_fp)
    res_fp.write(header + "\n")
    res_fp.write(dict_n)
    pp.pprint(results)
    res_fp.close()

    total_comb = len(seq_files)**2
    print len(results), "out of", total_comb, "seq combinations produced hits"
    print "Results saved to", res_fn


if __name__ == "__main__":
    main(sys.argv)
