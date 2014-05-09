#!/usr/bin/env python
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
A (minimal) Python version of bl2seq.
"""

import sys, optparse
from ncbi_toolkit import *

# possible values: both_rev, plus, minus, unknown, other, both
# both_rev, other: programs exit with "Aborted."
STRAND = strand.both  # NCBI bl2seq default

# 2006 default values
TBLASTX_OPTS_06 = {
    "DustFilteringLevel": 20,  # (-1)
    "DustFilteringWindow": 64,  # (-1)
    "RepeatFilteringDB": "humrep",  # (None)
    "DustFilteringLinker": 1  # (-1)
    }


def make_parser():
    """
    Build command line parser.

    @rtype: L{OptionParser<optparse.OptionParser>}
    @return: command line parser with options and defaults set
    """
    usage = "%prog QUERY_FILE TARGET_FILE [OPTIONS]"
    parser = optparse.OptionParser(usage)
    parser.set_description(__doc__.lstrip())
    parser.add_option("-l", type="int", dest="DbLength", metavar="INT",
                      help="database length")
    parser.add_option("-t", type="float", dest="EvalueThreshold",
                      metavar="FLOAT", help="evalue threshold for results")
    parser.add_option("-p", type="string", dest="prog", metavar="FLOAT",
                      help="blast program")
    parser.set_defaults(
        prog="eBlastn"
        )
    return parser


def get_program_help():
    """
    Build formatted help on blast program choices.

    @rtype: str
    @return: blast program help
    """
    lines = ["Available blast programs:"]
    lines.extend(['\t%2d: "%s"' % i for i in EProgram.values.iteritems()])
    return "\n".join(lines)


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


def main(argv):

    parser = make_parser()
    opt, args = parser.parse_args(argv)
    try:
        query_file = args[1]
        target_file = args[2]
    except IndexError:
        parser.print_help()
        sys.exit(2)

    # blast options
    options = {}
    for optname in "DbLength", "EvalueThreshold":
        optvalue = getattr(opt, optname)
        if optvalue is not None:
            options[optname] = optvalue
    if opt.prog:
        try:
            options['Program'] = getattr(EProgram, opt.prog)
        except AttributeError:
            print '"%s" is not a valid program' % opt.prog
            print get_program_help()
            sys.exit(2)
            
    if 'Program' not in options or options['Program'] == EProgram.eBlastn:
        options.setdefault('MatchReward', 1)
        #options.setdefault('XDropoff', 30.0)

    try:
        query_f, target_f = open(query_file), open(target_file)
    except IndexError:
        print "Usage: %s QUERY_FILE TARGET_FILE"
        sys.exit(2)

    seq_factory = blast_sseq_loc_from_fasta()
        
    query_str, target_str = query_f.read(), target_f.read()
    query_f.close()
    target_f.close()
    
    query_seq = seq_factory.make(query_str, STRAND, 0, 0, False)
    target_seq = seq_factory.make(target_str, STRAND, 0, 0, False)
    blast_engine = CBl2Seq(SSeqLoc(), SSeqLoc(), options['Program'])
    blast_engine.SetQuery(query_seq)
    os = blast_engine.GetOptions()
            
    for k in options:
        try:
            os[k] = options[k]
        except KeyError:
            print "Invalid Key: %s" % k

    blast_engine.SetSubject(target_seq)
    blast_engine.RunWithoutSeqalignGeneration()
    result = blast_engine.GetResults()
    print "\t".join(["query_start", "query_end", "subject_start",
                     "subject_end", "e_value", "bit_score"])
    hit_list = result[0]
    if hit_list is not None:
        hsp_list = hit_list[0]
        for i, hsp in enumerate(hsp_list):
            if options['Program'] == EProgram.eBlastn:
                progname = "blastn"
            else:
                progname = "tblastx"

            (q_start, q_end, s_start, s_end) = get_match_endpoints(
                hsp.query, hsp.subject, query_seq.length, target_seq.length,
                progname)
            
            print "\t".join(map(str, [q_start, q_end, s_start, s_end,
                                     hsp.evalue, hsp.bit_score]))
            

if __name__ == "__main__":
    main(sys.argv)
