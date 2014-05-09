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
Test blaster speed by comparison with NCBI blastall.
"""

import sys, os, time, random, optparse, logging
import ncbi_toolkit, BlastPython


GID_BASE = 1000
STRAND = ncbi_toolkit.strand.both


def make_seq(n):
    return ''.join([random.choice(['A', 'C', 'G', 'T']) for i in xrange(n)])


def make_fasta(seq, id, title):
    return '>gi|%d %s\n%s\n' % (id, title, seq)


def generate_db_fasta(fname, nseq, seq_len_ave, seq_len_sigma):
    fp = open(fname, 'w')
    for i in xrange(nseq):
        l = random.gauss(seq_len_ave, seq_len_sigma)
        l = max(int(l), int(seq_len_ave - 3*seq_len_sigma))
        s = make_seq(l)
        fp.write(make_fasta(s, GID_BASE + i, 'title %d' % i));
    fp.close()


def formatdb(dbname):
    os.system('formatdb -V -o T -p F -i %s' % dbname)


def str_db(dbname):
    fsio = BlastPython.FS_IO(dbname, 100000)
    fs   = BlastPython.fasta_stream_from_stream(fsio)
    of = open(dbname + '.str', 'w')
    for f in fs:
        h, t = f.split('\n')[0:2]
        i = (h.split()[0]).split("|")[1]
        of.write('%s %s\n' % (i, t))
    of.close()


def run_blastall(pname, dbname, qname, e):
    os.system('blastall -p %s -d %s -i %s -e %s -o blastall.out' % (
        pname, dbname, qname, e))


def run_blaster_stream(pname, query, subjects, e):
    kwds = {'Program' : pname}
    b = BlastPython.blaster(query, **kwds)
    counter = 0
    for subject in subjects:
        _, result = b.blast(subject)
        hit_list = result[0]
        if hit_list is not None:
            hsp_list = hit_list[0]
            for hsp in hsp_list:
                if hsp.evalue < e:
                    counter += 1
    return counter


def run_pyblast_fasta(dbname, qname, pname, e):
    tot_count = 0
    q_factory = BlastPython.seq_factory_from_fasta(STRAND)  
    q_fsio = BlastPython.FS_IO(qname, 1000000)
    q_fs = BlastPython.fasta_stream_from_stream(q_fsio)
    q_ss  = BlastPython.blast_seq_stream(q_factory, q_fs)
    db_factory = BlastPython.seq_factory_from_fasta(STRAND)
    db_fsio = BlastPython.FS_IO(dbname, 1000000)
    db_fs = BlastPython.fasta_stream_from_stream(db_fsio)
    db_ss  = BlastPython.blast_seq_stream(db_factory, db_fs)
    subjects = list(db_ss)
    for query in q_ss:
        tot_count += run_blaster_stream(pname, query, subjects, e)
    logging.info('run_pyblast_fasta: %d hits' % tot_count)


def run_pyblast_str(dbname, qname, pname, e):
    tot_count = 0
    q_factory = BlastPython.seq_factory_from_str(STRAND)  
    q_fsio = BlastPython.FS_IO(qname+'.str', 1000000)
    q_ss  = BlastPython.blast_seq_stream(q_factory, q_fsio)
    db_factory = BlastPython.seq_factory_from_str(STRAND)
    db_fsio = BlastPython.FS_IO(dbname+'.str', 1000000)
    db_ss  = BlastPython.blast_seq_stream(db_factory, db_fsio)
    subjects = list(db_ss)
    for i, query in enumerate(q_ss):
        tot_count += run_blaster_stream(pname, query, subjects, e)
    logging.info('run_pyblast_str: %d hits' % tot_count)


def run_pyblast_str_no_blast(dbname, qname, pname, e):
    fsio = BlastPython.FS_IO(dbname + '.str', 1000000)
    factory_str = BlastPython.seq_factory_from_str(STRAND)
    factory_fasta = BlastPython.seq_factory_from_fasta(STRAND)
    fasta = open(qname).read()
    query = factory_fasta.make(fasta)
    ss  = BlastPython.blast_seq_stream(factory_str, fsio)
    for s in ss:
        pass


def run_all_blast(dbname, qname, pname, e):
    formatdb(dbname)
    start = time.time()
    run_blastall(pname, dbname, qname, e)
    return (time.time() - start)


def run_all_pyblast_fasta(dbname, qname, pname, e):
    start = time.time()
    run_pyblast_fasta(dbname, qname, pname, e)
    return (time.time() - start)


def run_all_pyblast_str(dbname, qname, pname, e):
    start = time.time()
    run_pyblast_str(dbname, qname, pname, e)
    return (time.time() - start)


def run_all_pyblast_str_multiple(dbname, qname, pname, e):
    q_fsio = BlastPython.FS_IO(qname + '.str', 1000000)
    db_fsio = BlastPython.FS_IO(dbname + '.str', 1000000)
    factory = BlastPython.seq_factory_from_str(STRAND)
    queries = list(BlastPython.blast_seq_stream(factory, q_fsio))
    subjects = list(BlastPython.blast_seq_stream(factory, db_fsio))
    kwds = {'Program' : pname}
    b = BlastPython.blaster(queries, **kwds)
    start = time.time()
    _, res = b.blast(subjects)
    counter = 0
    for hit_list in res:
        if hit_list is None:
            continue
        for hsp_list in hit_list:
            for hsp in hsp_list:
                if hsp.evalue < e:
                    counter += 1
    logging.info('run_all_pyblast_str_multiple: %d hits' % counter)
    return (time.time() - start)


def make_parser():
    usage = "%prog [OPTIONS] NQUERIES NSUBJECTS SEQ_LEN_AVE SEQ_LEN_STD"
    parser = optparse.OptionParser(usage)
    parser.set_description(__doc__.lstrip())
    parser.add_option("-p", type="choice", dest="program", metavar="STRING",
                      help="blast program to use",
                      choices = ["blastn", "tblastx"])
    parser.add_option("-v", action="store_const", const=logging.INFO,
                      dest="loglevel", help="display verbose info")

    parser.set_defaults(program="blastn", loglevel=logging.WARNING)
    o = parser.get_option("-p")
    o.help += " %r" % o.choices
    return parser


def main(argv):

    parser = make_parser()
    opt, args = parser.parse_args(argv)

    logging.basicConfig(level=opt.loglevel, format='%(message)s')

    try:
        nseqs_q = int(args[1])
        nseqs_db = int(args[2])
        seq_len_ave = int(args[3])
        seq_len_std = int(args[4])
    except IndexError:
        parser.print_help()
        sys.exit(2)
    except ValueError, e:
        print "ERROR: all arguments beyond the first one must be integers"
        sys.exit(str(e))

    dbname = 'tdb'
    qname = 'query.fa'
    e = 1e10

    generate_db_fasta(qname, nseqs_q, seq_len_ave, seq_len_std)
    generate_db_fasta(dbname, nseqs_db, seq_len_ave, seq_len_std)
    str_db(dbname)
    str_db(qname)
    
    dt1 = run_all_blast(dbname, qname, opt.program, e)

    if opt.program == 'tblastx':
        pname = ncbi_toolkit.EProgram.eTblastx
    elif opt.program == 'blastn':
        pname = ncbi_toolkit.EProgram.eBlastn
    else:
        pname = ncbi_toolkit.EProgram.eTblastx

    dt2 = run_all_pyblast_fasta(dbname, qname, pname, e)
    dt3 = run_all_pyblast_str(dbname, qname, pname, e)
    dt4 = run_all_pyblast_str_multiple(dbname, qname, pname, e)
    logging.info("\t".join(["program", "nseqs_q", "nseqs_db", "seq_len_ave",
                            "seq_len_std", "blastall", "fasta", "str",
                            "str_multi"]))
    print "\t".join(map(str, [opt.program, nseqs_q, nseqs_db, seq_len_ave,
                              seq_len_std, dt1, dt2, dt3, dt4]))


if __name__ == "__main__":
    main(sys.argv)
