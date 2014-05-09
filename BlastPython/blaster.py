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
import ncbi_toolkit


class blaster(object):

    def __init__(self, query, **kw):
        # program type must be passed to the constructor to ensure correct
        # default settings for all other options
        program = kw.pop('Program', ncbi_toolkit.EProgram.eBlastn)
        self.query = query
        self.blast_engine = ncbi_toolkit.CBl2Seq(
            ncbi_toolkit.SSeqLoc(),
            ncbi_toolkit.SSeqLoc(),
            program
            )
        self.blast_engine.SetQuery(query)
        self.query_already_setup = False
        opts = self.blast_engine.SetOptions()
        for k in kw:
            opts[k] = kw[k]

    def get_options(self):
        return self.blast_engine.GetOptions()

    def set_options(self):
        return self.blast_engine.SetOptions()

    def blast(self, subject):
        """
        Blast on subject and return results.
        """
        self.blast_engine.SetSubject(subject)
        self.blast_engine.SetupSearch(self.query_already_setup)
        self.query_already_setup = True
        self.blast_engine.ScanDB()
        return subject, self.blast_engine.GetResults()
