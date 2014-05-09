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

import unittest

import test_sseq
import test_blaster
import test_streams

suite_sseq    = test_sseq.suite()
suite_blaster = test_blaster.suite()
suite_streams = test_streams.suite()


alltests = unittest.TestSuite((suite_sseq, suite_blaster, suite_streams))

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(alltests)
