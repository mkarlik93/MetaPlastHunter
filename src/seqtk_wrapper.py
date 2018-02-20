#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Michal Karlicki'
__copyright__ = 'Copyright 2018'
__credits__ = ['Michal Karlicki']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Michal Karlicki'
__email__ = 'michal.karlicki@gmail.com'
__status__ = 'Development'


import os
import sys
from helpers import *

class SeqtkError(BaseException):
    pass


class SeqtkRunner:

    """Wrapper for running seqtk."""
    def __init__(self,settings):
        #make sure kraken is installed
        if settings == False:
            self.checkForseqtk()
            self.path = ""
        else:
            self.path = Settings_loader(mode="seqtk").read_path()["seqtk"]

    def run_seqtk(self,orig_reads,reads,WorkDir,which_pair):
        os.chdir(WorkDir)
        path = self.path
        command = "%sseqtk subseq %s %s > %s" % (path,orig_reads, reads_1, "reads_"+reads_1+str(which_pair)+".fastq")
        os.system(command)


    def checkForseqtk(self):
        """Check to see if Kraken is on the system before we try to run it."""

        # Assume that a successful kraken -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['seqtk', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            print "  [Error] Make sure seqtk is on your system path or set usage to path in settings.txt"
            sys.exit()
