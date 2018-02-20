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
from settings import *

#TODO
#zintegrowanie z obecnymi skryptami
#sprawdzenie czy dziala
#Zbieranie logow

class KrakenError(BaseException):
    pass


class KrakenRunner:

    """Wrapper for running kraken."""

    def __init__(self,threads,settings):
        self.threads = threads
        # make sure kraken is installed
        if settings == False:
            self.checkForKraken()
            self.path = ""
        else:
            self.path = Settings_loader(mode="kraken").read_path()["kraken"]
            self.db = Settings_loader(mode="kraken").read_path()["kraken_db"]

    def run_classification(self,reads_1, reads_2,OutDir):
        path = self.path
        db = self.db
        threads = self.threads
        command = "%skraken -t %s --db %s --paired %s %s --out-fmt paired --fastq-output --classified-out classif > %skraken_out" % (path,str(self.threads),db,reads_1,reads_2,OutDir)
        os.system(command)

    def run_report(self,OutputDir):
        path = self.path
        db = self.db
        command_report = "kraken/kraken-report --db %s kraken_out > %skraken_report.txt" % (path, db,str(self.threads),OutputDir)
        os.system(command_report)

    def checkForKraken(self):
        """Check to see if Kraken is on the system before we try to run it."""

        # Assume that a successful kraken -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['kraken', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            print "  [Error] Make sure kraken is on your system path or set usage to path in settings.txt"
            sys.exit()
