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
import subprocess
import logging


#TODO
#zintegrowanie z obecnymi skryptami
#sprawdzenie czy dziala
#Zbieranie logow

class KrakenError(BaseException):
    pass


class KrakenRunner():

    """Wrapper for running kraken."""
    def __init__(self,reads_1, reads_2, threads, database_dir):
        self.logger = logging.getLogger()
        self.reads_1 = reads_1
        self.reads_2 = reads_2
        self.threads = threads
        self.database_dir = database_dir
        # make sure kraken is installed
        self.checkForKraken()


    def run_classification(self):
        pass


    def run_report(self):
        pass

    def checkForKraken(self):
        """Check to see if Kraken is on the system before we try to run it."""

        # Assume that a successful kraken -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['kraken', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  [Error] Make sure kraken is on your system path.")
            sys.exit()
