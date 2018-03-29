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

import glob
import logging
import os
import sys
from settings import *

#TODO

# Description
#

logger = logging.getLogger('src.krakenize')
logging.basicConfig(level=logging.INFO)

class KrakenError(BaseException):
    pass

class KrakenRunner:

    """ Params




    Wrapper for running kraken."""

    def __init__(self,threads,settings):
        self.threads = threads
        self.settings = settings

        if settings == None:

            self.checkForKraken()
            self.path = ""
        else:

            self.path = Settings_loader(mode="kraken",path=self.settings).read_path()["kraken"]
            self.db = Settings_loader(mode="kraken",path=self.settings).read_database()["kraken_db"]

    def run_classification(self,reads_1, reads_2,name):
        path = self.path
        db = self.db
        threads = self.threads
        command = "%skraken -t %s --db %s --paired %s %s --out-fmt paired --fastq-output --classified-out %s_classif > %s_kraken_out" % (path,str(self.threads),db,reads_1,reads_2,name,name)
        os.system(command)

#    def run_report(self,name):
#        path = self.path
#        db = self.db
#        command_report = "%sscripts/kraken-report --db %s %s_kraken_out > kraken_report.txt" % (path, db,name)
#        os.system(command_report)

    def checkForKraken(self):
        """Check to see if Kraken is on the system before we try to run it."""

        try:
            subprocess.call(['kraken', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            logger.error("Make sure kraken is on your system path or set usage to path in settings.txt")
            sys.exit()

class Pipeline_kraken:

    def __init__(self, list_sra, station_name,settings,threads):

        self.list_sra = list_sra
        self.station_name = station_name
        self.threads = threads
        self.settings = settings

    def kraken_not_multi(self):
        sra_ids = self.list_sra
        list_sra_ids = sra_ids.split(",")
        starting_dir = os.getcwd()
        kraken = KrakenRunner(self.threads,self.settings)
        for i in list_sra_ids:
            logger.info("reads prelimnary classification by kraken, for "+i)
            read_name_1 = i+"_1.fastq"
            read_name_2 = i+"_2.fastq"
            dir = "%s/%s/" % (self.station_name,i)
            os.chdir(dir)
            if len(glob.glob(read_name_1)) == 1:
                kraken.run_classification(read_name_1, read_name_2,i)
#                kraken.run_report(i)
                os.chdir(starting_dir)
                os.remove("%s/%s/%s" % (self.station_name,i,read_name_1))
                os.remove("%s/%s/%s" % (self.station_name,i,read_name_2))
            else:
                logger.error("There is no fastq files, consider option recalculate ")
                sys.exit()
                os.chdir(starting_dir)

    def run(self):
        self.kraken_not_multi()


if __name__ == "__main__":

    from time import gmtime, strftime
    import sys
    import os
    import argparse


    description = """

Version 1.0

This script was designed to use Kraken for metagenomic sequence classification.

If you have any questions, please do not hesitate to contact me
email address: michal.karlicki@gmail.com
"""

    epilog = """

"""


    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('sra_ids', metavar='sra_ids', type=str)
    parser.add_argument('station_name', metavar='station_name', type=str)
    parser.add_argument('threads', metavar='threads', type=int)
    parser.add_argument('settings', metavar='settings', type=str)


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    start = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    run = Pipeline_kraken(args.sra_ids,args.station_name).run()
    end = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    logger.info("Starting time: "+start)
    logger.info("Ending time: "+end)
